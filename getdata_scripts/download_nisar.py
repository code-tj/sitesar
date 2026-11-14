"""
Скачивание снимков NISAR по Душанбе (трек 012, восходящий, кадр 021) — команда Apex.

Что качает (по умолчанию):
  * GSLC: 2 файла — 21.07.2026 и 25.10.2026 (~12,9 ГБ каждый, ~26 ГБ всего) — когерентность.
  (GCOV и GUNW уже скачаны; их можно включить в настройке WHAT ниже.)

Запуск:
  pip install asf_search
  python download_nisar.py                  # в <repo>/../data/raw (вне репозитория)
  python download_nisar.py ~/files/raw      # или в свою папку (также $NISAR_RAW_DIR)
Логин/пароль NASA Earthdata спросит сам (или задайте переменные EDL_USER / EDL_PASS).
Обрыв связи — просто запустите ещё раз: докачает с места остановки.
Пока файл не докачан, он называется <имя>.h5.part — скрипты обработки его не трогают.
"""
import getpass
import os
import sys

import asf_search as asf

# ---------- настройки ----------
# куда сохранять: аргумент командной строки > $NISAR_RAW_DIR > <repo>/../data/raw
# (вне git-репозитория; $HONAINAV_DATA меняет корень data, как в h5_scripts)
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DATA = os.environ.get("HONAINAV_DATA", os.path.join(os.path.dirname(_REPO), "data"))
OUT_DIR = os.path.abspath(os.path.expanduser(
    sys.argv[1] if len(sys.argv) > 1 else os.environ.get("NISAR_RAW_DIR", os.path.join(_DATA, "raw"))))
WHAT = ["GSLC"]                      # GCOV/GUNW уже скачаны; при надобности: ["GSLC", "GCOV", "GUNW"]
TRACK_FRAME = "_012_A_021_"          # трек 012, восходящий, кадр 021
DUSHANBE = "POINT(68.78 38.56)"
GSLC_DATES = ("20260721T", "20261025T")  # пара для когерентности
CHUNK = 8 * 1024 * 1024
# --------------------------------


def size_of(r):
    """Ожидаемый размер файла в байтах (у NISAR bytes бывает словарём по файлам)."""
    b = r.properties.get("bytes")
    if isinstance(b, dict):
        name = r.properties["fileName"]
        v = b.get(name) or next(iter(b.values()), {})
        return int(v.get("bytes", 0)) if isinstance(v, dict) else int(v or 0)
    return int(b or 0)


def pick():
    res = asf.search(platform="NISAR", processingLevel=WHAT,
                     intersectsWith=DUSHANBE, start="2026-06-01")
    out = []
    for r in res:
        fid = r.properties["fileID"]
        if TRACK_FRAME not in fid:
            continue
        if "_GSLC_" in fid and not any(d in fid for d in GSLC_DATES):
            continue                      # из GSLC берём только 2 нужные даты
        out.append(r)
    order = {"GSLC": 0, "GCOV": 1, "GUNW": 2}
    out.sort(key=lambda r: (order.get(r.properties["processingLevel"], 9), r.properties["startTime"]))
    return out


def download(r, sess):
    url = r.properties["url"]
    name = r.properties["fileName"]
    path = os.path.join(OUT_DIR, name)
    part = path + ".part"                         # недокачанный файл — под другим именем
    total = size_of(r)
    if os.path.exists(path) and (not total or os.path.getsize(path) >= total):
        print(f"  уже есть: {name}")
        return
    have = os.path.getsize(part) if os.path.exists(part) else 0
    if total and have >= total:
        os.replace(part, path)
        print(f"  уже есть: {name}")
        return
    headers = {"Range": f"bytes={have}-"} if have else {}
    with sess.get(url, headers=headers, stream=True, timeout=120) as resp:
        if resp.status_code == 200 and have:      # сервер не умеет докачку — начинаем заново
            have = 0
        resp.raise_for_status()
        with open(part, "ab" if have else "wb") as f:
            done = have
            for chunk in resp.iter_content(CHUNK):
                f.write(chunk)
                done += len(chunk)
                if total:
                    sys.stdout.write(f"\r  {name[:60]}…  {done/1e9:5.2f} / {total/1e9:5.2f} ГБ")
                    sys.stdout.flush()
    if total and os.path.getsize(part) < total:
        raise IOError(f"соединение оборвалось на {os.path.getsize(part)/1e9:.2f} ГБ")
    os.replace(part, path)                        # переименовываем только целый файл
    print()


def main():
    # Под nohup вывод идёт в файл и без этого копится в буфере — лог выглядит пустым.
    sys.stdout.reconfigure(line_buffering=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.environ.get("EDL_USER") and not sys.stdin.isatty():
        sys.exit("Задайте EDL_USER и EDL_PASS: в фоне (nohup) логин спросить нельзя.")
    user = os.environ.get("EDL_USER") or input("Earthdata логин: ")
    pwd = os.environ.get("EDL_PASS") or getpass.getpass("Earthdata пароль: ")
    sess = asf.ASFSession().auth_with_creds(user, pwd)

    files = pick()
    gb = sum(size_of(r) for r in files) / 1e9
    print(f"Найдено файлов: {len(files)}, примерно {gb:.1f} ГБ")
    for r in files:
        print(f"  {r.properties['processingLevel']}  {r.properties['startTime'][:10]}  {r.properties['fileName']}")
    if not files:
        sys.exit("Ничего не найдено — проверьте интернет и логин.")

    for i, r in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {r.properties['processingLevel']} {r.properties['startTime'][:10]}")
        for attempt in range(3):
            try:
                download(r, sess)
                break
            except Exception as e:                # обрыв — пробуем ещё, докачка с места остановки
                print(f"\n  ошибка: {e}; повтор {attempt + 1}/3")
        else:
            print("  не удалось — запустите скрипт ещё раз позже")
    print("Готово. Файлы в папке:", os.path.abspath(OUT_DIR))


if __name__ == "__main__":
    main()
