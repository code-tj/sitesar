#!/usr/bin/env python3
"""
Вырезка территории Душанбе из HDF5-продуктов NASA-ISRO NISAR (L2 GCOV, L2 GUNW)
в лёгкие GeoTIFF для дальнейшей аналитики (SiteSAR, NASA Space Apps 2026).

Что делает:
  1. Находит в папке файлы NISAR_*.h5 (рекурсивно), разбирает имя гранулы:
     продукт, цикл, трек, направление орбиты (A/D), кадр, даты съёмки.
     Можно отфильтровать: --product GCOV --track 12 --direction A --frame 21.
  2. Внутри HDF5 находит сетки — группы с xCoordinates / yCoordinates
     (например /science/LSAR/GCOV/grids/frequencyA/ или
     /science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/).
  3. Берёт EPSG из датасета `projection` (обычно 32642 — UTM 42N), переводит
     bbox Душанбе из градусов в метры и находит окно пикселей.
  4. Читает с диска ТОЛЬКО это окно (h5py читает нужные чанки, а не 2 ГБ файла)
     и пишет многоканальный GeoTIFF (каналы = слои, напр. HHHH и HVHV),
     значения как в исходнике (gamma0 в линейных единицах, фаза в радианах).
  5. Пишет manifest.csv со списком всех результатов (дата, трек, слои, границы)
     — его удобно читать следующим скриптом (анализ по полигонам ЖК).

Результат по умолчанию — ВНЕ репозитория, в папке data рядом с ним
(<repo>/../data/raw_nisar, переопределяется -o или переменной HONAINAV_DATA):
  ../data/raw_nisar/dushanbe_gcov_20260709.tif         (HHHH, HVHV)
  ../data/raw_nisar/dushanbe_gunw_20260721_20260814.tif (unwrappedPhase, coherenceMagnitude)

Примеры:
  python3 crop_dushanbe.py /data/nisar
  python3 crop_dushanbe.py /data/nisar -o /mnt/big/raw_nisar --product GCOV --track 12 --direction A --frame 21
  python3 crop_dushanbe.py /data/nisar --layers HHHH HVHV numberOfLooks
  python3 crop_dushanbe.py /data/nisar --info          # показать структуру и окно, ничего не писать
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import rasterio
from pyproj import CRS, Transformer
from rasterio.transform import from_origin

# bbox Душанбе (lon_min, lat_min, lon_max, lat_max), WGS84 — как в описании проекта.
DUSHANBE_BBOX = (68.68, 38.48, 68.88, 38.62)
DEFAULT_EPSG = 32642  # UTM Zone 42N

# Объёмные данные храним вне git-репозитория: <repo>/../data (или $HONAINAV_DATA).
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("HONAINAV_DATA", REPO_ROOT.parent / "data")).resolve()
DEFAULT_OUTPUT = DATA_ROOT / "raw_nisar"

# Слои по умолчанию для каждого продукта (имя датасета или его конец пути).
DEFAULT_LAYERS = {
    "GCOV": ["HHHH", "HVHV"],
    "GUNW": ["unwrappedPhase", "coherenceMagnitude", "connectedComponents", "ionospherePhaseScreen"],
    "GSLC": ["HH", "HV"],
}
# Из каких групп брать слои по умолчанию (в GUNW есть ещё wrappedInterferogram,
# pixelOffsets с одноимёнными слоями — их берём только с флагом --all-groups).
DEFAULT_GROUP_HINT = {
    "GUNW": "unwrappedInterferogram",
}
DEFAULT_FREQUENCY = "A"

# Продукты, где нужна фаза: комплексные слои пишутся как complex64, а не |z|
# (из пары GSLC потом считается когерентность).
COMPLEX_PRODUCTS = {"GSLC", "RSLC"}

# NISAR_L2_PR_GCOV_025_012_A_021_2005_DHDH_M_20260709T003215_20260709T003250_...
GRANULE_RE = re.compile(
    r"NISAR_(?P<level>L\d)_(?P<proc>[A-Z]{2})_(?P<product>[A-Z0-9]{4})_"
    r"(?P<cycle>\d{3})_(?P<track>\d{3})_(?P<direction>[AD])_(?P<frame>\d{3})_")
DATE_RE = re.compile(r"(20\d{6})T\d{6}")

COORD_NAMES = ("xCoordinates", "yCoordinates")

log = logging.getLogger("crop_dushanbe")


# --------------------------------------------------------------------------- #
@dataclass
class Granule:
    path: Path
    product: str = "UNKNOWN"
    cycle: int | None = None
    track: int | None = None
    direction: str | None = None
    frame: int | None = None
    dates: list[str] = field(default_factory=list)

    @classmethod
    def from_path(cls, path: Path) -> "Granule":
        g = cls(path)
        m = GRANULE_RE.search(path.name)
        if m:
            g.product = m["product"]
            g.cycle, g.track, g.frame = int(m["cycle"]), int(m["track"]), int(m["frame"])
            g.direction = m["direction"]
        # уникальные даты в порядке появления (GCOV: 1 дата, GUNW: ref + sec)
        for d in DATE_RE.findall(path.name):
            if d not in g.dates:
                g.dates.append(d)
        return g

    def detect_product(self, f: h5py.File):
        """Если имя файла нестандартное — определить продукт по структуре."""
        if self.product != "UNKNOWN":
            return
        for prod in ("GCOV", "GUNW", "GSLC", "GOFF"):
            if f"/science/LSAR/{prod}" in f or f"/science/SSAR/{prod}" in f:
                self.product = prod
                return


@dataclass
class Grid:
    group: str            # путь группы с xCoordinates/yCoordinates
    x: np.ndarray         # центры пикселей, м
    y: np.ndarray
    epsg: int
    layers: dict          # относительное имя -> полный путь датасета

    @property
    def dx(self) -> float:
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        return float(self.y[1] - self.y[0])  # обычно < 0 (север вверху)


# --------------------------------------------------------------------------- #
def _read_epsg(f: h5py.File, group: str, fallback: int) -> int:
    """projection — скаляр с EPSG-кодом и/или атрибутом epsg_code."""
    g = f[group]
    while True:
        if "projection" in g:
            p = g["projection"]
            for key in ("epsg_code", "epsg"):
                if key in p.attrs:
                    return int(np.asarray(p.attrs[key]).ravel()[0])
            try:
                val = int(np.asarray(p[()]).ravel()[0])
                if 1000 < val < 100000:
                    return val
            except (TypeError, ValueError):
                pass
            if "spatial_ref" in p.attrs:
                wkt = p.attrs["spatial_ref"]
                wkt = wkt.decode() if isinstance(wkt, bytes) else str(wkt)
                epsg = CRS.from_wkt(wkt).to_epsg()
                if epsg:
                    return epsg
        if g.name == "/":
            break
        g = g.parent
    log.warning("%s: EPSG не найден в файле, используем %d", group, fallback)
    return fallback


def find_grids(f: h5py.File, frequency: str | None, fallback_epsg: int) -> list[Grid]:
    """Найти все группы, где лежат xCoordinates и yCoordinates."""
    groups = []

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and all(c in obj for c in COORD_NAMES):
            groups.append("/" + name)

    f.visititems(visit)
    # Служебные сетки калибровки/геометрии (…/metadata/…, шаг ~1 км) — не данные.
    groups = [g for g in groups if "/metadata/" not in g]
    if frequency:
        groups = [g for g in groups if f"frequency{frequency}" in g] or groups

    grids = []
    for gpath in groups:
        x = f[gpath]["xCoordinates"][()].astype(np.float64)
        y = f[gpath]["yCoordinates"][()].astype(np.float64)
        if x.ndim != 1 or y.ndim != 1 or x.size < 2 or y.size < 2:
            continue
        shape = (y.size, x.size)
        layers = {}

        def collect(name, obj, gpath=gpath, shape=shape):
            if isinstance(obj, h5py.Dataset) and obj.shape == shape:
                layers[name] = f"{gpath}/{name}"

        f[gpath].visititems(collect)
        # Не залезать в подгруппы, у которых своя сетка.
        nested = [g[len(gpath) + 1:] for g in groups if g.startswith(gpath + "/")]
        layers = {k: v for k, v in layers.items() if not any(k.startswith(n + "/") for n in nested)}
        if layers:
            grids.append(Grid(gpath, x, y, _read_epsg(f, gpath, fallback_epsg), layers))
    return grids


def select_layers(grid: Grid, wanted: list[str] | None, product: str) -> list[str]:
    """Сопоставить запрошенные имена (HHHH, unwrappedPhase, HH/unwrappedPhase) со слоями сетки."""
    if wanted and wanted != ["all"]:
        names = wanted
    elif wanted == ["all"]:
        names = None
    else:
        names = DEFAULT_LAYERS.get(product)

    real = {k: v for k, v in grid.layers.items()}
    if names is None:
        return sorted(real)
    out = []
    for n in names:
        hits = [k for k in real if k == n or k.endswith("/" + n)]
        out.extend(h for h in sorted(hits) if h not in out)
    return out


def bbox_to_window(grid: Grid, bbox, densify: int = 21):
    """bbox в градусах -> окно (r0, r1, c0, c1) в пикселях сетки."""
    lon_min, lat_min, lon_max, lat_max = bbox
    tr = Transformer.from_crs(4326, grid.epsg, always_xy=True)
    # Плотная обводка рамки: в UTM прямоугольник из градусов слегка искривлён.
    t = np.linspace(0, 1, densify)
    lons = np.concatenate([lon_min + t * (lon_max - lon_min), np.full(densify, lon_max),
                           lon_max - t * (lon_max - lon_min), np.full(densify, lon_min)])
    lats = np.concatenate([np.full(densify, lat_min), lat_min + t * (lat_max - lat_min),
                           np.full(densify, lat_max), lat_max - t * (lat_max - lat_min)])
    xs, ys = tr.transform(lons, lats)
    xmin, xmax, ymin, ymax = np.min(xs), np.max(xs), np.min(ys), np.max(ys)

    cols = np.where((grid.x >= xmin) & (grid.x <= xmax))[0]
    rows = np.where((grid.y >= ymin) & (grid.y <= ymax))[0]
    if cols.size == 0 or rows.size == 0:
        return None
    return int(rows.min()), int(rows.max()) + 1, int(cols.min()), int(cols.max()) + 1


# --------------------------------------------------------------------------- #
def _out_name(gr: Granule, grid: Grid, multi_grid: bool) -> str:
    dates = "_".join(gr.dates) if gr.dates else gr.path.stem
    name = f"dushanbe_{gr.product.lower()}_{dates}"
    if gr.direction == "D":  # нисходящая орбита — отдельный ряд, не путать с A
        name += f"_D{gr.track:03d}"
    if multi_grid:
        # .../frequencyA/unwrappedInterferogram/HH -> unwrappedInterferogram_HH
        rel = re.split(r"/frequency[AB]/?", grid.group)[-1] or grid.group.rsplit("/", 1)[-1]
        name += "_" + rel.strip("/").replace("/", "_")
    return name + ".tif"


def write_geotiff(path: Path, bands: list[np.ndarray], names: list[str], grid: Grid,
                  window, tags: dict, cog: bool):
    r0, r1, c0, c1 = window
    height, width = r1 - r0, c1 - c0
    # Координаты в NISAR — центры пикселей; угол = центр - полпикселя.
    transform = from_origin(grid.x[c0] - grid.dx / 2, grid.y[r0] - grid.dy / 2,
                            abs(grid.dx), abs(grid.dy))
    if grid.dy > 0:  # на всякий случай: сетка «юг вверху» -> переворачиваем
        bands = [b[::-1] for b in bands]
        transform = from_origin(grid.x[c0] - grid.dx / 2, grid.y[r1 - 1] + grid.dy / 2,
                                abs(grid.dx), abs(grid.dy))

    cplx = any(np.iscomplexobj(b) for b in bands)
    dtype = "complex64" if cplx else "float32"
    predictor = 1 if cplx else 3  # предиктор 3 только для вещественных float
    profile = dict(driver="GTiff", height=height, width=width, count=len(bands),
                   dtype=dtype, crs=CRS.from_epsg(grid.epsg).to_wkt(), transform=transform,
                   nodata=None if cplx else np.nan, compress="deflate", predictor=predictor,
                   tiled=width >= 256 and height >= 256, BIGTIFF="IF_SAFER")
    tmp = path.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", **profile) as dst:
        for i, (b, n) in enumerate(zip(bands, names), 1):
            dst.write(b.astype(dtype), i)
            dst.set_band_description(i, n)
        dst.update_tags(**tags)
    if cog:
        from rasterio.shutil import copy as rio_copy
        rio_copy(tmp, path, driver="COG", compress="deflate", predictor=predictor)
        tmp.unlink()
    else:
        os.replace(tmp, path)


def process(gr: Granule, out_dir: Path, bbox, layers_req, frequency, fallback_epsg,
            overwrite: bool, cog: bool, all_groups: bool = False,
            keep_complex: bool = False) -> list[dict]:
    rows = []
    with h5py.File(gr.path, "r") as f:
        gr.detect_product(f)
        grids = find_grids(f, frequency, fallback_epsg)
        if not grids:
            raise RuntimeError("не найдены xCoordinates/yCoordinates")

        plans = []
        for grid in grids:
            names = select_layers(grid, layers_req, gr.product)
            if names:
                plans.append((grid, names))
        hint = DEFAULT_GROUP_HINT.get(gr.product)
        if hint and not all_groups:
            plans = [p for p in plans if hint in p[0].group] or plans
        if not plans:
            raise RuntimeError("нет запрошенных слоёв; доступно: "
                               + ", ".join(sorted({k for g in grids for k in g.layers})))

        for grid, names in plans:
            window = bbox_to_window(grid, bbox)
            if window is None:
                log.warning("%s: bbox вне покрытия сетки %s", gr.path.name, grid.group)
                continue
            out = out_dir / _out_name(gr, grid, len(plans) > 1)
            r0, r1, c0, c1 = window
            if out.exists() and not overwrite:
                log.info("SKIP  %s: уже есть", out.name)
            else:
                bands, band_names = [], []
                for n in names:
                    ds = f[grid.layers[n]]
                    data = ds[r0:r1, c0:c1]  # читаются только нужные чанки
                    # Комплексные числа, которые h5py не распознал сам (напр. пары float16 r/i).
                    if data.dtype.names and len(data.dtype.names) == 2:
                        re_, im_ = data.dtype.names
                        data = data[re_].astype(np.float32) + 1j * data[im_].astype(np.float32)
                    if np.iscomplexobj(data):
                        if gr.product in COMPLEX_PRODUCTS or keep_complex:
                            data = data.astype(np.complex64)
                        else:
                            log.warning("%s: слой %s комплексный — пишем модуль |z|", gr.path.name, n)
                            data = np.abs(data).astype(np.float32)
                    else:
                        data = data.astype(np.float32)
                        fv = ds.attrs.get("_FillValue")
                        if fv is not None:
                            data[data == np.float32(np.asarray(fv).ravel()[0])] = np.nan
                    bands.append(data)
                    band_names.append(n)
                # Флаги обработки фазы GSLC (…FlatteningApplied) — нужны для интерферометрии.
                flags = {}

                def _flag(name, obj):
                    if isinstance(obj, h5py.Dataset) and name.endswith("FlatteningApplied"):
                        v = obj[()]
                        v = v.decode() if isinstance(v, bytes) else v
                        flags[name.rsplit("/", 1)[-1]] = str(v.item() if hasattr(v, "item") else v)

                f.visititems(_flag)
                tags = dict(**flags, granule=gr.path.name, product=gr.product, dates=",".join(gr.dates),
                            cycle=str(gr.cycle), track=str(gr.track), direction=str(gr.direction),
                            frame=str(gr.frame), source_group=grid.group, layers=",".join(band_names),
                            bbox_wgs84=",".join(map(str, bbox)),
                            window_rows=f"{r0}:{r1}", window_cols=f"{c0}:{c1}")
                out_dir.mkdir(parents=True, exist_ok=True)
                write_geotiff(out, bands, band_names, grid, window, tags, cog)
                log.info("OK    %s -> %s [%dx%d px, %s]", gr.path.name, out.name,
                         r1 - r0, c1 - c0, ",".join(band_names))
            rows.append(dict(
                tif=out.name, granule=gr.path.name, product=gr.product,
                date=gr.dates[0] if gr.dates else "",
                date2=gr.dates[1] if len(gr.dates) > 1 else "",
                cycle=gr.cycle, track=gr.track, direction=gr.direction, frame=gr.frame,
                epsg=grid.epsg, layers=",".join(names),
                width=c1 - c0, height=r1 - r0, pixel_m=f"{abs(grid.dx):g}x{abs(grid.dy):g}"))
    return rows


def info(gr: Granule, bbox, frequency, fallback_epsg):
    print(f"=== {gr.path.name}")
    print(f"    product={gr.product} cycle={gr.cycle} track={gr.track} dir={gr.direction} "
          f"frame={gr.frame} dates={gr.dates}")
    with h5py.File(gr.path, "r") as f:
        gr.detect_product(f)
        grids = find_grids(f, None, fallback_epsg)
        if not grids:
            print("    !! сетки (xCoordinates/yCoordinates) не найдены. Верхние группы:")
            f.visit(lambda n: print("      /" + n) if n.count("/") < 5 else None)
        for g in grids:
            w = bbox_to_window(g, bbox)
            print(f"    grid {g.group}  EPSG:{g.epsg}  {g.y.size}x{g.x.size}  "
                  f"pixel={abs(g.dx):g}x{abs(g.dy):g} m  window={w}")
            for k, v in sorted(g.layers.items()):
                print(f"      {k:40s} {f[v].dtype}")
            hint = DEFAULT_GROUP_HINT.get(gr.product)
            default = select_layers(g, None, gr.product) if not hint or hint in g.group else []
            print(f"      -> по умолчанию: {default}")


# --------------------------------------------------------------------------- #
LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"


def _init_logging(level: int):
    """Настройка логов; вызывается и в дочерних процессах пула."""
    logging.basicConfig(level=level, format=LOG_FORMAT, datefmt="%H:%M:%S", force=True)


def _worker(a):
    gr, *rest = a
    try:
        return gr.path.name, process(gr, *rest), None
    except Exception as e:  # noqa: BLE001
        return gr.path.name, [], f"{type(e).__name__}: {e}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Вырезать Душанбе из NISAR L2 HDF5 (GCOV/GUNW) в GeoTIFF.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("inputs", nargs="*", default=[str(DATA_ROOT / "raw")],
                    help="файлы .h5 и/или папки с ними (по умолчанию <repo>/../data/raw, "
                         "куда качает getdata_scripts/download_nisar.py)")
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT), help="папка для GeoTIFF")
    ap.add_argument("--bbox", nargs=4, type=float, default=DUSHANBE_BBOX,
                    metavar=("LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX"))
    ap.add_argument("--product", nargs="+", default=["GCOV", "GUNW", "GSLC"], help="какие продукты брать")
    ap.add_argument("--track", type=int, help="относительный трек, напр. 12")
    ap.add_argument("--frame", type=int, help="кадр, напр. 21")
    ap.add_argument("--direction", choices=["A", "D"], help="A — восходящая, D — нисходящая")
    ap.add_argument("--frequency", default=DEFAULT_FREQUENCY, choices=["A", "B"])
    ap.add_argument("--layers", nargs="+",
                    help="слои (имя или конец пути). 'all' — все. По умолчанию: "
                         + "; ".join(f"{k}: {' '.join(v)}" for k, v in DEFAULT_LAYERS.items()))
    ap.add_argument("--all-groups", action="store_true",
                    help="искать слои во всех группах (напр. и в wrappedInterferogram GUNW), "
                         "а не только в основной")
    ap.add_argument("--keep-complex", action="store_true",
                    help="комплексные слои писать как complex64 (для GSLC/RSLC включено всегда)")
    ap.add_argument("--epsg", type=int, default=DEFAULT_EPSG, help="EPSG, если в файле его нет")
    ap.add_argument("--cog", action="store_true", help="писать Cloud-Optimized GeoTIFF")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("-j", "--jobs", type=int, default=2, help="параллельных процессов")
    ap.add_argument("--info", action="store_true", help="только показать структуру и окно")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    _init_logging(log_level)
    bbox = tuple(args.bbox)

    files = []
    for inp in args.inputs:
        p = Path(inp)
        if p.is_dir():
            # *.h5.part (недокачанные файлы) отсекаются по расширению
            files += [q for q in p.rglob("*") if q.suffix.lower() in (".h5", ".hdf5")]
        elif p.is_file():
            files.append(p)
        else:
            log.warning("Не найдено: %s", inp)
    granules = []
    for fp in sorted(set(files)):
        gr = Granule.from_path(fp)
        if gr.product != "UNKNOWN" and gr.product not in args.product:
            continue
        if args.track is not None and gr.track not in (None, args.track):
            continue
        if args.frame is not None and gr.frame not in (None, args.frame):
            continue
        if args.direction and gr.direction not in (None, args.direction):
            continue
        granules.append(gr)
    if not granules:
        log.error("Подходящих файлов не найдено")
        return 1
    log.info("Файлов: %d, bbox=%s", len(granules), bbox)

    if args.info:
        for gr in granules:
            try:
                info(gr, bbox, args.frequency, args.epsg)
            except Exception as e:  # noqa: BLE001
                print(f"=== {gr.path.name}\n    !! {type(e).__name__}: {e}")
        return 0

    out_dir = Path(args.output)
    tasks = [(gr, out_dir, bbox, args.layers, args.frequency, args.epsg, args.overwrite, args.cog,
              args.all_groups, args.keep_complex)
             for gr in granules]
    rows, fails = [], 0
    if args.jobs > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(args.jobs, initializer=_init_logging, initargs=(log_level,)) as pool:
            results = [fut.result() for fut in as_completed([pool.submit(_worker, t) for t in tasks])]
    else:
        results = [_worker(t) for t in tasks]
    for name, r, err in results:
        if err:
            fails += 1
            log.error("FAIL  %s: %s", name, err)
        rows += r

    if rows:
        out_dir.mkdir(parents=True, exist_ok=True)
        manifest = out_dir / "manifest.csv"
        # Сливаем с прошлым запуском, чтобы manifest описывал всю папку.
        old = {}
        if manifest.exists():
            with open(manifest, newline="") as fh:
                old = {r["tif"]: r for r in csv.DictReader(fh)}
        old.update({r["tif"]: {k: str(v) for k, v in r.items()} for r in rows})
        # Убираем записи о файлах, которых больше нет в папке.
        old = {k: v for k, v in old.items() if (out_dir / Path(k).name).exists()}
        merged = sorted(old.values(), key=lambda r: (r["product"], r["direction"], r["date"], r["tif"]))
        with open(manifest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(merged)
        log.info("manifest: %s (%d записей)", manifest, len(merged))

    log.info("Готово: %d файлов, ошибок %d", len(tasks), fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
