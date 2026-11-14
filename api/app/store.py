"""
Загрузка результатов аналитики (analyze_complexes.py) из папки data и
подготовка их в форме, удобной фронтенду.

Источники (DATA_DIR, по умолчанию /data = <repo>/../data):
  analytics/complexes.geojson  — полигоны + метрики + вердикт
  analytics/timeseries.csv     — ряд GCOV по датам
  analytics/gunw_pairs.csv     — когерентность/смещение по парам GUNW
  analytics/city_stats.json    — даты, пороги, шум города
  static/city/layers.json      — общегородские растровые слои
Каталог ЖК (CATALOG_PATH, опционально) — названия/адреса; подмешивается по id,
чтобы обновлять подписи без перезапуска аналитики.

Файлы перечитываются автоматически, когда меняется время их изменения.
"""

from __future__ import annotations

import csv
import json
import math
import threading
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from pathlib import Path

FILES = ("analytics/complexes.geojson", "analytics/timeseries.csv", "analytics/gunw_pairs.csv",
         "analytics/city_stats.json", "static/city/layers.json")

# Вердикт аналитики -> статус для интерфейса.
STATUS = {
    "activity_detected": "active",
    "possible_activity": "possible",
    "no_changes": "no_changes",
    "insufficient_resolution": "insufficient",
}
STATUS_LABEL_RU = {
    "active": "Стройка идёт",
    "possible": "Возможна активность",
    "no_changes": "Изменений не видно",
    "insufficient": "Объект слишком мал для радара",
    "pending": "Нет данных анализа",
}
RISK_LABEL_RU = {"stable": "Грунт стабилен", "moderate_shift": "Умеренное смещение грунта",
                 "high_risk": "Значительное смещение грунта"}
TRIGGER_RU = {
    "gcov_sci": "яркость радара на участке выросла относительно района",
    "gcov_avi": "сигнал сильно скачет от снимка к снимку",
    "gunw_coherence": "когерентность (80 м) на участке ниже, чем вокруг",
    "gslc_coherence": "детальная когерентность (10×5 м) на участке ниже, чем вокруг",
    "gslc_sci": "детальная яркость (10×5 м) изменилась сильнее, чем у 95% кварталов",
}


def _num(v):
    """Строка из CSV -> int / float / None (пусто, nan, inf) / исходная строка."""
    if v is None or v == "":
        return None
    s = str(v).strip()
    if s.lstrip("-").isdigit():
        return int(s)
    try:
        f = float(s)
    except ValueError:
        return v
    return f if math.isfinite(f) else None


def _r(v, nd=3):
    return None if v is None else round(float(v), nd)


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return [{k: _num(v) if k not in ("complex_id", "obs_date", "date1", "date2", "granule_id",
                                          "crop_png_url") else v for k, v in r.items()}
                for r in csv.DictReader(fh)]


def _read_json(path: Path, default):
    if not path.exists():
        return default
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _centroid(geom: dict) -> list[float]:
    ring = geom["coordinates"][0] if geom["type"] == "Polygon" else geom["coordinates"][0][0]
    pts = ring[:-1] if ring[0] == ring[-1] else ring
    return [round(sum(p[0] for p in pts) / len(pts), 6), round(sum(p[1] for p in pts) / len(pts), 6)]


def _trend_per_month(dates: list[str], values: list) -> float | None:
    """Наклон линейной регрессии значения по времени, единиц в месяц (30.44 дня)."""
    pts = [(date.fromisoformat(d).toordinal(), v) for d, v in zip(dates, values) if v is not None]
    if len(pts) < 3:
        return None
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx == 0:
        return None
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx * 30.44


def _progress(values: list, threshold: float) -> list[float]:
    """Индекс «прогресса» 0..1 по датам для анимации на карте (не высота дома!).

    Накопленный рост очищенного сигнала относительно первой даты, нормированный
    на максимум; если рост меньше порога шума — нули (объект «стоит»).
    """
    vals = [v for v in values if v is not None]
    if not vals:
        return []
    base = vals[0]
    run, out = 0.0, []
    for v in values:
        if v is not None:
            run = max(run, v - base)
        out.append(run)
    top = max(out)
    if top < threshold:
        return [0.0] * len(out)
    return [round(x / top, 3) for x in out]


@dataclass
class Snapshot:
    complexes: dict = field(default_factory=dict)   # id -> {"properties", "geometry"}
    order: list = field(default_factory=list)
    series: dict = field(default_factory=dict)      # id -> list[row]
    pairs: dict = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    layers: list = field(default_factory=list)
    updated_at: str | None = None
    loaded: bool = False


class Store:
    def __init__(self, data_dir: Path, catalog_path: Path | None = None):
        self.data_dir = Path(data_dir)
        self.catalog_path = Path(catalog_path) if catalog_path else None
        self._sig = None
        self._snap = Snapshot()
        self._lock = threading.Lock()

    def _signature(self):
        paths = [self.data_dir / f for f in FILES]
        if self.catalog_path:
            paths.append(self.catalog_path)
        return tuple(p.stat().st_mtime if p.exists() else 0 for p in paths)

    def get(self) -> Snapshot:
        sig = self._signature()
        if sig != self._sig:
            with self._lock:
                if sig != self._sig:
                    self._snap = self._load()
                    self._sig = sig
        return self._snap

    # ------------------------------------------------------------------ #
    def _load(self) -> Snapshot:
        d = self.data_dir
        snap = Snapshot()
        fc = _read_json(d / "analytics/complexes.geojson", None)
        catalog = {}
        if self.catalog_path and self.catalog_path.exists():
            for f in _read_json(self.catalog_path, {}).get("features", []):
                pid = str((f.get("properties") or {}).get("id", ""))
                if pid:
                    catalog[pid] = f
        snap.stats = _read_json(d / "analytics/city_stats.json", {})
        snap.layers = _read_json(d / "static/city/layers.json", [])
        for r in _read_csv(d / "analytics/timeseries.csv"):
            snap.series.setdefault(str(r["complex_id"]), []).append(r)
        for r in _read_csv(d / "analytics/gunw_pairs.csv"):
            snap.pairs.setdefault(str(r["complex_id"]), []).append(r)
        for rows in snap.series.values():
            rows.sort(key=lambda r: r["obs_date"])
        for rows in snap.pairs.values():
            rows.sort(key=lambda r: r["date1"])

        feats = (fc or {}).get("features", [])
        seen = set()
        for f in feats:
            p = dict(f.get("properties") or {})
            cid = str(p.get("id"))
            if cid in catalog:  # подписи из каталога важнее (их правят чаще)
                p.update({k: v for k, v in catalog[cid]["properties"].items() if v not in (None, "")})
            snap.complexes[cid] = {"properties": p, "geometry": f["geometry"]}
            snap.order.append(cid)
            seen.add(cid)
        # ЖК из каталога, для которых анализа ещё нет, — показываем со статусом pending.
        for cid, f in catalog.items():
            if cid not in seen:
                snap.complexes[cid] = {"properties": dict(f["properties"]), "geometry": f["geometry"]}
                snap.order.append(cid)

        p = d / "analytics/complexes.geojson"
        if p.exists():
            snap.updated_at = datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")
        snap.loaded = fc is not None
        return snap


# --------------------------------------------------------------------------- #
# Представления для API
# --------------------------------------------------------------------------- #
def thresholds(snap: Snapshot) -> dict:
    return snap.stats.get("thresholds") or {}


def series_columns(snap: Snapshot, cid: str, base_url: str = "") -> dict:
    rows = snap.series.get(cid, [])
    col = lambda k, nd=2: [_r(r.get(k), nd) for r in rows]
    return {
        "dates": [r["obs_date"] for r in rows],
        "site_hh_db": col("raw_hh_db"),
        "background_hh_db": col("bg_hh_db"),
        "normalized_hh_db": col("norm_hh_db"),
        "site_hv_db": col("raw_hv_db"),
        "background_hv_db": col("bg_hv_db"),
        "normalized_hv_db": col("norm_hv_db"),
        "site_std_hh_db": col("std_hh_db"),
        "crop_png_url": [_abs(r.get("crop_png_url"), base_url) for r in rows],
    }


def _abs(url, base_url: str):
    if not url:
        return None
    return url if url.startswith("http") or not base_url else base_url.rstrip("/") + url


def trigger_ru(t: str, p: dict) -> str:
    """Фраза о признаке с учётом знака: яркость могла вырасти (каркас, бетон) или упасть
    (снос старых построек, расчистка, котлован)."""
    v = {"gcov_sci": p.get("sci_db"), "gslc_sci": p.get("gslc_sci_db")}.get(t)
    if v is not None:
        res = "10×5 м" if t == "gslc_sci" else "20 м"
        if v >= 0:
            return f"яркость радара ({res}) на участке выросла на {v:.1f} дБ относительно района — признак возведения конструкций"
        return (f"яркость радара ({res}) на участке упала на {abs(v):.1f} дБ относительно района — "
                "признак сноса старых построек, расчистки или котлована")
    return TRIGGER_RU.get(t, t)


def reasons_ru(p: dict, status: str) -> list[str]:
    out = []
    trig = [t for t in str(p.get("triggers") or "").split(",") if t]
    elev = [t for t in str(p.get("elevated") or "").split(",") if t]
    pct = p.get("activity_percentile")
    if pct is not None and status in ("active", "possible", "no_changes"):
        out.append(f"совокупность признаков сильнее, чем у {pct:.0f}% случайных участков города того же размера")
    if status == "active":
        out += [trigger_ru(t, p) for t in (trig or elev)]
    elif status == "possible":
        out.append("ни один признак отдельно не превысил порог, но несколько повышены одновременно: "
                   + "; ".join(trigger_ru(t, p) for t in elev))
        out.append("стоит перепроверить на следующих снимках или на месте")
    elif status == "no_changes":
        out.append("признаков активности сверх обычного для города не найдено")
    elif status == "insufficient":
        out.append(f"площадь {p.get('area_sqm', 0):.0f} м² — меньше 4 пикселей радара 20×20 м")
    if p.get("gslc_coh_site") is not None:
        out.append(f"детальная когерентность: {p['gslc_coh_site']:.2f} на участке против "
                   f"{p['gslc_coh_ring']:.2f} вокруг ({p.get('gslc_dates', '')})")
    if p.get("ground_shift_mm") is not None:
        out.append(f"смещение грунта вокруг: {p['ground_shift_mm']:+.1f} ± {p.get('ground_shift_sigma_mm') or 0:.1f} мм "
                   f"— {RISK_LABEL_RU.get(p.get('ground_risk_level'), '')}".rstrip(" —"))
    return out


def summary(snap: Snapshot, cid: str, base_url: str = "") -> dict:
    """Свойства ЖК для карты и списка (лёгкие, без полного ряда)."""
    c = snap.complexes[cid]
    p = c["properties"]
    status = STATUS.get(p.get("sar_verdict"), "pending")
    cols = series_columns(snap, cid)
    noise = (snap.stats.get("noise") or {}).get("avi_db", {}).get("p50", 0.3)
    use_gslc = p.get("gslc_coh_site") is not None
    trig = [t for t in str(p.get("triggers") or "").split(",") if t]
    return {
        "id": cid,
        "name": p.get("name") or f"ЖК {cid}",
        "district": p.get("district"),
        "address": p.get("address"),
        "description": p.get("description"),
        "registry_status": p.get("registry_status"),
        "reported_status": p.get("reported_status"),
        "status": status,
        "status_label": STATUS_LABEL_RU[status],
        "verdict": p.get("sar_verdict"),
        "confidence": p.get("confidence_score"),
        "activity_percentile": _r(p.get("activity_percentile"), 1),
        "triggers": trig,
        "elevated": [t for t in str(p.get("elevated") or "").split(",") if t],
        "center": [_r(p["center_lon"], 6), _r(p["center_lat"], 6)] if p.get("center_lon") is not None
        else _centroid(c["geometry"]),
        "area_sqm": _r(p.get("area_sqm"), 0),
        "pixels": {"gcov_20m": _r(p.get("nisar_pixels"), 1), "gslc_10x5m": p.get("gslc_pixels")},
        "signal": {
            "trend_db_per_month": _r(_trend_per_month(cols["dates"], cols["normalized_hh_db"]), 2),
            "sci_db": _r(p.get("sci_db"), 2),
            "avi_db": _r(p.get("avi_db"), 2),
            "sparkline_db": cols["normalized_hh_db"],
            "progress": _progress(cols["normalized_hh_db"], 2 * noise),
        },
        "coherence": {
            "site": _r(p.get("gslc_coh_site") if use_gslc else p.get("coh_site"), 3),
            "background": _r(p.get("gslc_coh_ring") if use_gslc else p.get("coh_ring"), 3),
            "source": "gslc_10x5m" if use_gslc else ("gunw_80m" if p.get("coh_site") is not None else None),
        },
        "ground": {
            "shift_mm": _r(p.get("ground_shift_mm"), 1),
            "sigma_mm": _r(p.get("ground_shift_sigma_mm"), 1),
            "z": _r(p.get("ground_shift_z"), 2),
            "risk": p.get("ground_risk_level"),
            "risk_label": RISK_LABEL_RU.get(p.get("ground_risk_level")),
        },
        "images": {
            "diff_png": _abs(p.get("diff_heatmap_url"), base_url),
            "gslc_coherence_png": _abs(p.get("gslc_coh_png_url"), base_url),
        },
    }


def detail(snap: Snapshot, cid: str, base_url: str = "") -> dict:
    c = snap.complexes[cid]
    p = c["properties"]
    s = summary(snap, cid, base_url)
    s["reasons"] = reasons_ru(p, s["status"])
    s["timeseries"] = series_columns(snap, cid, base_url)
    s["gunw_pairs"] = [{"date1": r["date1"], "date2": r["date2"], "coh_site": _r(r.get("coh_site")),
                        "coh_background": _r(r.get("coh_ring")), "shift_mm": _r(r.get("shift_mm"), 2)}
                       for r in snap.pairs.get(cid, [])]
    s["gslc"] = None if p.get("gslc_coh_site") is None else {
        "dates": p.get("gslc_dates"), "pixels": p.get("gslc_pixels"),
        "coh_site": _r(p.get("gslc_coh_site")), "coh_background": _r(p.get("gslc_coh_ring")),
        "coh_drop": _r(p.get("gslc_coh_drop")), "delta_db": _r(p.get("gslc_sci_db"), 2),
    }
    s["thresholds"] = thresholds(snap)
    # Всё, что посчитала аналитика, без переименований — для отладки и будущих полей.
    s["raw"] = {k: v for k, v in p.items() if k not in ("name", "district", "address", "description", "note")}
    return s


def city_layers(snap: Snapshot, base_url: str = "") -> list[dict]:
    out = []
    for l in snap.layers:
        w, s_, e, n = l["bbox_west"], l["bbox_south"], l["bbox_east"], l["bbox_north"]
        out.append({
            "id": l["layer_type"],
            "title": {"gcov_change_heatmap": "Изменения яркости радара",
                      "gunw_ground_shift_map": "Смещение грунта (InSAR)",
                      "gslc_coherence_map": "Когерентность 10×5 м"}.get(l["layer_type"], l["layer_type"]),
            "url": _abs(l["png_url"], base_url),
            "date_from": l.get("date_from"), "date_to": l.get("date_to"),
            "units": l.get("units"),
            "bbox": [w, s_, e, n],
            # Готово для MapLibre: map.addSource(id, {type: 'image', url, coordinates})
            "coordinates": [[w, n], [e, n], [e, s_], [w, s_]],
            "legend": {"min": l.get("vmin"), "max": l.get("vmax"), "colormap": l.get("cmap")},
        })
    return out


def haversine_m(lon1, lat1, lon2, lat2) -> float:
    r = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))
