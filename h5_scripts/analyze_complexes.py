#!/usr/bin/env python3
"""
Аналитика строительной активности по полигонам ЖК Душанбе на основе вырезок NISAR
(результат crop_dushanbe.py: ../data/raw_nisar/*.tif + manifest.csv).

Для каждого ЖК:
  1. Полигон переводится в метры (UTM растра). Строится:
       site — полигон, «размазанный» на --layover-m в сторону спутника
              (верхушки строящихся домов из-за layover смещаются к спутнику);
       ring — фоновое кольцо от --ring-inner до --ring-outer м вокруг ЖК
              (минус другие стройплощадки) — «спокойный квартал» для вычитания
              дождя/влажности/калибровки.
  2. По каждой дате GCOV:
       raw_db  = 10·log10(mean(gamma0 в site))
       bg_db   = 10·log10(median(gamma0 в ring))
       norm_db = raw_db − bg_db            (очищенный от погоды сигнал)
     для HH и HV.
  3. Индексы по временному ряду (norm HH):
       SCI = mean(последние k) − mean(первые k)   (структурный рост, дБ)
       AVI = mean |norm(t+1) − norm(t)|           (волатильность, дБ)
  4. По парам GUNW:
       когерентность в site и в ring (падение когерентности = земляные работы);
       смещение грунта в кольце, мм по лучу зрения (LOS), ОТНОСИТЕЛЬНО города:
       фаза − ионосфера − медиана фазы своей связной области, × (−λ/4π).
       Накопленное смещение = сумма по цепочке пар.
  4b. Если вырезаны два GSLC (пиксель 10×5 м) — когерентность и изменение
       яркости между первой и последней датой GSLC, в 8 раз детальнее GCOV
       (окно усреднения --gslc-window-m, по умолчанию 40 м).
  5. Вердикт: insufficient_resolution / activity_detected / possible_activity / no_changes.
     Каждый ЖК сравнивается с --calib-n случайными участками города ТОГО ЖЕ размера и с той же
     маской навала: признаки объединяются методом Фишера (по эмпирическим рангам), и
       activity_detected — совокупность сильнее, чем у --active-pct % (95) таких участков,
       possible_activity — сильнее, чем у --possible-pct % (90).
     Доля ложных срабатываний = 100 − перцентиль, одинакова для больших и маленьких ЖК.
     --threshold-rule возвращает старое правило «любой признак выше порога p95».
     + confidence 0..100 (по числу пикселей, числу дат и силе сигнала).
  6. Шум фона города: те же индексы для случайных площадок по всему городу —
     для калибровки порогов (city_stats.json).

Выход (по умолчанию в <repo>/../data):
  analytics/complexes.geojson   — полигоны + все метрики и вердикт (для БД / фронта)
  analytics/timeseries.csv      — ряд по датам для графика
  analytics/gunw_pairs.csv      — когерентность и смещение по каждой паре GUNW
  analytics/city_stats.json     — шум фона города и использованные пороги
  static/crops/jk_<id>_<date>_hh.png, jk_<id>_diff.png
  static/city/dushanbe_gcov_change.png, dushanbe_insar_shift.png + layers.json (границы WGS84)

Пример:
  python3 analyze_complexes.py geojson/complexes_dushanbe.geojson
  python3 analyze_complexes.py geojson/complexes_dushanbe.geojson --no-png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
import shapely
from shapely import affinity
from shapely.geometry import box as shp_box
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("HONAINAV_DATA", REPO_ROOT.parent / "data")).resolve()

WAVELENGTH_MM = 238.4  # NISAR L-band, 1.257 ГГц
PHASE_TO_MM = -WAVELENGTH_MM / (4 * math.pi)  # ≈ −18.97 мм/рад
# Азимут направления «на спутник» (от севера по часовой). NISAR смотрит вправо:
# восходящая орбита (курс ~350°) смотрит на восток → спутник на западе (~260°);
# нисходящая (курс ~190°) смотрит на запад → спутник на востоке (~100°).
TOWARD_SAT_AZ = {"A": 260.0, "D": 100.0}

log = logging.getLogger("analyze")


# --------------------------------------------------------------------------- #
@dataclass
class Scene:
    path: Path
    product: str
    date: str
    date2: str
    cycle: str
    granule: str
    direction: str
    bands: dict
    transform: rasterio.Affine
    crs: object
    shape: tuple
    tags: dict = field(default_factory=dict)

    @property
    def pixel_area(self) -> float:
        return abs(self.transform.a * self.transform.e)

    @property
    def pixel_m(self) -> float:
        return math.sqrt(self.pixel_area)


@dataclass
class Complex:
    cid: str
    props: dict
    geom_wgs: object
    poly: object = None   # UTM
    site: object = None
    ring: object = None
    metrics: dict = field(default_factory=dict)
    series: list = field(default_factory=list)
    pairs: list = field(default_factory=list)
    calib: dict = field(default_factory=dict, repr=False)


def iso(d: str) -> str:
    return f"{d[:4]}-{d[4:6]}-{d[6:8]}" if d and len(d) == 8 else d


def to_db(x):
    return 10.0 * np.log10(np.maximum(x, 1e-6))


# --------------------------------------------------------------------------- #
def load_scenes(raw_dir: Path, product: str) -> list[Scene]:
    manifest = raw_dir / "manifest.csv"
    if manifest.exists():
        with open(manifest, newline="") as fh:
            rows = [r for r in csv.DictReader(fh) if r["product"] == product]
    else:
        rows = [dict(tif=p.name, product=product, date=p.stem.split("_")[2], date2="",
                     cycle="", granule="", direction="A")
                for p in sorted(raw_dir.glob(f"dushanbe_{product.lower()}_*.tif"))]
    scenes = []
    for r in sorted(rows, key=lambda r: (r["date"], r.get("date2", ""))):
        p = Path(r["tif"])
        p = p if p.is_absolute() else raw_dir / p.name
        if not p.exists():
            log.warning("нет файла %s", p)
            continue
        with rasterio.open(p) as src:
            bands = {}
            for i in range(1, src.count + 1):
                a = src.read(i)  # комплексные (GSLC) оставляем комплексными — нужна фаза
                bands[(src.descriptions[i - 1] or f"b{i}").split("/")[-1]] = \
                    a.astype(np.complex64) if np.iscomplexobj(a) else a.astype(np.float32)
            scenes.append(Scene(p, product, r["date"], r.get("date2", ""), r.get("cycle", ""),
                                r.get("granule", ""), r.get("direction") or "A", bands,
                                src.transform, src.crs, (src.height, src.width), dict(src.tags())))
    return scenes


def dedupe_by_date(scenes: list[Scene]) -> list[Scene]:
    """Одна сцена на дату — самая большая по площади (отсекает служебные мелкие сетки)."""
    best: dict = {}
    for sc in scenes:
        if sc.date not in best or sc.shape[0] * sc.shape[1] > best[sc.date].shape[0] * best[sc.date].shape[1]:
            best[sc.date] = sc
    return [best[d] for d in sorted(best)]


def load_complexes(path: Path) -> list[Complex]:
    fc = json.load(open(path, encoding="utf-8"))
    out = []
    for i, f in enumerate(fc["features"], 1):
        props = f.get("properties") or {}
        cid = str(props.get("id") or props.get("ext_id") or i)
        geom = shape(f["geometry"])
        if geom.geom_type not in ("Polygon", "MultiPolygon"):
            log.warning("ЖК %s: тип %s пропущен", cid, geom.geom_type)
            continue
        out.append(Complex(cid, props, geom))
    return out


def floors_of(props: dict) -> int | None:
    """Этажность: поле max_floors, иначе максимум «NN эт.» из описания."""
    v = props.get("max_floors")
    if v not in (None, ""):
        try:
            return int(v)
        except (TypeError, ValueError):
            pass
    nums = [int(n) for n in re.findall(r"(\d+)\s*эт", str(props.get("description") or props.get("floors") or ""))]
    return max(nums) if nums else None


def make_site(poly, layover_m: float, dx: float, dy: float):
    """«Размазать» полигон в сторону спутника — сумма Минковского с отрезком длины layover_m."""
    steps = max(1, int(math.ceil(layover_m / 5)))
    return unary_union([affinity.translate(poly, dx * layover_m * k / steps, dy * layover_m * k / steps)
                        for k in range(steps + 1)])


def build_geometries(cx: list[Complex], crs, layover_m: float, toward_az: float,
                     ring_inner: float, ring_outer: float, incidence_deg: float = 40.0, floor_m: float = 3.2,
                     by_floors: bool = False):
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    to_utm = lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1]))
    dx, dy = math.sin(math.radians(toward_az)), math.cos(math.radians(toward_az))
    for c in cx:
        c.poly = shapely.transform(c.geom_wgs, to_utm)
        # Отражение от верха здания высотой H смещается к спутнику на H / tan(угла падения):
        # 25 этажей ≈ 80 м → ~95 м при 40°. На реальных данных Душанбе (14.11.2026) широкая маска
        # РАЗМЫЛА сигнал строек соседними кварталами, поэтому по умолчанию выключено (--layover-by-floors).
        fl = floors_of(c.props)
        lay = max(layover_m, fl * floor_m / math.tan(math.radians(incidence_deg))) if (fl and by_floors) else layover_m
        c.metrics["max_floors"], c.metrics["layover_m"] = fl, lay
        c.site = make_site(c.poly, lay, dx, dy)
    all_sites = unary_union([c.site.buffer(20) for c in cx])
    for c in cx:
        c.ring = c.poly.buffer(ring_outer).difference(c.poly.buffer(ring_inner)).difference(all_sites)


_mask_cache: dict = {}


def mask_for(geom, sc: Scene) -> np.ndarray:
    """Булева маска геометрии на сетке сцены (растеризуется только окно вокруг объекта)."""
    key = (geom.wkb, tuple(sc.transform), sc.shape)
    if key in _mask_cache:
        return _mask_cache[key]
    m = np.zeros(sc.shape, bool)
    t = sc.transform
    minx, miny, maxx, maxy = geom.bounds
    cols = sorted(((minx - t.c) / t.a, (maxx - t.c) / t.a))
    rows = sorted(((miny - t.f) / t.e, (maxy - t.f) / t.e))
    c0, c1 = max(0, int(math.floor(cols[0])) - 1), min(sc.shape[1], int(math.ceil(cols[1])) + 1)
    r0, r1 = max(0, int(math.floor(rows[0])) - 1), min(sc.shape[0], int(math.ceil(rows[1])) + 1)
    if c1 > c0 and r1 > r0:
        wt = t * rasterio.Affine.translation(c0, r0)
        g = [mapping(geom)]
        sub = geometry_mask(g, out_shape=(r1 - r0, c1 - c0), transform=wt, invert=True)
        if not sub.any():  # объект меньше пикселя — берём все задетые пиксели
            sub = geometry_mask(g, out_shape=(r1 - r0, c1 - c0), transform=wt, invert=True, all_touched=True)
        m[r0:r1, c0:c1] = sub
    _mask_cache[key] = m
    return m


# --------------------------------------------------------------------------- #
def series_indices(norm: np.ndarray) -> tuple[float, float]:
    v = norm[np.isfinite(norm)]
    if v.size < 2:
        return float("nan"), float("nan")
    k = 2 if v.size >= 4 else 1
    sci = float(np.mean(v[-k:]) - np.mean(v[:k]))
    avi = float(np.mean(np.abs(np.diff(v))))
    return sci, avi


def analyze_gcov(c: Complex, scenes: list[Scene]):
    for sc in scenes:
        ms, mr = mask_for(c.site, sc), mask_for(c.ring, sc)
        row = dict(complex_id=c.cid, obs_date=iso(sc.date), cycle_num=sc.cycle, granule_id=sc.granule,
                   n_site_px=int(ms.sum()), n_ring_px=int(mr.sum()))
        for pol, band in (("hh", "HHHH"), ("hv", "HVHV")):
            if band not in sc.bands:
                continue
            a = sc.bands[band]
            s, r = a[ms], a[mr]
            s, r = s[np.isfinite(s) & (s > 0)], r[np.isfinite(r) & (r > 0)]
            if s.size == 0 or r.size == 0:
                continue
            row[f"raw_{pol}_db"] = float(to_db(s.mean()))
            row[f"bg_{pol}_db"] = float(to_db(np.median(r)))
            row[f"norm_{pol}_db"] = row[f"raw_{pol}_db"] - row[f"bg_{pol}_db"]
            row[f"std_{pol}_db"] = float(np.std(to_db(s)))
        c.series.append(row)
    norm = np.array([r.get("norm_hh_db", np.nan) for r in c.series])
    c.metrics["sci_db"], c.metrics["avi_db"] = series_indices(norm)
    c.metrics["n_dates"] = int(np.isfinite(norm).sum())


def dedupe_gunw(scenes: list[Scene]) -> list[Scene]:
    """Одна сцена на пару дат: с unwrappedPhase и с наибольшим числом слоёв."""
    best: dict = {}
    for sc in scenes:
        key = (sc.date, sc.date2)
        score = ("unwrappedPhase" in sc.bands, len(sc.bands))
        if key not in best or score > best[key][0]:
            best[key] = (score, sc)
    out = [v[1] for _, v in sorted(best.items())]
    if len(out) < len(scenes):
        log.info("GUNW: %d файлов → %d пар (дубликаты отброшены)", len(scenes), len(out))
    return out


def remove_plane(ph: np.ndarray, sc: Scene, min_coh: float = 0.4) -> tuple[np.ndarray, float]:
    """Вычесть плоскость a + b·x + c·y, подобранную по надёжным пикселям города.

    Убирает плавные градиенты (остаток ионосферы, тропосфера, орбита), которые
    на 20 км дают миллиметры «ложного» смещения между районами.
    Возвращает фазу и модуль градиента, рад/км.
    """
    rows, cols = np.indices(ph.shape)
    km = sc.pixel_m / 1000.0
    ok = np.isfinite(ph)
    coh = sc.bands.get("coherenceMagnitude")
    if coh is not None:
        good = ok & (np.nan_to_num(coh) >= min_coh)
        ok = good if good.sum() > 100 else ok
    if ok.sum() < 10:
        return ph, float("nan")
    a = np.column_stack([np.ones(ok.sum()), cols[ok] * km, rows[ok] * km])
    coef, *_ = np.linalg.lstsq(a, ph[ok], rcond=None)
    plane = coef[0] + coef[1] * cols * km + coef[2] * rows * km
    return ph - plane, float(math.hypot(coef[1], coef[2]))


def prepare_gunw(sc: Scene, detrend: bool = True) -> np.ndarray | None:
    """Фаза → смещение LOS (мм) относительно медианы своей связной области."""
    ph = sc.bands.get("unwrappedPhase")
    if ph is None:
        return None
    ph = ph.astype(np.float64).copy()
    tag = f"GUNW {iso(sc.date)}→{iso(sc.date2)}"
    stats = [f"phase {np.isfinite(ph).mean():.0%}"]
    iono = sc.bands.get("ionospherePhaseScreen")
    if iono is not None:
        ok = np.isfinite(iono)
        stats.append(f"iono {ok.mean():.0%} (med {np.nanmedian(iono) if ok.any() else float('nan'):.2f} рад)")
        if ok.mean() >= 0.5:
            ph = np.where(ok, ph - iono, np.nan)
        else:
            log.warning("%s: ionospherePhaseScreen почти пуст — не вычитаем", tag)
    cc = sc.bands.get("connectedComponents")
    if cc is not None:
        cc = np.nan_to_num(cc).astype(np.int64)
        labels = np.unique(cc[(cc > 0) & np.isfinite(ph)])
        stats.append(f"cc>0 {(cc > 0).mean():.0%}, областей {labels.size}")
        if labels.size:
            ph[cc == 0] = np.nan
            for lab in labels:
                m = cc == lab
                ph[m] -= np.nanmedian(ph[m])
        else:
            log.warning("%s: connectedComponents без валидных областей — вычитаем общую медиану", tag)
            ph -= np.nanmedian(ph)
    else:
        ph -= np.nanmedian(ph)
    if detrend:
        ph, grad = remove_plane(ph, sc)
        stats.append(f"плоскость {grad * abs(PHASE_TO_MM):.2f} мм/км")
    stats.append(f"итог {np.isfinite(ph).mean():.0%}")
    log.info("%s: %s", tag, ", ".join(stats))
    return ph * PHASE_TO_MM


def analyze_gunw(c: Complex, scenes: list[Scene], disp: dict):
    total, n = 0.0, 0
    for sc in scenes:
        ms, mr = mask_for(c.site, sc), mask_for(c.ring, sc)
        row = dict(complex_id=c.cid, date1=iso(sc.date), date2=iso(sc.date2))
        coh = sc.bands.get("coherenceMagnitude")
        if coh is not None:
            row["coh_site"] = float(np.nanmean(coh[ms])) if ms.any() else float("nan")
            row["coh_ring"] = float(np.nanmean(coh[mr])) if mr.any() else float("nan")
        d = disp.get(sc.path)
        if d is not None and mr.any() and np.isfinite(d[mr]).any():
            row["shift_mm"] = float(np.nanmedian(d[mr]))
            total += row["shift_mm"]
            n += 1
        c.pairs.append(row)
    cs = [p.get("coh_site", np.nan) for p in c.pairs]
    cr = [p.get("coh_ring", np.nan) for p in c.pairs]
    c.metrics["coh_site"] = float(np.nanmean(cs)) if np.isfinite(cs).any() else None
    c.metrics["coh_ring"] = float(np.nanmean(cr)) if np.isfinite(cr).any() else None
    if c.metrics["coh_site"] is not None and c.metrics["coh_ring"] is not None:
        c.metrics["coh_drop"] = c.metrics["coh_ring"] - c.metrics["coh_site"]
    c.metrics["ground_shift_mm"] = total if n else None
    c.metrics["n_pairs"] = n


# --------------------------------------------------------------------------- #
def _moving_avg(a: np.ndarray, w: int, axis: int) -> np.ndarray:
    """Скользящее среднее окна w вдоль оси (края — повтор крайних значений)."""
    if w <= 1:
        return a
    a = np.moveaxis(a, axis, -1)
    left = w // 2
    ap = np.pad(a, [(0, 0)] * (a.ndim - 1) + [(left, w - 1 - left)], mode="edge")
    c = np.cumsum(ap, axis=-1, dtype=np.complex128 if np.iscomplexobj(a) else np.float64)
    c = np.concatenate([np.zeros_like(c[..., :1]), c], axis=-1)
    return np.moveaxis((c[..., w:] - c[..., :-w]) / w, -1, axis)


def boxcar(a: np.ndarray, wy: int, wx: int) -> np.ndarray:
    return _moving_avg(_moving_avg(a, wy, 0), wx, 1)


def _is_true(v) -> bool | None:
    if v is None:
        return None
    return str(v).strip().lower() in ("true", "1", "yes")


def prepare_gslc(scenes: list[Scene], window_m: float) -> dict | None:
    """Когерентность и мощность по паре GSLC (первая и последняя дата)."""
    if len(scenes) < 2:
        if scenes:
            log.info("GSLC: только одна дата (%s) — для когерентности нужна пара", iso(scenes[0].date))
        return None
    s1, s2 = scenes[0], scenes[-1]
    if s1.transform != s2.transform or s1.shape != s2.shape:
        log.warning("GSLC %s и %s на разных сетках — пропускаем", iso(s1.date), iso(s2.date))
        return None
    for sc in (s1, s2):
        flat = _is_true(sc.tags.get("topographicFlatteningApplied"))
        if flat is False:
            log.warning("GSLC %s: topographicFlatteningApplied=False — в интерферограмме останутся полосы "
                        "рельефа, когерентность будет занижена", iso(sc.date))
    rx, ry = abs(s1.transform.a), abs(s1.transform.e)
    wx, wy = max(1, int(round(window_m / rx))), max(1, int(round(window_m / ry)))
    g = dict(s1=s1, s2=s2, wx=wx, wy=wy, coh={})
    for pol in ("HH", "HV"):
        if pol in s1.bands and pol in s2.bands:
            a = np.nan_to_num(s1.bands[pol])
            b = np.nan_to_num(s2.bands[pol])
            num = boxcar(a * np.conj(b), wy, wx)
            den = np.sqrt(boxcar(np.abs(a) ** 2, wy, wx) * boxcar(np.abs(b) ** 2, wy, wx))
            with np.errstate(all="ignore"):
                g["coh"][pol] = np.where(den > 0, np.abs(num) / den, np.nan).astype(np.float32)
    if "HH" not in g["coh"]:
        log.warning("GSLC: нет слоя HH в обеих датах — пропускаем")
        return None
    g["p1"] = (np.abs(s1.bands["HH"]) ** 2).astype(np.float32)
    g["p2"] = (np.abs(s2.bands["HH"]) ** 2).astype(np.float32)
    log.info("GSLC %s→%s: пиксель %gx%g м, окно %dx%d px (~%g м), когерентность HH медиана %.2f; "
             "flattening: topo=%s, ellipsoid=%s",
             iso(s1.date), iso(s2.date), rx, ry, wx, wy, window_m, float(np.nanmedian(g["coh"]["HH"])),
             s1.tags.get("topographicFlatteningApplied", "?"), s1.tags.get("ellipsoidalFlatteningApplied", "?"))
    return g


def analyze_gslc(c: Complex, g: dict):
    s1 = g["s1"]
    ms, mr = mask_for(c.site, s1), mask_for(c.ring, s1)
    m = c.metrics
    m["gslc_dates"] = f"{iso(s1.date)}→{iso(g['s2'].date)}"
    m["gslc_pixels"] = int(ms.sum())
    if not ms.any() or not mr.any():
        return
    coh = g["coh"]["HH"]
    m["gslc_coh_site"] = float(np.nanmean(coh[ms]))
    m["gslc_coh_ring"] = float(np.nanmean(coh[mr]))
    m["gslc_coh_drop"] = m["gslc_coh_ring"] - m["gslc_coh_site"]
    norm = []
    for p in (g["p1"], g["p2"]):
        norm.append(float(to_db(np.nanmean(p[ms])) - to_db(np.nanmedian(p[mr]))))
    m["gslc_norm_hh_db"] = norm
    m["gslc_sci_db"] = norm[1] - norm[0]


def gslc_noise(g: dict, patch_m: float, ring_inner: float, ring_outer: float,
               n: int = 400, seed: int = 2) -> dict | None:
    """Падение когерентности и изменение яркости GSLC на случайных площадках города."""
    s1 = g["s1"]
    rx, ry = abs(s1.transform.a), abs(s1.transform.e)
    h, w = s1.shape
    px, py = max(1, int(round(patch_m / rx))), max(1, int(round(patch_m / ry)))
    rin_x, rin_y = int(math.ceil((ring_inner + patch_m / 2) / rx)), int(math.ceil((ring_inner + patch_m / 2) / ry))
    rout_x, rout_y = int(math.ceil((ring_outer + patch_m / 2) / rx)), int(math.ceil((ring_outer + patch_m / 2) / ry))
    if 2 * rout_y + 1 >= h or 2 * rout_x + 1 >= w:
        return None
    rng = np.random.default_rng(seed)
    coh, p1, p2 = g["coh"]["HH"], g["p1"], g["p2"]
    drops, scis = [], []
    for _ in range(n):
        r = int(rng.integers(rout_y, h - rout_y))
        c = int(rng.integers(rout_x, w - rout_x))
        site = (slice(r - py // 2, r - py // 2 + py), slice(c - px // 2, c - px // 2 + px))
        win = (slice(r - rout_y, r + rout_y + 1), slice(c - rout_x, c + rout_x + 1))
        inner = np.ones((2 * rout_y + 1, 2 * rout_x + 1), bool)
        inner[rout_y - rin_y:rout_y + rin_y + 1, rout_x - rin_x:rout_x + rin_x + 1] = False
        with np.errstate(all="ignore"):
            drops.append(np.nanmean(coh[win][inner]) - np.nanmean(coh[site]))
            n1 = to_db(np.nanmean(p1[site])) - to_db(np.nanmedian(p1[win][inner]))
            n2 = to_db(np.nanmean(p2[site])) - to_db(np.nanmedian(p2[win][inner]))
        scis.append(abs(n2 - n1))
    pct = lambda a: {f"p{q}": round(float(np.nanpercentile(a, q)), 3) for q in (50, 75, 90, 95)}
    return dict(patch_px=[py, px], coh_drop=pct(np.array(drops)), abs_sci_db=pct(np.array(scis)))


# --------------------------------------------------------------------------- #
def city_noise(scenes: list[Scene], patch_px: int, ring_in_px: int, ring_out_px: int,
               n: int = 400, seed: int = 0) -> dict:
    """Те же SCI/AVI для случайных площадок по городу — оценка фонового шума."""
    rng = np.random.default_rng(seed)
    h, w = scenes[0].shape
    half = patch_px // 2
    stack = np.stack([to_db(np.where(sc.bands["HHHH"] > 0, sc.bands["HHHH"], np.nan)) for sc in scenes])
    lin = np.stack([sc.bands["HHHH"] for sc in scenes])
    scis, avis = [], []
    for _ in range(n):
        r = int(rng.integers(ring_out_px, h - ring_out_px))
        col = int(rng.integers(ring_out_px, w - ring_out_px))
        site = lin[:, r - half:r - half + patch_px, col - half:col - half + patch_px]
        ring = lin[:, r - ring_out_px:r + ring_out_px + 1, col - ring_out_px:col + ring_out_px + 1].copy()
        c0 = ring_out_px - ring_in_px
        ring[:, c0:c0 + 2 * ring_in_px + 1, c0:c0 + 2 * ring_in_px + 1] = np.nan
        with np.errstate(all="ignore"):
            norm = to_db(np.nanmean(site, axis=(1, 2))) - to_db(np.nanmedian(ring, axis=(1, 2)))
        sci, avi = series_indices(norm)
        if np.isfinite(avi):
            scis.append(abs(sci))
            avis.append(avi)
    del stack
    scis, avis = np.array(scis), np.array(avis)
    pct = lambda a: {f"p{q}": round(float(np.percentile(a, q)), 3) for q in (50, 75, 90, 95)}
    return dict(n_patches=int(avis.size), patch_px=patch_px, abs_sci_db=pct(scis), avi_db=pct(avis))


def gunw_noise(scenes: list[Scene], disp: dict, patch_m: float, ring_inner: float,
               ring_outer: float, n: int = 400, seed: int = 1) -> dict | None:
    """Случайные площадки по городу на сетке GUNW: падение когерентности и
    накопленное смещение в кольце — фон, с которым сравниваем ЖК."""
    if not scenes:
        return None
    rng = np.random.default_rng(seed)
    px = scenes[0].pixel_m
    h, w = scenes[0].shape
    p = max(1, int(round(patch_m / px)))
    rin = int(math.ceil((ring_inner + patch_m / 2) / px))
    rout = int(math.ceil((ring_outer + patch_m / 2) / px))
    if 2 * rout + 1 >= min(h, w):
        return None
    cohs = [sc.bands.get("coherenceMagnitude") for sc in scenes]
    maps = [disp.get(sc.path) for sc in scenes]
    drops, shifts = [], []
    for _ in range(n):
        r = int(rng.integers(rout, h - rout))
        c = int(rng.integers(rout, w - rout))
        sl_site = (slice(r - p // 2, r - p // 2 + p), slice(c - p // 2, c - p // 2 + p))
        ring = np.zeros((h, w), bool)
        ring[r - rout:r + rout + 1, c - rout:c + rout + 1] = True
        ring[r - rin:r + rin + 1, c - rin:c + rin + 1] = False
        cs, cr, tot, k = [], [], 0.0, 0
        with np.errstate(all="ignore"):
            for coh, d in zip(cohs, maps):
                if coh is not None:
                    cs.append(np.nanmean(coh[sl_site]))
                    cr.append(np.nanmean(coh[ring]))
                if d is not None and np.isfinite(d[ring]).any():
                    tot += float(np.nanmedian(d[ring]))
                    k += 1
        if cs and np.isfinite(cs).any() and np.isfinite(cr).any():
            drops.append(np.nanmean(cr) - np.nanmean(cs))
        if k:
            shifts.append(tot)
    out = dict(patch_px=p, ring_px=[rin, rout])
    pct = lambda a: {f"p{q}": round(float(np.percentile(a, q)), 3) for q in (50, 75, 90, 95)}
    if drops:
        out["coh_drop"] = pct(np.array(drops))
    if shifts:
        sh = np.array(shifts)
        out["shift_mm"] = pct(np.abs(sh))
        out["shift_sigma_mm"] = round(float(1.4826 * np.median(np.abs(sh - np.median(sh)))), 3)
    return out


INDICATORS = ("gcov_sci", "gcov_avi", "gunw_coherence", "gslc_coherence", "gslc_sci")


def indicator_values(m: dict) -> dict:
    """Значения признаков активности (больше = сильнее признак)."""
    raw = {"gcov_sci": abs(m["sci_db"]) if m.get("sci_db") is not None else None,
           "gcov_avi": m.get("avi_db"),
           "gunw_coherence": m.get("coh_drop"),
           "gslc_coherence": m.get("gslc_coh_drop"),
           "gslc_sci": abs(m["gslc_sci_db"]) if m.get("gslc_sci_db") is not None else None}
    return {k: float(v) for k, v in raw.items() if v is not None and np.isfinite(v)}


def virtual_complexes(n: int, ref: Scene, side_m: float, layover_m: float, toward_az: float,
                      ring_inner: float, ring_outer: float, seed: int = 7) -> list[Complex]:
    """Случайные квадраты типичного размера ЖК по всему городу — «нулевая гипотеза»."""
    rng = np.random.default_rng(seed)
    w, s_, e, n_ = rasterio.transform.array_bounds(ref.shape[0], ref.shape[1], ref.transform)
    mg = ring_outer + side_m + layover_m
    dx, dy = math.sin(math.radians(toward_az)), math.cos(math.radians(toward_az))
    out = []
    for i in range(n):
        x, y = rng.uniform(w + mg, e - mg), rng.uniform(s_ + mg, n_ - mg)
        poly = shp_box(x - side_m / 2, y - side_m / 2, x + side_m / 2, y + side_m / 2)
        v = Complex(f"v{i}", {}, None, poly=poly, site=make_site(poly, layover_m, dx, dy))
        v.ring = poly.buffer(ring_outer).difference(poly.buffer(ring_inner)).difference(v.site.buffer(20))
        out.append(v)
    return out


def calibrate(virtuals: list[Complex], gcov, gunw, disp, gs, thr: dict) -> dict:
    """Те же расчёты на виртуальных ЖК: распределения признаков, оценка Фишера,
    доля ложных срабатываний правила active."""
    vals = {k: [] for k in INDICATORS}
    per = []
    for v in virtuals:
        analyze_gcov(v, gcov)
        if gunw:
            analyze_gunw(v, gunw, disp)
        if gs:
            analyze_gslc(v, gs)
        iv = indicator_values(v.metrics)
        per.append(iv)
        for k, x in iv.items():
            vals[k].append(x)
        _mask_cache.clear()  # маски виртуальных ЖК больше не нужны
    noise = {k: np.sort(np.array(a)) for k, a in vals.items() if len(a) >= 20}
    fisher = np.array([fisher_score(iv, noise) for iv in per])
    fpr_active = float(np.mean([any(iv.get(k, -np.inf) > t for k, t in thr.items() if t) for iv in per]))
    return dict(noise=noise, fisher=np.sort(fisher), n=len(per), fpr_active=fpr_active)


def pvalue(noise: np.ndarray, x: float) -> float:
    """Доля виртуальных ЖК с признаком не слабее x (эмпирический p-value)."""
    return (1 + (len(noise) - np.searchsorted(noise, x, side="left"))) / (len(noise) + 1)


def fisher_score(iv: dict, noise: dict) -> float:
    """Метод Фишера: −2·Σ ln p. Растёт, когда несколько признаков повышены одновременно."""
    return float(-2 * sum(math.log(pvalue(noise[k], x)) for k, x in iv.items() if k in noise))


def classify(c: Complex, min_px: float, t_sci: float, t_avi: float, t_coh: float,
             noise_avi: float, n_total_dates: int, t_gcoh: float | None = None, t_gsci: float | None = None,
             possible_pct: float = 90.0, active_pct: float = 95.0, threshold_rule: bool = False):
    m = c.metrics
    npx = m["nisar_pixels"]
    if npx < min_px or not np.isfinite(m.get("avi_db", np.nan)):
        m["sar_verdict"], m["confidence_score"] = "insufficient_resolution", 0
        return
    ratios = {"gcov_sci": abs(m["sci_db"]) / t_sci, "gcov_avi": m["avi_db"] / t_avi}
    if m.get("coh_drop") is not None:
        ratios["gunw_coherence"] = max(m["coh_drop"], 0) / t_coh
    if m.get("gslc_coh_drop") is not None and t_gcoh:
        ratios["gslc_coherence"] = max(m["gslc_coh_drop"], 0) / t_gcoh
    if m.get("gslc_sci_db") is not None and t_gsci:
        ratios["gslc_sci"] = abs(m["gslc_sci_db"]) / t_gsci
    # Какие признаки превысили порог — для объяснения вердикта на фронте/в боте.
    m["triggers"] = ",".join(k for k, v in ratios.items() if v > 1.0)
    strength = max(ratios.values())
    pct = m.get("activity_percentile")
    f_pix = min(1.0, math.sqrt(npx / 25.0))
    f_t = min(1.0, m["n_dates"] / max(n_total_dates, 1))
    if threshold_rule or pct is None:
        # Старое правило: любой признак выше своего порога p95.
        active, possible = strength > 1.0, False
        f_sig = min(1.0, 0.5 + 0.5 * (strength - 1.0)) if active else min(1.0, 0.5 + 0.5 * (1.0 - strength))
    else:
        # Калиброванное правило: совокупность признаков против участков того же размера.
        active = pct >= active_pct
        possible = not active and pct >= possible_pct
        if active:
            f_sig = min(1.0, 0.5 + 0.5 * (pct - active_pct) / max(100 - active_pct, 1))
        elif possible:
            f_sig = 0.5 + 0.5 * (pct - possible_pct) / max(active_pct - possible_pct, 1)
        else:
            f_sig = 0.5 + 0.5 * (possible_pct - pct) / max(possible_pct, 1)
    m["sar_verdict"] = "activity_detected" if active else "possible_activity" if possible else "no_changes"
    m["confidence_score"] = int(round(100 * f_pix * f_t * f_sig))
    m["snr"] = m["avi_db"] / noise_avi if noise_avi else None


def ground_risk(mm, sigma=None):
    """Смещение меньше 5 мм или в пределах 2σ фона города — stable."""
    if mm is None:
        return None
    a = abs(mm)
    if a < 5 or (sigma and a < 2 * sigma):
        return "stable"
    return "moderate_shift" if a < 15 else "high_risk"


# --------------------------------------------------------------------------- #
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _window(geom, sc: Scene, margin: float):
    minx, miny, maxx, maxy = geom.buffer(margin).bounds
    win = from_bounds(minx, miny, maxx, maxy, sc.transform).round_offsets().round_lengths()
    r0, c0 = max(0, win.row_off), max(0, win.col_off)
    r1, c1 = min(sc.shape[0], r0 + win.height), min(sc.shape[1], c0 + win.width)
    t = sc.transform
    extent = (t.c + c0 * t.a, t.c + c1 * t.a, t.f + r1 * t.e, t.f + r0 * t.e)
    return (slice(r0, r1), slice(c0, c1)), extent


def _outline(ax, geom, **kw):
    for g in getattr(geom, "geoms", [geom]):
        if g.is_empty:
            continue
        x, y = g.exterior.xy
        ax.plot(x, y, **kw)


def render_complex(c: Complex, scenes: list[Scene], out_dir: Path, vmin: float, vmax: float,
                   margin: float) -> dict:
    plt = _plt()
    urls = {}
    for sc, row in zip(scenes, c.series):
        sl, ext = _window(c.poly, sc, margin)
        img = to_db(sc.bands["HHHH"][sl])
        fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
        ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax, extent=ext, interpolation="nearest")
        _outline(ax, c.site, color="yellow", lw=1, ls="--")
        _outline(ax, c.poly, color="red", lw=1.5)
        ax.set_title(f"ЖК {c.cid} · {row['obs_date']} · HH", fontsize=9)
        ax.set_axis_off()
        name = f"jk_{c.cid}_{sc.date}_hh.png"
        fig.savefig(out_dir / name, bbox_inches="tight")
        plt.close(fig)
        row["crop_png_url"] = f"/static/crops/{name}"
    # Карта разницы: (последние k − первые k) − сдвиг фона.
    k = 2 if len(scenes) >= 4 else 1
    sl, ext = _window(c.poly, scenes[0], margin)
    first = to_db(np.mean([s.bands["HHHH"][sl] for s in scenes[:k]], axis=0))
    last = to_db(np.mean([s.bands["HHHH"][sl] for s in scenes[-k:]], axis=0))
    bg_shift = np.nanmean([r.get("bg_hh_db", np.nan) for r in c.series[-k:]]) - \
        np.nanmean([r.get("bg_hh_db", np.nan) for r in c.series[:k]])
    diff = last - first - bg_shift
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    im = ax.imshow(diff, cmap="RdBu_r", vmin=-6, vmax=6, extent=ext, interpolation="nearest")
    _outline(ax, c.poly, color="black", lw=1.5)
    ax.set_title(f"ЖК {c.cid} · Δ HH {iso(scenes[-1].date)} − {iso(scenes[0].date)}, дБ", fontsize=9)
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    name = f"jk_{c.cid}_diff.png"
    fig.savefig(out_dir / name, bbox_inches="tight")
    plt.close(fig)
    urls["diff_heatmap_url"] = f"/static/crops/{name}"
    return urls


def render_gslc_crop(c: Complex, g: dict, out_dir: Path, margin: float) -> dict:
    """Кроп когерентности GSLC вокруг ЖК (тёмное = поверхность изменилась)."""
    plt = _plt()
    sl, ext = _window(c.poly, g["s1"], margin)
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    im = ax.imshow(g["coh"]["HH"][sl], cmap="magma", vmin=0, vmax=1, extent=ext, interpolation="nearest")
    _outline(ax, c.poly, color="cyan", lw=1.5)
    ax.set_title(f"ЖК {c.cid} · когерентность GSLC {iso(g['s1'].date)} → {iso(g['s2'].date)}", fontsize=8)
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    name = f"jk_{c.cid}_gslc_coh.png"
    fig.savefig(out_dir / name, bbox_inches="tight")
    plt.close(fig)
    return {"gslc_coh_png_url": f"/static/crops/{name}"}


def render_city_layer(arr: np.ndarray, sc: Scene, path: Path, cmap: str, vmin: float, vmax: float) -> dict:
    """Перепроецировать в WGS84 и сохранить PNG с прозрачностью для ImageOverlay."""
    from rasterio.warp import Resampling, calculate_default_transform, reproject
    plt = _plt()
    h, w = sc.shape
    bounds = rasterio.transform.array_bounds(h, w, sc.transform)
    dst_t, dw, dh = calculate_default_transform(sc.crs, "EPSG:4326", w, h, *bounds)
    dst = np.full((dh, dw), np.nan, dtype=np.float32)
    reproject(arr.astype(np.float32), dst, src_transform=sc.transform, src_crs=sc.crs,
              dst_transform=dst_t, dst_crs="EPSG:4326", resampling=Resampling.bilinear,
              src_nodata=np.nan, dst_nodata=np.nan)
    norm = np.clip((dst - vmin) / (vmax - vmin), 0, 1)
    rgba = plt.get_cmap(cmap)(norm)
    rgba[..., 3] = np.where(np.isfinite(dst), 0.75, 0.0)
    plt.imsave(path, rgba)
    west, north = dst_t.c, dst_t.f
    east, south = west + dw * dst_t.a, north + dh * dst_t.e
    return dict(png_url=f"/static/city/{path.name}", bbox_west=west, bbox_south=south,
                bbox_east=east, bbox_north=north, vmin=vmin, vmax=vmax, cmap=cmap)


# --------------------------------------------------------------------------- #
def _clean(v):
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else round(float(v), 4)
    if isinstance(v, np.integer):
        return int(v)
    return v


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: _clean(r.get(k)) for k in keys} for r in rows])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Аналитика стройактивности по полигонам ЖК (NISAR GCOV/GUNW).",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("geojson", help="полигоны ЖК (FeatureCollection, WGS84)")
    ap.add_argument("--raw", default=str(DATA_ROOT / "raw_nisar"), help="папка с tif и manifest.csv")
    ap.add_argument("-o", "--output", default=str(DATA_ROOT), help="корень для analytics/ и static/")
    ap.add_argument("--layover-m", type=float, help="размазать полигон к спутнику, м (по умолч. = 1 пиксель)")
    ap.add_argument("--ring-inner", type=float, default=100.0, help="внутренний радиус фонового кольца, м")
    ap.add_argument("--ring-outer", type=float, default=400.0, help="внешний радиус фонового кольца, м")
    ap.add_argument("--min-pixels", type=float, default=4.0, help="меньше — insufficient_resolution")
    ap.add_argument("--t-sci", type=float, help="порог |SCI|, дБ (по умолч. p95 шума города)")
    ap.add_argument("--t-avi", type=float, help="порог AVI, дБ (по умолч. p95 шума города)")
    ap.add_argument("--t-coh", type=float, help="порог падения когерентности (по умолч. p95 шума города)")
    ap.add_argument("--plan-thresholds", action="store_true",
                    help="фиксированные пороги из плана: SCI 1.8, AVI 1.2, coh 0.25")
    ap.add_argument("--auto-thresholds", action="store_true", help=argparse.SUPPRESS)  # теперь по умолчанию
    ap.add_argument("--no-detrend", action="store_true", help="не вычитать плоскость из фазы GUNW")
    ap.add_argument("--gslc-window-m", type=float, default=40.0, help="окно усреднения когерентности GSLC, м")
    ap.add_argument("--layover-by-floors", action="store_true",
                    help="расширять маску к спутнику по этажности (H / tan угла падения); по умолч. — 1 пиксель")
    ap.add_argument("--incidence-deg", type=float, default=40.0,
                    help="угол падения радара для --layover-by-floors, градусы")
    ap.add_argument("--floor-m", type=float, default=3.2, help="высота этажа, м")
    ap.add_argument("--calib-n", type=int, default=200,
                    help="число случайных участков того же размера для калибровки каждого ЖК")
    ap.add_argument("--active-pct", type=float, default=95.0,
                    help="activity_detected, если совокупность признаков сильнее, чем у N%% таких участков")
    ap.add_argument("--possible-pct", type=float, default=90.0,
                    help="possible_activity, если сильнее, чем у N%% (но ниже --active-pct)")
    ap.add_argument("--threshold-rule", action="store_true",
                    help="старое правило: active, если любой признак выше порога p95")
    ap.add_argument("--t-gslc-coh", type=float, help="порог падения когерентности GSLC (по умолч. p95 шума)")
    ap.add_argument("--t-gslc-sci", type=float, help="порог |ΔдБ| GSLC (по умолч. p95 шума)")
    ap.add_argument("--no-png", action="store_true", help="не рисовать PNG")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    raw = Path(args.raw)
    gcov = load_scenes(raw, "GCOV")
    gunw = dedupe_gunw(load_scenes(raw, "GUNW"))
    gslc = dedupe_by_date(load_scenes(raw, "GSLC"))
    if not gcov:
        log.error("нет GCOV tif в %s — сначала запустите crop_dushanbe.py", raw)
        return 1
    if any(s.transform != gcov[0].transform or s.shape != gcov[0].shape for s in gcov):
        log.warning("сетки GCOV различаются — маски считаются для каждого файла отдельно")
    ref = gcov[0]
    direction = ref.direction if ref.direction in TOWARD_SAT_AZ else "A"
    layover = args.layover_m if args.layover_m is not None else ref.pixel_m
    log.info("GCOV: %d дат (%s … %s), пиксель %.0f м; GUNW: %d пар; орбита %s",
             len(gcov), iso(gcov[0].date), iso(gcov[-1].date), ref.pixel_m, len(gunw), direction)

    cx = load_complexes(Path(args.geojson))
    build_geometries(cx, ref.crs, layover, TOWARD_SAT_AZ[direction], args.ring_inner, args.ring_outer,
                     args.incidence_deg, args.floor_m, args.layover_by_floors)
    log.info("ЖК: %d", len(cx))

    # Шум фона города — площадка типичного размера ЖК.
    med_px = float(np.median([c.poly.area for c in cx])) / ref.pixel_area
    patch = max(3, int(round(math.sqrt(med_px))))
    noise = city_noise(gcov, patch, int(math.ceil(args.ring_inner / ref.pixel_m)) + patch // 2,
                       int(math.ceil(args.ring_outer / ref.pixel_m)) + patch // 2)
    disp = {sc.path: prepare_gunw(sc, detrend=not args.no_detrend) for sc in gunw}
    patch_m = math.sqrt(float(np.median([c.poly.area for c in cx])))
    gnoise = gunw_noise(gunw, disp, patch_m, args.ring_inner, args.ring_outer)
    gs = prepare_gslc(gslc, args.gslc_window_m)
    snoise = gslc_noise(gs, patch_m, args.ring_inner, args.ring_outer) if gs else None

    # Пороги: явно заданные > из плана (--plan-thresholds) > p95 шума города.
    plan = dict(sci=1.8, avi=1.2, coh=0.25)
    auto = dict(sci=noise["abs_sci_db"]["p95"], avi=noise["avi_db"]["p95"],
                coh=(gnoise or {}).get("coh_drop", {}).get("p95", plan["coh"]))
    base = plan if args.plan_thresholds else auto
    t_sci = args.t_sci if args.t_sci is not None else base["sci"]
    t_avi = args.t_avi if args.t_avi is not None else base["avi"]
    t_coh = args.t_coh if args.t_coh is not None else base["coh"]
    t_gcoh = args.t_gslc_coh if args.t_gslc_coh is not None else (
        plan["coh"] if args.plan_thresholds or not snoise else snoise["coh_drop"]["p95"])
    t_gsci = args.t_gslc_sci if args.t_gslc_sci is not None else (
        plan["sci"] if args.plan_thresholds or not snoise else snoise["abs_sci_db"]["p95"])
    shift_sigma = (gnoise or {}).get("shift_sigma_mm")
    log.info("Шум города GCOV (площадка %dx%d px): AVI p50=%.2f p95=%.2f дБ, |SCI| p50=%.2f p95=%.2f дБ",
             patch, patch, noise["avi_db"]["p50"], noise["avi_db"]["p95"],
             noise["abs_sci_db"]["p50"], noise["abs_sci_db"]["p95"])
    if gnoise:
        log.info("Шум города GUNW: падение когерентности p50=%.2f p95=%.2f, смещение кольца σ=%.1f мм (|.| p95=%.1f)",
                 gnoise.get("coh_drop", {}).get("p50", float("nan")), gnoise.get("coh_drop", {}).get("p95", float("nan")),
                 shift_sigma or float("nan"), gnoise.get("shift_mm", {}).get("p95", float("nan")))
    if snoise:
        log.info("Шум города GSLC (площадка %dx%d px): падение когерентности p50=%.2f p95=%.2f, "
                 "|ΔдБ| p50=%.2f p95=%.2f", snoise["patch_px"][0], snoise["patch_px"][1],
                 snoise["coh_drop"]["p50"], snoise["coh_drop"]["p95"],
                 snoise["abs_sci_db"]["p50"], snoise["abs_sci_db"]["p95"])
    log.info("Пороги (%s): SCI=%.2f дБ, AVI=%.2f дБ, падение когерентности=%.2f%s",
             "план" if args.plan_thresholds else "p95 шума города", t_sci, t_avi, t_coh,
             f"; GSLC: когерентность={t_gcoh:.2f}, ΔдБ={t_gsci:.2f}" if gs else "")
    out = Path(args.output)
    (out / "analytics").mkdir(parents=True, exist_ok=True)
    crops, city = out / "static" / "crops", out / "static" / "city"
    if not args.no_png:
        crops.mkdir(parents=True, exist_ok=True)
        city.mkdir(parents=True, exist_ok=True)
        allhh = to_db(np.stack([s.bands["HHHH"] for s in gcov]))
        vmin, vmax = np.nanpercentile(allhh, [2, 98])

    for c in cx:
        cen = c.geom_wgs.centroid
        c.metrics.update(area_sqm=c.poly.area, nisar_pixels=c.poly.area / ref.pixel_area,
                         center_lat=cen.y, center_lon=cen.x, pixel_m=ref.pixel_m,
                         first_obs_date=iso(gcov[0].date), last_obs_date=iso(gcov[-1].date))
        analyze_gcov(c, gcov)
        if gunw:
            analyze_gunw(c, gunw, disp)
        if gs:
            analyze_gslc(c, gs)

    # Калибровка: каждый ЖК против случайных участков города того же размера и маски навала.
    thr = dict(gcov_sci=t_sci, gcov_avi=t_avi, gunw_coherence=t_coh,
               gslc_coherence=t_gcoh if gs else None, gslc_sci=t_gsci if gs else None)
    calibs: dict = {}
    for c in cx:
        side, lay = math.sqrt(c.poly.area), c.metrics["layover_m"]
        key = (round(side / 10) * 10, round(lay / 10) * 10)  # похожие ЖК делят калибровку
        if key not in calibs:
            calibs[key] = calibrate(virtual_complexes(args.calib_n, ref, max(key[0], ref.pixel_m), max(key[1], 1.0),
                                                      TOWARD_SAT_AZ[direction], args.ring_inner, args.ring_outer,
                                                      seed=len(calibs) + 7), gcov, gunw, disp, gs, thr)
            log.info("Калибровка %d участков %.0f×%.0f м (навал %.0f м): старое правило «любой признак > p95» "
                     "срабатывает на %.0f%% случайных участков", calibs[key]["n"], key[0], key[0], key[1],
                     100 * calibs[key]["fpr_active"])
        c.calib = calibs[key]
    for c in cx:
        calib = c.calib
        iv = indicator_values(c.metrics)
        score = fisher_score(iv, calib["noise"])
        c.metrics["activity_score"] = score
        c.metrics["activity_percentile"] = float(100 * np.searchsorted(calib["fisher"], score, side="right")
                                                 / len(calib["fisher"]))
        # Признаки выше p75 города — «повышены», даже если ниже порога.
        c.metrics["elevated"] = ",".join(k for k, x in iv.items()
                                         if k in calib["noise"] and pvalue(calib["noise"][k], x) <= 0.25)
        c.metrics["fpr_threshold_rule"] = calib["fpr_active"]

    for c in cx:
        classify(c, args.min_pixels, t_sci, t_avi, t_coh, noise["avi_db"]["p50"], len(gcov),
                 t_gcoh if gs else None, t_gsci if gs else None,
                 args.possible_pct, args.active_pct, args.threshold_rule)
        shift = c.metrics.get("ground_shift_mm")
        c.metrics["ground_shift_sigma_mm"] = shift_sigma
        c.metrics["ground_shift_z"] = shift / shift_sigma if shift is not None and shift_sigma else None
        c.metrics["ground_risk_level"] = ground_risk(shift, shift_sigma)
        if not args.no_png and c.metrics["sar_verdict"] != "insufficient_resolution":
            c.metrics.update(render_complex(c, gcov, crops, vmin, vmax, margin=200))
            if gs:
                c.metrics.update(render_gslc_crop(c, gs, crops, margin=200))
        m = c.metrics
        log.info("ЖК %-4s %6.0f м² %5.1f px | SCI %+5.2f AVI %4.2f дБ | coh site/ring %s | сдвиг %s мм %s | %s (%d%%)",
                 c.cid, m["area_sqm"], m["nisar_pixels"], m.get("sci_db", np.nan), m.get("avi_db", np.nan),
                 f"{m['coh_site']:.2f}/{m['coh_ring']:.2f}" if m.get("coh_site") is not None else "—",
                 (f"{m['ground_shift_mm']:+.1f}" + (f"±{shift_sigma:.1f}" if shift_sigma else ""))
                 if m.get("ground_shift_mm") is not None else "—",
                 f"({m['ground_risk_level']}, z={m['ground_shift_z']:+.1f})" if m.get("ground_shift_z") is not None else "",
                 m["sar_verdict"], m["confidence_score"])
        if m.get("gslc_coh_site") is not None:
            log.info("        GSLC %s: %d px | coh site/ring %.2f/%.2f | ΔдБ %+.2f | признаки: %s",
                     m["gslc_dates"], m["gslc_pixels"], m["gslc_coh_site"], m["gslc_coh_ring"],
                     m["gslc_sci_db"], m.get("triggers") or "—")
        log.info("        навал %.0f м (%s эт.) | совокупность признаков сильнее, чем у %.0f%% кварталов | "
                 "повышены: %s", m["layover_m"], m.get("max_floors") or "?", m["activity_percentile"],
                 m.get("elevated") or "—")

    layers = []
    if not args.no_png:
        k = 2 if len(gcov) >= 4 else 1
        first = to_db(np.mean([s.bands["HHHH"] for s in gcov[:k]], axis=0))
        last = to_db(np.mean([s.bands["HHHH"] for s in gcov[-k:]], axis=0))
        change = (last - first) - np.nanmedian(last - first)
        layers.append(dict(layer_type="gcov_change_heatmap", date_from=iso(gcov[0].date), date_to=iso(gcov[-1].date),
                           **render_city_layer(change, ref, city / "dushanbe_gcov_change.png", "RdBu_r", -6, 6)))
        maps = [d for d in disp.values() if d is not None]
        if maps and all(d.shape == maps[0].shape for d in maps):
            cum = np.sum(maps, axis=0)
            lim = float(np.nanpercentile(np.abs(cum), 98)) or 10.0
            layers.append(dict(layer_type="gunw_ground_shift_map", date_from=iso(gunw[0].date),
                               date_to=iso(gunw[-1].date2), units="mm LOS, relative to city",
                               **render_city_layer(cum, gunw[0], city / "dushanbe_insar_shift.png",
                                                   "RdBu", -lim, lim)))
        if gs:
            layers.append(dict(layer_type="gslc_coherence_map", date_from=iso(gs["s1"].date),
                               date_to=iso(gs["s2"].date),
                               **render_city_layer(gs["coh"]["HH"], gs["s1"], city / "dushanbe_gslc_coherence.png",
                                                   "magma", 0.0, 1.0)))
        json.dump(layers, open(city / "layers.json", "w"), indent=2)

    # --- Выход ---
    feats = []
    for c in cx:
        props = {**c.props, "id": c.cid, **{k: _clean(v) for k, v in c.metrics.items()}}
        feats.append(dict(type="Feature", properties=props, geometry=mapping(c.geom_wgs)))
    json.dump(dict(type="FeatureCollection", features=feats),
              open(out / "analytics" / "complexes.geojson", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    write_csv(out / "analytics" / "timeseries.csv", [r for c in cx for r in c.series])
    write_csv(out / "analytics" / "gunw_pairs.csv", [r for c in cx for r in c.pairs])
    json.dump(dict(noise=noise, gunw_noise=gnoise, gslc_noise=snoise,
                   calibration=dict(rule="threshold" if args.threshold_rule else "calibrated",
                                    n_virtual_per_size=args.calib_n, active_percentile=args.active_pct,
                                    possible_percentile=args.possible_pct,
                                    sizes=[dict(side_m=k[0], layover_m=k[1], fpr_threshold_rule=round(v["fpr_active"], 3))
                                           for k, v in calibs.items()],
                                    incidence_deg=args.incidence_deg, floor_m=args.floor_m),
                   gslc=dict(dates=[iso(gs["s1"].date), iso(gs["s2"].date)], window_px=[gs["wy"], gs["wx"]],
                             thresholds=dict(coh_drop=t_gcoh, sci_db=t_gsci),
                             flattening=dict(topo=gs["s1"].tags.get("topographicFlatteningApplied"),
                                             ellipsoid=gs["s1"].tags.get("ellipsoidalFlatteningApplied")))
                   if gs else None,
                   thresholds=dict(source="plan" if args.plan_thresholds else "city_noise_p95",
                                   sci_db=t_sci, avi_db=t_avi, coh_drop=t_coh, min_pixels=args.min_pixels),
                   gunw_detrend=not args.no_detrend,
                   pixel_m=ref.pixel_m, layover_m=layover, ring_m=[args.ring_inner, args.ring_outer],
                   orbit_direction=direction, gcov_dates=[iso(s.date) for s in gcov],
                   gunw_pairs=[[iso(s.date), iso(s.date2)] for s in gunw],
                   gunw_layers=sorted(gunw[0].bands) if gunw else []),
              open(out / "analytics" / "city_stats.json", "w"), indent=2, ensure_ascii=False)
    log.info("Готово: %s/analytics, %s/static", out, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
