"""
SiteSAR API — данные спутникового анализа строек Душанбе для фронтенда и бота.

Запуск: cd api && docker compose up -d --build   →   http://localhost:8080/docs
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import RedirectResponse
from fastapi.staticfiles import StaticFiles

from . import schemas
from .store import STATUS_LABEL_RU, Store, city_layers, detail, haversine_m, summary

DATA_DIR = Path(os.environ.get("DATA_DIR", "/data"))
CATALOG_PATH = os.environ.get("CATALOG_PATH") or None
PUBLIC_URL = os.environ.get("PUBLIC_URL", "").rstrip("/")  # напр. https://api.example.tj
CORS_ORIGINS = [o.strip() for o in os.environ.get("CORS_ORIGINS", "*").split(",") if o.strip()]

store = Store(DATA_DIR, Path(CATALOG_PATH) if CATALOG_PATH else None)

app = FastAPI(
    title="SiteSAR API",
    version="1.0.0",
    description="Результаты анализа снимков NISAR по новостройкам Душанбе: статус стройки по радару, "
                "временные ряды сигнала, когерентность, смещение грунта и растровые слои для карты.",
)
app.add_middleware(CORSMiddleware, allow_origins=CORS_ORIGINS, allow_methods=["GET"], allow_headers=["*"])
app.add_middleware(GZipMiddleware, minimum_size=1000)


@app.middleware("http")
async def no_cache_api(request: Request, call_next):
    """Данные меняются после каждого пересчёта аналитики — не даём браузеру/прокси их кэшировать."""
    response = await call_next(request)
    if request.url.path.startswith("/api/") or request.url.path == "/health":
        response.headers["Cache-Control"] = "no-cache, must-revalidate"
    return response
# PNG-кропы и слои города, которые рисует analyze_complexes.py.
app.mount("/static", StaticFiles(directory=DATA_DIR / "static", check_dir=False), name="static")


def base_url(request: Request) -> str:
    """Абсолютный адрес API для ссылок на картинки (фронтенд живёт на другом порту)."""
    return PUBLIC_URL or str(request.base_url).rstrip("/")


def _get(cid: str):
    snap = store.get()
    if cid not in snap.complexes:
        raise HTTPException(404, f"ЖК '{cid}' не найден")
    return snap


@app.get("/", include_in_schema=False)
def root():
    return RedirectResponse("/docs")


@app.get("/health", response_model=schemas.Health, tags=["service"])
def health():
    snap = store.get()
    return {"data_loaded": snap.loaded, "complexes": len(snap.order), "updated_at": snap.updated_at}


@app.get("/api/v1/meta", response_model=schemas.Meta, tags=["meta"],
         summary="Даты снимков, пороги, легенда и сводка по статусам")
def meta():
    snap = store.get()
    st = snap.stats
    by_status: dict[str, int] = {}
    for cid in snap.order:
        s = summary(snap, cid)["status"]
        by_status[s] = by_status.get(s, 0) + 1
    gslc = st.get("gslc") or {}
    return {
        "updated_at": snap.updated_at,
        "orbit": {"A": "ascending", "D": "descending"}.get(st.get("orbit_direction")),
        "dates": st.get("gcov_dates", []),
        "gunw_pairs": st.get("gunw_pairs", []),
        "gslc_dates": gslc.get("dates", []),
        "pixel_m": st.get("pixel_m"),
        "thresholds": {**st.get("thresholds", {}), **({"gslc": gslc.get("thresholds")} if gslc else {})},
        "statuses": STATUS_LABEL_RU,
        "stats": {"total": len(snap.order), "by_status": by_status},
    }


@app.get("/api/v1/complexes", response_model=schemas.FeatureCollection, tags=["complexes"],
         summary="Все ЖК — GeoJSON для карты и списка")
def complexes(request: Request,
              status: Optional[list[schemas.Status]] = Query(None, description="Фильтр по статусу (можно несколько)"),
              q: Optional[str] = Query(None, description="Поиск по названию, адресу, району"),
              bbox: Optional[str] = Query(None, description="west,south,east,north")):
    snap = store.get()
    b = base_url(request)
    box = [float(x) for x in bbox.split(",")] if bbox else None
    feats = []
    for cid in snap.order:
        s = summary(snap, cid, b)
        if status and s["status"] not in status:
            continue
        if q and q.lower() not in " ".join(str(s.get(k) or "") for k in ("name", "address", "district")).lower():
            continue
        if box and not (box[0] <= s["center"][0] <= box[2] and box[1] <= s["center"][1] <= box[3]):
            continue
        feats.append({"type": "Feature", "id": cid, "geometry": snap.complexes[cid]["geometry"], "properties": s})
    return {"type": "FeatureCollection", "features": feats}


@app.get("/api/v1/complexes/{cid}", response_model=schemas.ComplexDetail, tags=["complexes"],
         summary="Карточка ЖК: вывод, причины, ряд по датам, когерентность, грунт")
def complex_detail(cid: str, request: Request):
    snap = _get(cid)
    return detail(snap, cid, base_url(request))


@app.get("/api/v1/complexes/{cid}/timeseries", response_model=schemas.Timeseries, tags=["complexes"],
         summary="Только ряд по датам — столбцами, сразу для Chart.js")
def complex_timeseries(cid: str, request: Request):
    snap = _get(cid)
    return detail(snap, cid, base_url(request))["timeseries"]


@app.get("/api/v1/layers", response_model=list[schemas.CityLayer], tags=["map"],
         summary="Общегородские растровые слои (PNG + углы для MapLibre)")
def layers(request: Request):
    return city_layers(store.get(), base_url(request))


@app.get("/api/v1/nearby", response_model=list[schemas.NearbyItem], tags=["complexes"],
         summary="Ближайшие ЖК к точке — для геолокации в Telegram-боте")
def nearby(request: Request,
           lat: float = Query(..., ge=-90, le=90), lon: float = Query(..., ge=-180, le=180),
           radius_m: float = Query(1000, gt=0, le=50000), limit: int = Query(5, ge=1, le=50)):
    snap = store.get()
    b = base_url(request)
    out = []
    for cid in snap.order:
        s = summary(snap, cid, b)
        d = haversine_m(lon, lat, s["center"][0], s["center"][1])
        if d <= radius_m:
            out.append({"distance_m": round(d, 1), "complex": s})
    return sorted(out, key=lambda x: x["distance_m"])[:limit]
