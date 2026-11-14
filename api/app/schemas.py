"""Схемы ответов API — по ним строится документация /docs."""

from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field

Status = Literal["active", "possible", "no_changes", "insufficient", "pending"]


class Signal(BaseModel):
    trend_db_per_month: Optional[float] = Field(
        None, description="Наклон очищенного сигнала HH (участок минус фон), дБ в месяц")
    sci_db: Optional[float] = Field(None, description="Structural Change Index: последние 2 даты минус первые 2, дБ")
    avi_db: Optional[float] = Field(None, description="Activity Volatility Index: средний скачок между датами, дБ")
    sparkline_db: list[Optional[float]] = Field(
        default_factory=list, description="Очищенный сигнал по датам meta.dates — для мини-графика")
    progress: list[float] = Field(
        default_factory=list,
        description="0..1 по датам: накопленный рост сигнала для анимации на карте (это НЕ высота дома)")


class Coherence(BaseModel):
    site: Optional[float] = Field(None, description="Средняя когерентность на участке (0..1; ниже = поверхность менялась)")
    background: Optional[float] = Field(None, description="То же в кольце 100–400 м вокруг")
    source: Optional[Literal["gslc_10x5m", "gunw_80m"]] = None


class Ground(BaseModel):
    shift_mm: Optional[float] = Field(None, description="Накопленное смещение грунта в кольце вокруг ЖК, мм по лучу зрения, относительно города")
    sigma_mm: Optional[float] = Field(None, description="Погрешность (разброс по случайным кварталам), мм")
    z: Optional[float] = Field(None, description="shift_mm / sigma_mm; |z| < 2 — в пределах шума")
    risk: Optional[Literal["stable", "moderate_shift", "high_risk"]] = None
    risk_label: Optional[str] = None


class Images(BaseModel):
    diff_png: Optional[str] = Field(None, description="Карта изменения яркости вокруг ЖК (последние − первые даты)")
    gslc_coherence_png: Optional[str] = Field(None, description="Кроп когерентности 10×5 м вокруг ЖК")


class Pixels(BaseModel):
    gcov_20m: Optional[float] = None
    gslc_10x5m: Optional[int] = None


class ComplexSummary(BaseModel):
    id: str
    name: str
    district: Optional[str] = None
    address: Optional[str] = None
    description: Optional[str] = None
    registry_status: Optional[str] = Field(None, description="Статус в реестре разрешений (green / none / …)")
    reported_status: Optional[str] = Field(None, description="Статус по проверке на месте (active / stalled / completed) — для сверки с радаром")
    status: Status = Field(description="Вывод радара для интерфейса")
    status_label: str
    verdict: Optional[str] = Field(None, description="Исходный вердикт аналитики")
    confidence: Optional[int] = Field(None, description="Достоверность вывода, 0..100")
    activity_percentile: Optional[float] = Field(
        None, description="Совокупность признаков сильнее, чем у N% случайных участков того же размера "
                          "(≥95 → active, ≥90 → possible)")
    triggers: list[str] = Field(default_factory=list, description="Какие признаки по отдельности превысили порог p95")
    elevated: list[str] = Field(default_factory=list, description="Какие признаки повышены (выше p75 города)")
    center: list[float] = Field(description="[lon, lat]")
    area_sqm: Optional[float] = None
    pixels: Pixels
    signal: Signal
    coherence: Coherence
    ground: Ground
    images: Images


class Timeseries(BaseModel):
    dates: list[str]
    site_hh_db: list[Optional[float]] = Field(description="Сигнал HH на участке, дБ")
    background_hh_db: list[Optional[float]] = Field(description="Медиана HH в фоновом кольце, дБ")
    normalized_hh_db: list[Optional[float]] = Field(description="Участок минус фон, дБ (без влияния погоды)")
    site_hv_db: list[Optional[float]]
    background_hv_db: list[Optional[float]]
    normalized_hv_db: list[Optional[float]]
    site_std_hh_db: list[Optional[float]]
    crop_png_url: list[Optional[str]] = Field(description="Снимок HH вокруг ЖК на каждую дату")


class GunwPair(BaseModel):
    date1: str
    date2: str
    coh_site: Optional[float] = None
    coh_background: Optional[float] = None
    shift_mm: Optional[float] = None


class Gslc(BaseModel):
    dates: Optional[str] = None
    pixels: Optional[int] = None
    coh_site: Optional[float] = None
    coh_background: Optional[float] = None
    coh_drop: Optional[float] = None
    delta_db: Optional[float] = None


class ComplexDetail(ComplexSummary):
    reasons: list[str] = Field(default_factory=list, description="Почему такой вывод — готовые фразы для карточки/бота")
    timeseries: Timeseries
    gunw_pairs: list[GunwPair]
    gslc: Optional[Gslc] = None
    thresholds: dict[str, Any] = Field(default_factory=dict)
    raw: dict[str, Any] = Field(default_factory=dict, description="Все поля аналитики без изменений")


class Feature(BaseModel):
    type: Literal["Feature"] = "Feature"
    id: str
    geometry: dict[str, Any]
    properties: ComplexSummary


class FeatureCollection(BaseModel):
    type: Literal["FeatureCollection"] = "FeatureCollection"
    features: list[Feature]


class CityLayer(BaseModel):
    id: str
    title: str
    url: str
    date_from: Optional[str] = None
    date_to: Optional[str] = None
    units: Optional[str] = None
    bbox: list[float] = Field(description="[west, south, east, north]")
    coordinates: list[list[float]] = Field(description="Углы для MapLibre image source: СЗ, СВ, ЮВ, ЮЗ")
    legend: dict[str, Any]


class Stats(BaseModel):
    total: int
    by_status: dict[str, int]


class Meta(BaseModel):
    updated_at: Optional[str] = Field(None, description="Когда аналитика последний раз обновила результаты (UTC)")
    satellite: str = "NASA-ISRO NISAR (L-band)"
    orbit: Optional[str] = None
    dates: list[str] = Field(description="Даты снимков GCOV — ось X графиков и шкала времени")
    gunw_pairs: list[list[str]] = []
    gslc_dates: list[str] = []
    pixel_m: Optional[float] = None
    thresholds: dict[str, Any] = {}
    statuses: dict[str, str] = Field(description="Статус -> подпись для легенды")
    stats: Stats


class NearbyItem(BaseModel):
    distance_m: float
    complex: ComplexSummary


class Health(BaseModel):
    status: Literal["ok"] = "ok"
    data_loaded: bool
    complexes: int
    updated_at: Optional[str] = None
