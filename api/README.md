# SiteSAR API (FastAPI)

Отдаёт фронтенду и Telegram-боту результаты `h5_scripts/analyze_complexes.py`:
статус стройки по радару, ряд сигнала по датам, когерентность, смещение грунта,
PNG-кропы и общегородские слои. БД не нужна — API читает файлы из `<repo>/../data`
и сам перечитывает их, когда аналитика их обновит.

## Запуск
```bash
cd api
docker compose up -d --build        # http://localhost:8080/docs — интерактивная документация
docker compose logs -f              # логи
```
Внешний адрес: **https://apex-api.code.tj** (nginx → контейнер :8080). Фронтенд ходит туда по умолчанию,
на localhost — в `http://localhost:8080`; переопределить можно `?api=…` в адресе страницы.

Переменные (необязательные): `HONAINAV_DATA` — папка data (по умолчанию `../../data`),
`PUBLIC_URL` — внешний адрес API для ссылок на картинки (по умолчанию `https://apex-api.code.tj`;
для локального запуска: `PUBLIC_URL= docker compose up -d`), `CORS_ORIGINS` — разрешённые домены (по умолчанию `*`).

Каталог ЖК — `geojson/complexes_dushanbe.geojson` (id, название, адрес, статус по реестру и по проверке).
Аналитику запускать на нём же, чтобы id совпадали:
```bash
python3 h5_scripts/analyze_complexes.py geojson/complexes_dushanbe.geojson
```

## Эндпоинты
| Метод | Что отдаёт | Для чего на фронте |
|---|---|---|
| `GET /api/v1/meta` | даты снимков, пары GUNW, даты GSLC, пороги, подписи статусов, сводка по статусам | шкала времени, плашки «строятся / на паузе», легенда |
| `GET /api/v1/complexes` | GeoJSON FeatureCollection: полигоны + краткие свойства | слой на карте, список, поиск. Фильтры: `?status=active&status=possible`, `?q=Рудаки`, `?bbox=w,s,e,n` |
| `GET /api/v1/complexes/{id}` | карточка: всё из списка + `reasons`, `timeseries`, `gunw_pairs`, `gslc`, `thresholds`, `raw` | правая панель ЖК, график, слайдер снимков |
| `GET /api/v1/complexes/{id}/timeseries` | только ряд по датам, столбцами | Chart.js |
| `GET /api/v1/layers` | общегородские PNG + `coordinates` для MapLibre | переключатель слоёв |
| `GET /api/v1/nearby?lat=&lon=&radius_m=1000` | ближайшие ЖК с расстоянием | бот: «отправить геопозицию» |
| `GET /static/...` | PNG-кропы и слои | `<img>`, MapLibre image source |
| `GET /health` | загружены ли данные | мониторинг |

## Формат
Свойства ЖК (`properties` в списке и корень карточки):
```jsonc
{
  "id": "zhk-001", "name": "ЖК «Бульвар Скандинавия» (Sebiston)",
  "district": "Сино", "address": "ул. Панджруд, 91–92 мкр", "description": "блоки В 22 эт., Г 24 эт.",
  "registry_status": "green",          // реестр разрешений
  "reported_status": "active",         // проверка на месте — для сверки «радар vs реальность»
  "status": "active",                  // active | possible | no_changes | insufficient | pending
  "status_label": "Стройка идёт",
  "confidence": 73,                    // 0..100
  "activity_percentile": 99.5,         // сильнее, чем у N% случайных участков того же размера (≥95 active, ≥90 possible)
  "elevated": ["gunw_coherence", "gslc_coherence"],  // признаки выше p75 города (объяснение для possible)
  "triggers": ["gunw_coherence", "gslc_coherence", "gslc_sci"],
  "center": [68.7484, 38.5901],        // [lon, lat]
  "area_sqm": 14388, "pixels": {"gcov_20m": 36.0, "gslc_10x5m": 377},
  "signal": {
    "trend_db_per_month": 0.31,        // наклон очищенного сигнала
    "sci_db": 1.17, "avi_db": 0.54,
    "sparkline_db": [ ... ],           // по meta.dates
    "progress": [0, 0.1, ..., 1]       // 0..1 для анимации роста на карте (НЕ высота дома)
  },
  "coherence": {"site": 0.26, "background": 0.54, "source": "gslc_10x5m"},
  "ground": {"shift_mm": 2.3, "sigma_mm": 3.1, "z": 0.7, "risk": "stable", "risk_label": "Грунт стабилен"},
  "images": {"diff_png": "http://…/static/crops/jk_zhk-001_diff.png", "gslc_coherence_png": "…"}
}
```
Ряд (`timeseries`) — столбцами, индекс = дата: `dates`, `site_hh_db`, `background_hh_db`,
`normalized_hh_db` (участок − фон), то же для HV, `crop_png_url` (снимок на каждую дату).

## Как заменить демо-данные во фронтенде (`frontend/index.html`)
| Сейчас (демо) | Из API |
|---|---|
| `GEO` | `GET /api/v1/complexes` |
| `it.on` (по `truth`) | `properties.status === 'active'` (или `'possible'`) |
| `DATES` | `meta.dates` → формат `ДД.ММ` |
| `d.site`, `d.bg` | `timeseries.site_hh_db`, `timeseries.background_hh_db` |
| `d.trend` | `signal.trend_db_per_month` |
| `d.coh` | `coherence.site` |
| `d.ground` | `ground.shift_mm` (± `ground.sigma_mm`) |
| `d.h` (высота для 3D) | `signal.progress[i] × этажность × 3.2 м` — только визуализация |
| текст «почему» | `reasons` из карточки |

```js
const API = 'https://apex-api.code.tj/api/v1';
const [meta, fc] = await Promise.all([fetch(`${API}/meta`).then(r => r.json()),
                                      fetch(`${API}/complexes`).then(r => r.json())]);
const card = await fetch(`${API}/complexes/zhk-001`).then(r => r.json());
// MapLibre: слой когерентности
const L = (await fetch(`${API}/layers`).then(r => r.json())).find(l => l.id === 'gslc_coherence_map');
map.addSource(L.id, {type: 'image', url: L.url, coordinates: L.coordinates});
map.addLayer({id: L.id, type: 'raster', source: L.id, paint: {'raster-opacity': 0.7}});
```
