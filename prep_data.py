"""HK Wind Sketch: build the small per-area data files the page loads.

For each demo area (a ~600 m square) this fetches, from free government open data:
  * LandsD "Building" dataset (CSDI WFS, landsd_rcd_1637211194312_35158):
    footprints with BaseHeight / TopHeight (mPD), storeys, block type;
  * PlanD "Site Wind Availability Data" (500 m grid, RAMS simulation): the
    16-direction frequency table at 200 m for the grid cell at the area centre,
    for annual, summer and winter.
and writes app/static/data/<area>.json (local metres, x east, y south, origin at
the square's top-left corner).

Only the few features inside each square are requested. The 26 MB PlanD index
CSV lives in C:/hk-wind-sketch-cache (never on G:). Re-run any time:
    python prep_data.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "app" / "static" / "data"
CACHE = Path(r"C:/hk-wind-sketch-cache")
WFS = ("https://portal.csdi.gov.hk/server/services/common/"
       "landsd_rcd_1637211194312_35158/MapServer/WFSServer")
PLAND_CSV_URL = "https://www.pland.gov.hk/pland_en/resources/site_wind/open_data/site_wind_open_data.csv"
SIZE_M = 600.0       # simulated square
FETCH_PAD_M = 20.0   # fetch a little wider so edge blocks are complete

AREAS = [
    {"id": "kwun_tong", "name_en": "Kwun Tong", "name_zh": "觀塘",
     "lat": 22.3112, "lon": 114.2243,
     "note_en": "Former industrial grid, now office towers on deep podiums. Streets run about 40 degrees off north.",
     "note_zh": "前工業區街網，而家多數係高座寫字樓加大平台，街道同正北差大約 40 度。"},
    {"id": "central", "name_en": "Central", "name_zh": "中環",
     "lat": 22.2815, "lon": 114.1578,
     "note_en": "Tall, closely spaced towers on a slope; narrow streets between Queen's Road and Des Voeux Road.",
     "note_zh": "斜坡上密集高樓，皇后大道同德輔道之間街道窄。"},
    {"id": "mong_kok", "name_en": "Mong Kok", "name_zh": "旺角",
     "lat": 22.3188, "lon": 114.1697,
     "note_en": "Dense tenement grid either side of Nathan Road: long street walls, few gaps.",
     "note_zh": "彌敦道兩旁密集唐樓街網，街牆長，空隙少。"},
    {"id": "tsuen_wan", "name_en": "Tsuen Wan", "name_zh": "荃灣",
     "lat": 22.3712, "lon": 114.1150,
     "note_en": "New-town mix of estates, malls and older streets near the MTR.",
     "note_zh": "新市鎮：屋苑、商場同舊街混合，近港鐵站。"},
]

M_PER_DEG_LAT = 110574.0


def m_per_deg_lon(lat: float) -> float:
    return 111320.0 * math.cos(math.radians(lat))


def get(url: str, tries: int = 3, timeout: int = 90) -> bytes:
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "hk-wind-sketch-prep/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def fetch_buildings(a: dict) -> list[dict]:
    half = SIZE_M / 2 + FETCH_PAD_M
    dlat = half / M_PER_DEG_LAT
    dlon = half / m_per_deg_lon(a["lat"])
    s, n = a["lat"] - dlat, a["lat"] + dlat
    w, e = a["lon"] - dlon, a["lon"] + dlon
    feats, start = [], 0
    while True:
        url = (f"{WFS}?service=wfs&version=2.0.0&request=GetFeature&typeNames=csdi:Building"
               f"&outputFormat=GEOJSON&count=1000&startIndex={start}"
               f"&bbox={s:.6f},{w:.6f},{n:.6f},{e:.6f},urn:ogc:def:crs:EPSG::4326")
        d = json.loads(get(url).decode("utf-8"))
        batch = d.get("features", [])
        feats += batch
        if len(batch) < 1000:
            break
        start += 1000
        time.sleep(1)
    return feats


def to_local(a: dict, lon: float, lat: float) -> tuple[float, float]:
    half = SIZE_M / 2
    x = (lon - a["lon"]) * m_per_deg_lon(a["lat"]) + half
    y = -(lat - a["lat"]) * M_PER_DEG_LAT + half       # y grows southwards (screen)
    return round(x, 2), round(y, 2)


def convert(a: dict, feats: list[dict]) -> tuple[list[dict], dict]:
    out, types, skipped = [], {}, 0
    for f in feats:
        p = f.get("properties") or {}
        g = f.get("geometry") or {}
        if p.get("Status") not in (None, "Active"):
            skipped += 1
            continue
        base, top = p.get("BaseHeight"), p.get("TopHeight")
        if base is None or top is None:
            skipped += 1
            continue
        h = float(top) - float(base)
        polys = g.get("coordinates") or []
        if g.get("type") == "Polygon":
            polys = [polys]
        rings = []
        for poly in polys:
            if not poly:
                continue
            ring = [to_local(a, c[0], c[1]) for c in poly[0]]   # outer ring only
            rings.append(ring)
        if not rings:
            skipped += 1
            continue
        t = p.get("BuildingBlockType") or ""
        types[t] = types.get(t, 0) + 1
        out.append({
            "id": p.get("OBJECTID"),
            "name": (p.get("BuildingNameEN") or "").strip(),
            "name_zh": (p.get("BuildingNameTC") or "").strip(),
            "type": t,
            "base": round(float(base), 1),
            "top": round(float(top), 1),
            "h": round(h, 1),
            "storeys": p.get("Storeys"),
            "rings": rings,
        })
    return out, {"block_types": types, "skipped": skipped}


def pland_rows() -> list[dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / "site_wind_open_data.csv"
    if not p.exists():
        p.write_bytes(get(PLAND_CSV_URL, timeout=300))
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def site_wind(a: dict, rows: list[dict]) -> dict:
    # NB the CSV's MIN_LAT/MIN_LNG columns are swapped (MIN_LAT holds longitude).
    hit = {}
    for r in rows:
        try:
            lng0, lat0 = float(r["MIN_LAT"]), float(r["MIN_LNG"])
            lng1, lat1 = float(r["MAX_LAT"]), float(r["MAX_LNG"])
        except (KeyError, ValueError):
            continue
        if lng0 <= a["lon"] < lng1 and lat0 <= a["lat"] < lat1 and r["HEIGHT"].strip() == "200":
            hit[r["RANGE_TYPE"].strip()] = r
    res = {"height_m": 200, "seasons": {}}
    for season in ("annual", "summer", "winter"):
        r = hit.get(season)
        if not r:
            continue
        url = r["FILE_URL"].strip()
        txt = get(url).decode("utf-8", "replace").splitlines()
        dirs = txt[0].split("\t")[2:18]
        sums = [float(v) for v in txt[1].split("\t")[2:18]]
        res["grid_id"] = r["GRID_ID"]
        res["cell_hk1980"] = [float(r["MIN_X"]), float(r["MIN_Y"]), float(r["MAX_X"]), float(r["MAX_Y"])]
        res["seasons"][season] = {"dirs": dirs, "freq": sums, "source": url}
        time.sleep(0.5)
    return res


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = pland_rows()
    index = []
    for a in AREAS:
        feats = fetch_buildings(a)
        blds, stats = convert(a, feats)
        wind = site_wind(a, rows)
        data = {
            "id": a["id"], "name_en": a["name_en"], "name_zh": a["name_zh"],
            "note_en": a["note_en"], "note_zh": a["note_zh"],
            "centre_wgs84": [a["lat"], a["lon"]], "size_m": SIZE_M,
            "buildings": blds, "site_wind": wind,
            "sources": {
                "buildings": "Lands Department, HKSAR Government: Building dataset (CSDI Portal, landsd_rcd_1637211194312_35158), via DATA.GOV.HK",
                "wind": "Planning Department, HKSAR Government: Site Wind Availability Data (CSDI / DATA.GOV.HK)",
                "fetched": date.today().isoformat(),
            },
            "stats": stats,
        }
        (OUT / f"{a['id']}.json").write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
        size = (OUT / f"{a['id']}.json").stat().st_size
        print(f"{a['id']}: {len(feats)} features -> {len(blds)} blocks, {stats}, wind seasons {list(wind['seasons'])}, {size/1024:.0f} KB")
        index.append({"id": a["id"], "name_en": a["name_en"], "name_zh": a["name_zh"]})
        time.sleep(1)
    (OUT / "areas.json").write_text(json.dumps(index, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
