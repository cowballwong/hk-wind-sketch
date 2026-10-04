"""HK Wind Sketch v2: compact index of PlanD's Site Wind Availability grid.

The page fetches the 200 m wind tables for any site straight from pland.gov.hk
(CORS is open), but it needs to know which table file belongs to which 500 m
cell. PlanD's index CSV (26 MB, cached at C:/hk-wind-sketch-cache) says; this
packs it into a small regular grid: HK1980 origin, 500 m step, and for every
cell the table name such as "d_01664" ("" = no cell, e.g. far out at sea).

    python prep_pland_grid.py [out.json]   (default src/data/pland_grid.json;
                                            in the repo: app/static/data/pland_grid.json)
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CSV = Path(r"C:/hk-wind-sketch-cache/site_wind_open_data.csv")
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "src" / "data" / "pland_grid.json"


def main() -> None:
    rows = list(csv.DictReader(CSV.open(encoding="utf-8-sig", newline="")))
    cells: dict[tuple[float, float], set[str]] = {}
    for r in rows:
        m = re.search(r"/([a-z]+_\d+)_lev200\.txt", r["FILE_URL"])
        if not m or r["HEIGHT"].strip() != "200":
            continue
        cells.setdefault((float(r["MIN_X"]), float(r["MIN_Y"])), set()).add(m.group(1))
    bad = {k: v for k, v in cells.items() if len(v) != 1}
    assert not bad, f"seasons disagree on table number: {list(bad.items())[:3]}"
    x0 = min(k[0] for k in cells); y0 = min(k[1] for k in cells)
    cols = int(round((max(k[0] for k in cells) - x0) / 500)) + 1
    nrows = int(round((max(k[1] for k in cells) - y0) / 500)) + 1
    ids = [""] * (cols * nrows)
    for (x, y), v in cells.items():
        c = int(round((x - x0) / 500)); r = int(round((y - y0) / 500))
        ids[r * cols + c] = next(iter(v))
    prefixes = sorted({next(iter(v)).split("_")[0] for v in cells.values()})
    sample = f"prefixes {prefixes}"
    out = {
        "source": "Planning Department, HKSAR Government: Site Wind Availability Data (index CSV via DATA.GOV.HK)",
        "crs": "HK1980 Grid (EPSG:2326), cell lower-left corners", "x0": x0, "y0": y0, "step": 500,
        "cols": cols, "rows": nrows,
        "url": "https://www.pland.gov.hk/pland_en/resources/site_wind/site_wind/table/{season}/{id}_lev200.txt",
        "ids": ids,
    }
    OUT.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    print(f"{len(cells)} cells, grid {cols} x {nrows}, {OUT.stat().st_size/1024:.0f} KB; sample {sample}")


if __name__ == "__main__":
    main()
