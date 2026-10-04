"""HK Wind Sketch v3: cut the LandsD 5 m DTM into small tiles the page can fetch.

    python prep_terrain.py [out_dir]     default C:/hk-wind-sketch-v2/app/static/data/dtm

Source: Lands Department, Digital Terrain Model 5 m grid of the whole territory
(DATA.GOV.HK dataset hk-landsd-openmap-5m-grid-dtm, CSDI landsd_rcd_1638158088368_93806),
https://www.landsd.gov.hk/landsd_psi_data/SMO/data/Whole_HK_DTM_5m.zip, cached at
C:/hk-wind-terrain-cache (the zip holds one ESRI ASCII grid, HK1980 Grid, metres PD,
whole metres, sea stored as 0).

Output: 1 km x 1 km tiles of 200 x 200 cells (5 m, native resolution), int16 metres,
little endian, rows from north to south, compressed with zlib (the browser inflates
them with DecompressionStream('deflate')). Only tiles that contain land are written;
a missing tile means sea (0 mPD). index.json lists them. Never run on G: (C: only).
"""
from __future__ import annotations

import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(r"C:/hk-wind-terrain-cache/Whole_HK_DTM_5m.asc")
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:/hk-wind-sketch-v2/app/static/data/dtm")
TILE = 200          # cells per tile side (1 km)
NPY = SRC.with_suffix(".npy")


def load() -> tuple[dict, np.ndarray]:
    hdr = {}
    with SRC.open("r") as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    if NPY.exists():
        return hdr, np.load(NPY)
    t = time.time()
    df = pd.read_csv(SRC, sep=r"\s+", header=None, skiprows=6, dtype=np.float32, engine="c")
    a = df.to_numpy()
    assert a.shape == (int(hdr["nrows"]), int(hdr["ncols"])), a.shape
    assert np.all(a == np.round(a)), "expected whole metres"
    a = a.astype(np.int16)
    np.save(NPY, a)
    print(f"parsed {SRC.name} in {time.time() - t:.0f} s: {a.shape}, {a.min()}..{a.max()} m")
    return hdr, a


def main() -> None:
    assert not str(OUT).lower().startswith("g:"), "write tiles on C:, not the Google Drive"
    hdr, a = load()
    nrows, ncols = a.shape
    nodata = int(hdr.get("nodata_value", -9999))
    a = np.where(a == nodata, 0, a).astype(np.int16)
    tcols, trows = -(-ncols // TILE), -(-nrows // TILE)
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.bin"):
        old.unlink()
    tiles, raw, packed = [], 0, 0
    top_n = hdr["yllcorner"] + nrows * hdr["cellsize"]      # north edge of the grid
    for tj in range(trows):              # tj counts from the north edge
        for ti in range(tcols):
            blk = np.zeros((TILE, TILE), np.int16)
            src = a[tj * TILE:(tj + 1) * TILE, ti * TILE:(ti + 1) * TILE]
            blk[:src.shape[0], :src.shape[1]] = src
            if not (blk > 0).any():
                continue
            b = blk.astype("<i2").tobytes()
            z = zlib.compress(b, 9)
            (OUT / f"{ti}_{tj}.bin").write_bytes(z)
            tiles.append(f"{ti}_{tj}")
            raw += len(b); packed += len(z)
    idx = {
        "source": "Lands Department, HKSAR Government: Digital Terrain Model 5 m grid (DATA.GOV.HK hk-landsd-openmap-5m-grid-dtm), derived from DOP5000 (Dec 2015 release); accuracy about +/-5 m; includes elevated roads and vegetation height",
        "crs": "HK1980 Grid (EPSG:2326), heights in metres above Principal Datum",
        "x0": hdr["xllcorner"], "north": top_n, "cell": hdr["cellsize"], "tile": TILE,
        "cols": tcols, "rows": trows,
        "order": "tile {ti}_{tj}: ti from the west edge x0, tj from the north edge; int16 little endian metres, rows north to south, zlib",
        "tiles": tiles,
    }
    (OUT / "index.json").write_text(json.dumps(idx, separators=(",", ":")), encoding="utf-8")
    sizes = [p.stat().st_size for p in OUT.glob("*.bin")]
    print(f"{len(tiles)} land tiles of {tcols * trows}; raw {raw / 1e6:.1f} MB, written {packed / 1e6:.1f} MB; "
          f"largest {max(sizes) / 1024:.0f} KB, median {sorted(sizes)[len(sizes) // 2] / 1024:.0f} KB; index {(OUT / 'index.json').stat().st_size / 1024:.0f} KB")
    # spot checks against values LandsD's other services gave (DTM file values: Peak 549 at E833034.52 N815227.378)
    def at(E, N):
        i = int((E - hdr["xllcorner"]) // hdr["cellsize"]); j = int((top_n - N) // hdr["cellsize"])
        return int(a[j, i])
    print("spot: Peak", at(833034.52, 815227.378))


if __name__ == "__main__":
    main()
