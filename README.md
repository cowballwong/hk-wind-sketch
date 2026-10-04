# HK Wind Sketch · 香港風環境草圖

See how wind moves between buildings at street level in Hong Kong, live in your browser.

Search any place or address in Hong Kong (or click the map, or pick one of four ready-made areas: Kwun Tong, Central, Mong Kok, Tsuen Wan). Your browser fetches the government's building footprints for a 600 m square there and a small wind simulation runs on a flat slice about 2 m above the ground. Pick a wind direction, draw a new block, move it, rotate it, or take an existing building away, and watch the street wind change, in plan or in a 3D view with the buildings at their real heights.

**This is a qualitative 2D sketch for teaching and early design feel.** It is not an Air Ventilation Assessment (AVA) under HPLB/ETWB TC No. 1/06, gives no wind velocity ratio (VR) figures, no comfort rating and no pass or fail, and says nothing about wind loads. See [DISCLAIMER.md](DISCLAIMER.md).

## Use it

- **Online:** https://anzonwong.com/hk-wind-sketch/
- **Locally:** serve the folder with any static server (for example `python -m http.server`) and open `index.html`. The page needs the internet for the government data services, the base map and the pinned CDN libraries (Leaflet 1.9.4, three.js 0.160.0).
- **As a small Flask app:** `pip install flask`, then `python app/app.py` and open http://127.0.0.1:5058/wind.

## How it works

- **Solver:** a D2Q9 lattice Boltzmann method (BGK collision with a Smagorinsky turbulence model), halfway bounce-back at walls, a uniform incoming wind applied through a thin sponge layer round the edge of the 600 m square. Grid of 240, 300 or 400 cells across (2.5, 2 or 1.5 m per cell).
- **Where it runs:** WebGL2 fragment shaders with float render targets; it asks the browser for the low-power GPU. If WebGL2 float targets are missing it falls back to plain JavaScript on a coarser 150 x 150 grid.
- **Kind to laptops:** capped at 30 frames per second, stops when the tab is hidden, and pauses itself after 90 seconds with no input.
- **Any site:** the page asks the CSDI web feature service for LandsD buildings in the square (bounding-box query, 1,000 per page), the PlanD 200 m wind tables for the 500 m cell at its centre (`app/static/data/pland_grid.json` maps cells to table files; rebuild with `prep_pland_grid.py`), and searches with the LandsD location search (map.gov.hk) plus the Address Lookup Service. All these services answer browser requests directly; there is no server of our own.
- **3D view:** three.js, loaded only when you open it. Blocks are extruded from footprint to TopHeight minus BaseHeight on a flat site model (no terrain); the coloured layer is the same 2D slice draped at 2 m, and the streaks follow it. It does not simulate 3D flow.
- **Height:** for the wind, used only to decide what stands in the 2 m slice. A block at least 2.5 m tall from its own base blocks it; open-sided structures do not. A 10-storey and a 60-storey tower look the same to a 2D slice.

## Data

| What | Source | Licence |
|---|---|---|
| Building footprints, base and top height (mPD), storeys, block type | Lands Department, HKSAR Government: Building dataset, CSDI Portal (`landsd_rcd_1637211194312_35158`, WFS) | DATA.GOV.HK Terms and Conditions |
| Wind direction frequency, annual / summer / winter, 200 m level, 500 m grid | Planning Department, HKSAR Government: Site Wind Availability Data | DATA.GOV.HK Terms and Conditions |
| Base map and labels (XYZ tiles) | Map from Lands Department (CSDI Map API) | CSDI Map API terms: LandsD logo and copyright notice on the map |
| Place search | Lands Department location search (map.gov.hk) | DATA.GOV.HK / CSDI terms |
| Address search | Address Lookup Service (ALS), HKSAR Government | ALS terms |

The data is provided by the Government as is, with no warranty. `prep_data.py` rebuilds the four ready-made area files from the live services (only the few hundred buildings inside each square are requested). Coordinates in the JSON are local metres, x east and y south, from the top-left corner of the square.

To add an area, add an entry to `AREAS` in `prep_data.py`, run it, and add a button in `app/templates/wind/_body.html`.

## Limitations

- 2D only: no terrain, trees, footbridges, covered walkways or podium decks; open ground floors are not known; air cannot go over roofs, so speeds in narrow gaps are exaggerated compared with real 3D flow.
- Downwash from tall towers is not represented.
- The incoming wind is uniform and has the same strength for every direction; the PlanD wind rose only suggests which directions to try.
- Cells are 1.5 to 2.5 m, so lanes under about two cells wide are poorly resolved.
- Ready-made areas are a snapshot (October 2026); searched sites use whatever the LandsD service returns that day.
- If a government service is down or slow, the site will not load; press Try again.

## Licence

Code: MIT (see [LICENSE](LICENSE)). Data: see the table above; keep the attribution when you reuse it.
