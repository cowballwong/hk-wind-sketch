# HK Wind Sketch · 香港風環境草圖

See how wind moves over the terrain and around the buildings of any place in Hong Kong, live in your browser, in 3D.

Search any place or address in Hong Kong (or click the map, or pick one of four ready-made areas: Kwun Tong, Central, Mong Kok, Tsuen Wan). Your browser fetches the government's buildings and ground levels for a 600 m square there and a small 3D wind simulation runs on your graphics card. Pick a wind direction, draw your own building outline corner by corner, edit it, or take an existing building away, and watch the wind change: 3D streamlines over the roofs and round the towers, a wind colour on the ground about 2 m up, and a vertical section along the wind.

**This is a qualitative teaching sketch on a coarse grid, for early design feel.** It is not an Air Ventilation Assessment (AVA) under HPLB/ETWB TC No. 1/06, gives no wind velocity ratio (VR) figures, no comfort rating and no pass or fail, and says nothing about wind loads. The 3D speeds are approximate. See [DISCLAIMER.md](DISCLAIMER.md).

## Use it

- **Online:** https://anzonwong.com/hk-wind-sketch/
- **Locally:** serve the folder with any static server (for example `python -m http.server`) and open `index.html`. The page needs the internet for the government data services, the base map and the pinned CDN libraries (Leaflet 1.9.4, three.js 0.160.0). The terrain tiles are in this repository.
- **As a small Flask app:** `pip install flask`, then `python app/app.py` and open http://127.0.0.1:5058/wind.

## How it works

- **3D solver** (`app/static/solver3d.js`): lattice Boltzmann, D3Q19, regularised BGK collision with a Smagorinsky eddy viscosity, pull streaming, half-way bounce-back on the ground and buildings, a free-slip lid, and the incoming wind applied as an equilibrium with a power-law profile (exponent 0.25, reaching the reference speed at 200 m, PlanD's level) on the four sides through a thin sponge layer, so any of the 16 directions works on one grid. Own code, written for this project.
- **Grid:** 128 x 128 cells across the 600 m square at Normal (4.7 m cells; Light 96, Fine 160), and as many layers as the terrain and buildings need. Total cells are capped (about 1.15 million at Normal); on steep or very tall sites the cells get bigger so the hill and the towers still fit under the lid.
- **Geometry:** a 2.5D height field. Each ground column is solid from the terrain up to the tallest roof standing on it (LandsD TopHeight in mPD). Drawn buildings stand on the lowest ground under their outline. Adding, moving or removing a building updates the height field on the GPU without restarting the flow.
- **Where it runs:** WebGPU compute shaders when available; otherwise WebGL2 fragment shaders with the layers tiled into a 2D atlas and five float render targets. Chrome on Windows ignores WebGPU's low-power request, so if WebGPU would pick a different graphics card from the low-power WebGL2 one (for example a laptop's discrete GPU), the page uses WebGL2 instead. Both back ends give the same numbers. The step rate is measured once per site and the solver uses about half of each frame.
- **Fallback:** on devices with neither, and by default on phones, the page runs the v2 sketch: a D2Q9 slice about 2 m above the ground (WebGL2, or plain JavaScript on a coarse grid), and says so under the picture. Phones can switch to 3D.
- **Drawing:** click corner by corner; double-click, click the first corner, press Enter or right-click to finish; press and drag for a rectangle; Esc cancels, Backspace takes back a corner. Snapping to corners of drawn buildings, right angles to the previous edge (Shift flips it) and a 1 m grid. Select a drawn building to move it, drag corners, drag an edge's midpoint to add a corner, right-click a corner to delete it, or press Delete.
- **Display:** the ground colour is the speed in the first air cell above the terrain (about 2 to 3 m up) compared with the incoming wind at the same height; streamlines start from a rake along the upwind edge at 3 to 100 m above the ground; particles drift with the 3D field; the section is a vertical slice along the wind.
- **Kind to laptops:** capped at 30 frames per second, the solver takes about half of each frame, it stops when the tab is hidden and pauses after 90 seconds with no input.
- **Any site:** the page asks the CSDI web feature service for LandsD buildings in the square, the PlanD 200 m wind tables for the 500 m cell at its centre (`app/static/data/pland_grid.json`, rebuilt with `prep_pland_grid.py`), and the terrain tiles under the square from `app/static/data/dtm/` (rebuilt with `prep_terrain.py`). Search uses the LandsD location search (map.gov.hk) and the Address Lookup Service. There is no server of our own.

## Data

| What | Source | Licence |
|---|---|---|
| Building footprints, base and top height (mPD), storeys, block type | Lands Department, HKSAR Government: Building dataset, CSDI Portal (`landsd_rcd_1637211194312_35158`, WFS) | DATA.GOV.HK Terms and Conditions |
| Ground levels | Lands Department, HKSAR Government: Digital Terrain Model 5 m grid of the whole territory (DATA.GOV.HK `hk-landsd-openmap-5m-grid-dtm`, CSDI `landsd_rcd_1638158088368_93806`), derived from the December 2015 DOP5000 release, about +/-5 m. Cut into 1,602 tiles of 1 km (200 x 200 cells, int16 metres, zlib), 23 MB, in `app/static/data/dtm/` | DATA.GOV.HK Terms and Conditions |
| Wind direction frequency, annual / summer / winter, 200 m level, 500 m grid | Planning Department, HKSAR Government: Site Wind Availability Data | DATA.GOV.HK Terms and Conditions |
| Base map and labels (XYZ tiles) | Map from Lands Department (CSDI Map API) | CSDI Map API terms: LandsD logo and copyright notice on the map |
| Place search | Lands Department location search (map.gov.hk) | DATA.GOV.HK / CSDI terms |
| Address search | Address Lookup Service (ALS), HKSAR Government | ALS terms |

The data is provided by the Government as is, with no warranty. `prep_data.py` rebuilds the four ready-made area files from the live services. Coordinates in the JSON are local metres, x east and y south, from the top-left corner of the square.

## Limitations

- Coarse: cells are 4 to 6 m, so lanes and gaps narrower than about two cells close up, and small buildings can disappear or grow.
- Simple solids: no open ground floors, footbridges, covered walkways or trees as separate objects. The terrain model includes elevated roads and tree canopy, which become solid ground.
- The incoming wind profile is assumed (power law, same strength for every direction); the PlanD wind rose only suggests which directions to try. The sides of the square are pulled towards the incoming wind, so the flow near the edges is not free.
- Towers or hills above the lid (about 190 to 700 m up, depending on the site) are cut off there; the page says so when it happens.
- Not checked against a wind tunnel or a reference CFD case. Density stays within about 2% of its starting value in the runs we checked, which is fine for a sketch and not for numbers.
- Ready-made areas are a snapshot (October 2026); searched sites use whatever the LandsD service returns that day; the terrain model is from 2015.
- If a government service is down or slow, the site will not load; press Try again.

## Licence

Code: MIT (see [LICENSE](LICENSE)). Data: see the table above; keep the attribution when you reuse it.
