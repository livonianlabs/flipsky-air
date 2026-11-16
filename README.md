# FlipSky air quality

Air quality forecast tiles for the FlipSky weather widget.

Twice a day a GitHub Action takes the latest 00 UTC run of the **CAMS European air quality forecasts**
(ensemble, surface, hourly, 0–72 h), cuts it into 1° tiles and publishes them on GitHub Pages
(`t/<lat>_<lon>.bin.gz`, format described in `cams_tiles.py`).

Contains modified Copernicus Atmosphere Monitoring Service information (CC-BY 4.0).
Neither the European Commission nor ECMWF is responsible for any use that may be made of the
Copernicus information or data it contains.

## Pollen and bathing places

The same job cuts the CAMS pollen forecast (alder, birch, grasses, mugwort, olive, ragweed) into `p/<lat>_<lon>.bin.gz`.

`b/<lat>_<lon>.json` lists the EU's bathing places (name, position, kind of water, water quality class) from the
Bathing Water Directive register: European Environment Agency data via EMODnet Human Activities (CC BY 4.0).
`beaches.py` rebuilds them; run it by hand once a year and commit `beaches/`.

## Sea temperature

`sea_tiles.py` cuts the Copernicus Marine global (about 9 km) and Baltic (about 2 km) sea temperature forecasts, today and
9 days ahead, into `s/` (5° tiles) and `sb/` (1° tiles). It needs the `CMEMS_USER` and `CMEMS_PASS` secrets (a free
Copernicus Marine account) and is skipped without them. Generated using E.U. Copernicus Marine Service Information.

## The rest of the world's air

`global_air.py` cuts the CAMS global atmospheric composition forecast (about 0.4°, every 3 hours, 0-72 h) into `ga/`
(5° tiles, format version 2 in the script). The app uses them where the European tiles have nothing.

## The chance of rain worldwide

`ecmwf_rain.py` takes ECMWF's open ensemble probabilities (tpg1: at least 1 mm in 24 hours, 10 days, 0.25°) by byte
ranges and cuts them into `rp/` (5° tiles). Source: ECMWF open data, CC BY 4.0.

## Wind gusts worldwide

`ecmwf_gust.py` takes 10fg (the strongest gust in each 3 hours, 0.25°, 3-72 h) from ECMWF's open high-resolution run
and cuts it into `gu/` (5° tiles). Source: ECMWF open data, CC BY 4.0.

## The world's beaches

`beaches_world.py` adds the world's named beaches from OpenStreetMap (via Overture Maps' base theme) to `beaches/`,
after `beaches.py`: about 37,000 outside the EU register. Beaches © OpenStreetMap contributors (ODbL), via Overture Maps.

## Licences of the published data

- `b/` (bathing places): the EU register's part is EEA data via EMODnet (CC BY 4.0); the beaches from OpenStreetMap
  form a database derived from OpenStreetMap and are made available under the Open Database License (ODbL 1.0),
  © OpenStreetMap contributors, via Overture Maps.
- `t/`, `p/`, `ga/`: Contains modified Copernicus Atmosphere Monitoring Service information.
- `s/`, `sb/`, `sal/`, `sc/`: Generated using E.U. Copernicus Marine Service Information.
- `rp/`, `gu/`: ECMWF open data, CC BY 4.0.
- `st/`: NOAA GHCN-Daily, public domain.
