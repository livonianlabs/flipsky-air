# FlipSky air quality

Air quality forecast tiles for the FlipSky weather widget.

Twice a day a GitHub Action takes the latest 00 UTC run of the **CAMS European air quality forecasts**
(ensemble, surface, hourly, 0–72 h), cuts it into 1° tiles and publishes them on GitHub Pages
(`t/<lat>_<lon>.bin.gz`, format described in `cams_tiles.py`).

Contains modified Copernicus Atmosphere Monitoring Service information (CC-BY 4.0).
Neither the European Commission nor ECMWF is responsible for any use that may be made of the
Copernicus information or data it contains.
