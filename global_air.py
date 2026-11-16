"""Air quality for the rest of the world from the CAMS global composition forecast, beside the European tiles.

usage: global_air.py <out dir> [--date YYYY-MM-DD] [--nc file.zip]
  Today's 00 UTC run (yesterday's if not out yet) of "CAMS global atmospheric composition forecasts" (about 0.4°,
  every 3 hours, 0-72 h): PM2.5 and PM10 at the surface, ozone, NO2 and SO2 on the lowest model level (137).
  Writes <out>/ga/<lat>_<lon>.bin.gz, 5-degree tiles (lower-left corner, multiples of 5); the app uses them where the
  European tiles have nothing (outside Europe).

Tile format (little-endian), gzipped:
  'AQ' 2 (version), int32 base (Unix hours of the first step), uint8 steps, uint8 every (hours between steps),
  uint8 rows, uint8 cols, float32 lat of row 0, float32 lon of col 0, float32 step (degrees; rows go north, cols east),
  then uint16 [steps][5][rows][cols]: PM2.5, PM10, O3, NO2, SO2 in 0.1 ug/m3, 65535 = no data.
Gases come as mass mixing ratios (kg/kg) and are turned into ug/m3 with the air's density near the ground (1.2 kg/m3).

Contains modified Copernicus Atmosphere Monitoring Service information.
"""
import datetime as dt
import gzip
import os
import struct
import sys
import tempfile
import zipfile

import numpy as np

DATASET = 'cams-global-atmospheric-composition-forecasts'
VARS = ['particulate_matter_2.5um', 'particulate_matter_10um', 'ozone', 'nitrogen_dioxide', 'sulphur_dioxide']
# The file's names, the app's order, and the factor to ug/m3.
SHORT = [('pm2p5', 1e9), ('pm10', 1e9), ('go3', 1.2e9), ('no2', 1.2e9), ('so2', 1.2e9)]
EVERY = 3
STEPS = 25  # 0-72 h


def fetch(date, target):
    import cdsapi
    kw = {}
    if os.environ.get('ADS_KEY'):
        kw = {'url': os.environ.get('ADS_URL', 'https://ads.atmosphere.copernicus.eu/api'), 'key': os.environ['ADS_KEY']}
    cdsapi.Client(**kw).retrieve(DATASET, {
        'variable': VARS, 'model_level': ['137'],
        'date': [f'{date}/{date}'], 'type': ['forecast'], 'time': ['00:00'],
        'leadtime_hour': [str(h) for h in range(0, STEPS * EVERY, EVERY)],
        'data_format': 'netcdf_zip',
    }, target)


def main():
    args = sys.argv[1:]
    out = args[0]
    date = args[args.index('--date') + 1] if '--date' in args else dt.datetime.now(dt.timezone.utc).date().isoformat()
    path = args[args.index('--nc') + 1] if '--nc' in args else None
    if path is None:
        path = os.path.join(tempfile.mkdtemp(), 'global.zip')
        try:
            fetch(date, path)
        except Exception as e:
            if '--date' in args:
                raise
            print('today not available:', e)
            date = (dt.date.fromisoformat(date) - dt.timedelta(days=1)).isoformat()
            fetch(date, path)
    import netCDF4
    d = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as z:
        z.extractall(d)
        datasets = [netCDF4.Dataset(os.path.join(d, n)) for n in z.namelist() if n.endswith('.nc')]

    def var(short):
        for ds in datasets:
            if short in ds.variables:
                return ds, ds.variables[short]
        raise KeyError(short)

    ds, _ = var('pm2p5')
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    lon = np.where(lon > 180, lon - 360, lon)
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    lat, lon = lat[la_o], lon[lo_o]
    step = float(round(abs(lat[1] - lat[0]), 4))
    cube = np.full((STEPS, len(SHORT), len(lat), len(lon)), 65535, dtype=np.uint16)
    for k, (short, factor) in enumerate(SHORT):
        _, v = var(short)
        a = np.ma.filled(np.ma.masked_invalid(np.squeeze(np.ma.array(v[:]))), np.nan).astype(np.float64) * factor
        a = a[:STEPS][:, la_o][:, :, lo_o]
        cube[:a.shape[0], k] = np.where(np.isfinite(a), np.clip(np.round(a * 10), 0, 65534), 65535).astype(np.uint16)
    base = int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp() // 3600)

    folder = os.path.join(out, 'ga')
    os.makedirs(folder, exist_ok=True)
    count = 0
    for la in range(-90, 90, 5):
        rows = np.nonzero((lat >= la) & (lat < la + 5))[0]
        if len(rows) == 0:
            continue
        for lo in range(-180, 180, 5):
            cols = np.nonzero((lon >= lo) & (lon < lo + 5))[0]
            if len(cols) == 0:
                continue
            block = cube[:, :, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            head = b'AQ' + struct.pack('<BiBBBBfff', 2, base, STEPS, EVERY, len(rows), len(cols), float(lat[rows[0]]), float(lon[cols[0]]), step)
            with gzip.open(os.path.join(folder, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.astype('<u2').tobytes())
            count += 1
    print('global air tiles', count, 'grid', len(lat), 'x', len(lon), 'step', step, 'base', base)


if __name__ == '__main__':
    main()
