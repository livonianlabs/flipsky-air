"""Air quality for FlipSky from Copernicus CAMS: the European forecast cut into small tiles the app downloads.

usage: cams_tiles.py <out dir> [--date YYYY-MM-DD] [--nc file.nc]
  Fetches today's 00 UTC run (yesterday's if today's isn't out yet) of "CAMS European air quality forecasts" (ensemble, surface, hourly, 0-72 h)
  from the Atmosphere Data Store (key in ~/.cdsapirc or ADS_KEY / ADS_URL), or reads --nc.
  Writes <out>/t/<lat>_<lon>.bin.gz, one per 1-degree square (lower-left corner, whole degrees), and
  <out>/index.json.

Tile format (little-endian), gzipped:
  'AQ' 1 (version), int32 base (Unix hours of lead time 0), uint8 hours, uint8 rows, uint8 cols,
  float32 lat of row 0, float32 lon of col 0, float32 step (degrees; rows go north, cols east),
  then uint16 [hours][5][rows][cols]: PM2.5, PM10, O3, NO2, SO2 in 0.1 ug/m3, 65535 = no data.

Contains modified Copernicus Atmosphere Monitoring Service information.
"""
import datetime as dt
import gzip
import json
import os
import struct
import sys
import tempfile
import zipfile

import numpy as np

DATASET = 'cams-europe-air-quality-forecasts'
# The app's order; the names CAMS gives them in the request and in the file.
POLLUTANTS = [
    ('particulate_matter_2.5um', 'pm2p5'),
    ('particulate_matter_10um', 'pm10'),
    ('ozone', 'o3'),
    ('nitrogen_dioxide', 'no2'),
    ('sulphur_dioxide', 'so2'),
]
# Three days ahead: if a daily run fails, the app still has today and tomorrow.
HOURS = 73


def fetch(date, target):
    import cdsapi
    kw = {}
    if os.environ.get('ADS_KEY'):
        kw = {'url': os.environ.get('ADS_URL', 'https://ads.atmosphere.copernicus.eu/api'), 'key': os.environ['ADS_KEY']}
    client = cdsapi.Client(**kw)
    client.retrieve(DATASET, {
        'variable': [p for p, _ in POLLUTANTS],
        'model': ['ensemble'],
        'level': ['0'],
        'date': [f'{date}/{date}'],
        'type': ['forecast'],
        'time': ['00:00'],
        'leadtime_hour': [str(h) for h in range(HOURS)],
        'data_format': 'netcdf_zip',
    }, target)


def open_nc(path):
    """The netCDF file, from the zip the ADS hands out or as is."""
    import netCDF4
    if zipfile.is_zipfile(path):
        tmp = tempfile.mkdtemp()
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.endswith('.nc')]
            z.extractall(tmp, names)
        return [netCDF4.Dataset(os.path.join(tmp, n)) for n in names]
    return [netCDF4.Dataset(path)]


def var_named(datasets, short):
    for ds in datasets:
        for name, v in ds.variables.items():
            if name.startswith(short) and v.ndim >= 3:
                return ds, v
    raise KeyError(short)


def main():
    args = sys.argv[1:]
    out = args[0]
    date = args[args.index('--date') + 1] if '--date' in args else dt.datetime.now(dt.timezone.utc).date().isoformat()
    nc = args[args.index('--nc') + 1] if '--nc' in args else None
    if nc is None:
        nc = os.path.join(tempfile.mkdtemp(), 'cams.zip')
        try:
            fetch(date, nc)
        except Exception as e:
            if '--date' in args:
                raise
            print('today not available:', e)
            date = (dt.date.fromisoformat(date) - dt.timedelta(days=1)).isoformat()
            fetch(date, nc)
    datasets = open_nc(nc)

    ds, first = var_named(datasets, POLLUTANTS[0][1])
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    lon = np.where(lon > 180, lon - 360, lon)
    # Rows from south to north, columns from west to east.
    lat_order = np.argsort(lat)
    lon_order = np.argsort(lon)
    lat, lon = lat[lat_order], lon[lon_order]
    step = float(round(abs(lat[1] - lat[0]), 4))

    cube = np.full((HOURS, len(POLLUTANTS), len(lat), len(lon)), 65535, dtype=np.uint16)
    for k, (_, short) in enumerate(POLLUTANTS):
        _, v = var_named(datasets, short)
        a = np.ma.filled(np.ma.masked_invalid(np.squeeze(np.ma.array(v[:]))), np.nan).astype(np.float32)
        a = a[:HOURS][:, lat_order][:, :, lon_order]
        q = np.where(np.isfinite(a), np.clip(np.round(a * 10), 0, 65534), 65535).astype(np.uint16)
        cube[:q.shape[0], k] = q
    base = int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp() // 3600)

    os.makedirs(os.path.join(out, 't'), exist_ok=True)
    count = 0
    for la in range(int(np.floor(lat[0])), int(np.ceil(lat[-1]))):
        rows = np.nonzero((lat >= la) & (lat < la + 1))[0]
        if len(rows) == 0:
            continue
        for lo in range(int(np.floor(lon[0])), int(np.ceil(lon[-1]))):
            cols = np.nonzero((lon >= lo) & (lon < lo + 1))[0]
            if len(cols) == 0:
                continue
            block = cube[:, :, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            head = b'AQ' + struct.pack('<BiBBBfff', 1, base, HOURS, len(rows), len(cols), float(lat[rows[0]]), float(lon[cols[0]]), step)
            with gzip.open(os.path.join(out, 't', f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.astype('<u2').tobytes())
            count += 1
    with open(os.path.join(out, 'index.json'), 'w') as f:
        json.dump({
            'base': base, 'hours': HOURS, 'south': float(lat[0]), 'north': float(lat[-1]), 'west': float(lon[0]), 'east': float(lon[-1]),
            'step': step, 'tiles': count, 'source': 'Contains modified Copernicus Atmosphere Monitoring Service information ' + date[:4],
        }, f, indent=1)
    print('tiles', count, 'grid', len(lat), 'x', len(lon), 'step', step, 'base', base)


if __name__ == '__main__':
    main()
