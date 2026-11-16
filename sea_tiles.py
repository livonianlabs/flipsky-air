"""Sea surface temperature for FlipSky from Copernicus Marine: today and the days ahead, cut into tiles the app downloads.

usage: sea_tiles.py <out dir> [--days N]
  Needs CMEMS_USER and CMEMS_PASS (a free Copernicus Marine account); without them it says so and does nothing.
  Global: cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m (about 9 km, daily means, 10 days ahead) -> <out>/s/<lat>_<lon>.bin.gz,
    5-degree tiles (lower-left corner, multiples of 5).
  Baltic: cmems_mod_bal_phy_anfc_P1D-m (about 2 km) -> <out>/sb/<lat>_<lon>.bin.gz, 1-degree tiles; the app prefers them.

Tile format (little-endian), gzipped:
  'SW' 1 (version), int32 base (Unix days of day 0), uint8 days, uint8 rows, uint8 cols,
  float32 lat of row 0, float32 lon of col 0, float32 lat step, float32 lon step (degrees; rows go north, cols east),
  then uint8 [days][rows][cols]: (°C + 5) * 5, so 0.2° steps from -5 to 45.8; 255 = land or no data.
Tiles that are all land are not written.

Generated using E.U. Copernicus Marine Service Information.
"""
import datetime as dt
import gzip
import json
import os
import struct
import sys
import tempfile

import numpy as np

SETS = [
    # dataset, folder, tile size (degrees)
    ('cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m', 's', 5),
    ('cmems_mod_bal_phy_anfc_P1D-m', 'sb', 1),
]


def fetch(dataset, start, end, target):
    import copernicusmarine
    copernicusmarine.subset(
        dataset_id=dataset, variables=['thetao'], minimum_depth=0, maximum_depth=1.5,
        start_datetime=start.isoformat() + 'T00:00:00', end_datetime=end.isoformat() + 'T00:00:00',
        username=os.environ['CMEMS_USER'], password=os.environ['CMEMS_PASS'],
        output_directory=os.path.dirname(target), output_filename=os.path.basename(target), overwrite=True,
    )


def cut(path, folder, size, base):
    import netCDF4
    ds = netCDF4.Dataset(path)
    v = ds.variables['thetao']
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    a = np.ma.filled(np.ma.masked_invalid(np.ma.array(v[:])), np.nan).astype(np.float32)
    # time, depth, lat, lon -> days, lat, lon (the top layer).
    if a.ndim == 4:
        a = a[:, 0]
    lo_order = np.argsort(np.where(lon > 180, lon - 360, lon))
    lon = np.where(lon > 180, lon - 360, lon)[lo_order]
    la_order = np.argsort(lat)
    lat = lat[la_order]
    a = a[:, la_order][:, :, lo_order]
    q = np.where(np.isfinite(a), np.clip(np.round((a + 5) * 5), 0, 254), 255).astype(np.uint8)
    lat_step = float(round(abs(lat[1] - lat[0]), 6))
    lon_step = float(round(abs(lon[1] - lon[0]), 6))
    os.makedirs(folder, exist_ok=True)
    count = 0
    for la in range(int(np.floor(lat[0] / size) * size), int(np.ceil(lat[-1])), size):
        rows = np.nonzero((lat >= la) & (lat < la + size))[0]
        if len(rows) == 0:
            continue
        for lo in range(int(np.floor(lon[0] / size) * size), int(np.ceil(lon[-1])), size):
            cols = np.nonzero((lon >= lo) & (lon < lo + size))[0]
            if len(cols) == 0:
                continue
            block = q[:, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            if (block == 255).all():
                continue
            assert block.shape[1] < 256 and block.shape[2] < 256, block.shape
            head = b'SW' + struct.pack('<BiBBBffff', 1, base, block.shape[0], block.shape[1], block.shape[2],
                                       float(lat[rows[0]]), float(lon[cols[0]]), lat_step, lon_step)
            with gzip.open(os.path.join(folder, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.tobytes())
            count += 1
    return count, a.shape


def main():
    args = sys.argv[1:]
    out = args[0]
    days = int(args[args.index('--days') + 1]) if '--days' in args else 10
    if not os.environ.get('CMEMS_USER') or not os.environ.get('CMEMS_PASS'):
        print('sea: no Copernicus Marine account (CMEMS_USER, CMEMS_PASS): skipped')
        return
    today = dt.datetime.now(dt.timezone.utc).date()
    base = (today - dt.date(1970, 1, 1)).days
    summary = {}
    for dataset, folder, size in SETS:
        try:
            nc = os.path.join(tempfile.mkdtemp(), folder + '.nc')
            fetch(dataset, today, today + dt.timedelta(days=days - 1), nc)
            n, shape = cut(nc, os.path.join(out, folder), size, base)
            summary[folder] = {'dataset': dataset, 'tiles': n, 'shape': list(shape)}
            print(folder, n, 'tiles', shape)
        except Exception as e:
            print(folder, 'failed:', e)
    with open(os.path.join(out, 'sea.json'), 'w') as f:
        json.dump({'base': base, 'days': days, 'sets': summary,
                   'source': 'Generated using E.U. Copernicus Marine Service Information'}, f, indent=1)


if __name__ == '__main__':
    main()
