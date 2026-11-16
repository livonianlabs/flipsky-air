"""The sea's usual temperature by month, for trips planned further ahead than the forecast: ten years of Copernicus
Marine's global reanalysis, averaged month by month, cut into tiles once a year.

usage: sea_climate.py <out dir> [--from YYYY] [--to YYYY]
  Needs CMEMS_USER and CMEMS_PASS. Monthly means of cmems_mod_glo_phy_my_0.083deg_P1M-m (top layer), the years given
  (default 2012-2021), averaged per calendar month, coarsened 3 × 3 to 0.25°, written as <out>/<lat>_<lon>.bin.gz in
  10-degree tiles.

Tile format (little-endian), gzipped:
  'SC' 1, uint8 rows, uint8 cols, float32 lat of row 0, float32 lon of col 0, float32 step (rows go north, cols east),
  then uint8 [12][rows][cols]: (°C + 5) × 5, January first; 255 = land or no data.

Generated using E.U. Copernicus Marine Service Information.
"""
import gzip
import os
import struct
import sys
import tempfile

import numpy as np

DATASET = 'cmems_mod_glo_phy_my_0.083deg_P1M-m'


def main():
    args = sys.argv[1:]
    out = args[0]
    y0 = int(args[args.index('--from') + 1]) if '--from' in args else 2012
    y1 = int(args[args.index('--to') + 1]) if '--to' in args else 2021
    import copernicusmarine
    import netCDF4
    nc = os.path.join(tempfile.mkdtemp(), 'sst.nc')
    copernicusmarine.subset(
        dataset_id=DATASET, variables=['thetao'], minimum_depth=0, maximum_depth=1,
        start_datetime=f'{y0}-01-01T00:00:00', end_datetime=f'{y1}-12-31T00:00:00',
        username=os.environ['CMEMS_USER'], password=os.environ['CMEMS_PASS'],
        output_directory=os.path.dirname(nc), output_filename=os.path.basename(nc), overwrite=True,
    )
    ds = netCDF4.Dataset(nc)
    v = ds.variables['thetao']
    times = netCDF4.num2date(ds.variables['time'][:], ds.variables['time'].units)
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    lon = np.where(lon > 180, lon - 360, lon)
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    lat, lon = lat[la_o], lon[lo_o]
    r, c = (len(lat) // 3) * 3, (len(lon) // 3) * 3
    months = []
    for m in range(1, 13):
        idx = [i for i, t in enumerate(times) if t.month == m]
        a = np.ma.filled(np.ma.masked_invalid(np.ma.array(v[idx])), np.nan).astype(np.float32)
        if a.ndim == 4:
            a = a[:, 0]
        with np.errstate(all='ignore'):
            mean = np.nanmean(a, axis=0)[la_o][:, lo_o]
            coarse = np.nanmean(mean[:r, :c].reshape(r // 3, 3, c // 3, 3), axis=(1, 3))
        months.append(coarse)
        print('month', m, len(idx))
    cube = np.stack(months)
    clat = lat[:r].reshape(-1, 3).mean(axis=1)
    clon = lon[:c].reshape(-1, 3).mean(axis=1)
    step = float(round(clat[1] - clat[0], 6))
    q = np.where(np.isfinite(cube), np.clip(np.round((cube + 5) * 5), 0, 254), 255).astype(np.uint8)
    os.makedirs(out, exist_ok=True)
    count = 0
    for la in range(-90, 90, 10):
        rows = np.nonzero((clat >= la) & (clat < la + 10))[0]
        if len(rows) == 0:
            continue
        for lo in range(-180, 180, 10):
            cols = np.nonzero((clon >= lo) & (clon < lo + 10))[0]
            if len(cols) == 0:
                continue
            block = q[:, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            if (block == 255).all():
                continue
            head = b'SC' + struct.pack('<BBBfff', 1, block.shape[1], block.shape[2], float(clat[rows[0]]), float(clon[cols[0]]), step)
            with gzip.open(os.path.join(out, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.tobytes())
            count += 1
    print('sea climate', y0, y1, 'tiles', count)


if __name__ == '__main__':
    main()
