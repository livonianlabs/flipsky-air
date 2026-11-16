"""Sea surface salinity for FlipSky's water card: a year's mean from Copernicus Marine, cut into tiles once a year.

usage: salinity.py <out dir> [--year YYYY]
  Needs CMEMS_USER and CMEMS_PASS. Takes the monthly means of the global analysis (cmems_mod_glo_phy-so_anfc_0.083deg_P1M-m,
  top layer) for the last full year, averages them, coarsens 3 × 3 to 0.25°, and writes <out>/<lat>_<lon>.bin.gz,
  10-degree tiles (lower-left corner, multiples of 10).

Tile format (little-endian), gzipped:
  'SL' 1, uint8 rows, uint8 cols, float32 lat of row 0, float32 lon of col 0, float32 step (degrees; rows go north,
  cols east), then uint8 [rows][cols]: salinity × 5 (0.2 ‰ steps up to 50.8 ‰), 255 = land or no data.
Salinity barely changes at a place from year to year, so this is run by hand (the salinity workflow) and committed.

Generated using E.U. Copernicus Marine Service Information.
"""
import datetime as dt
import gzip
import os
import struct
import sys
import tempfile

import numpy as np

DATASET = 'cmems_mod_glo_phy-so_anfc_0.083deg_P1M-m'


def main():
    args = sys.argv[1:]
    out = args[0]
    year = int(args[args.index('--year') + 1]) if '--year' in args else dt.date.today().year - 1
    import copernicusmarine
    import netCDF4
    nc = os.path.join(tempfile.mkdtemp(), 'so.nc')
    copernicusmarine.subset(
        dataset_id=DATASET, variables=['so'], minimum_depth=0, maximum_depth=1.5,
        start_datetime=f'{year}-01-01T00:00:00', end_datetime=f'{year}-12-31T00:00:00',
        username=os.environ['CMEMS_USER'], password=os.environ['CMEMS_PASS'],
        output_directory=os.path.dirname(nc), output_filename=os.path.basename(nc), overwrite=True,
    )
    ds = netCDF4.Dataset(nc)
    a = np.ma.filled(np.ma.masked_invalid(np.ma.array(ds.variables['so'][:])), np.nan).astype(np.float32)
    if a.ndim == 4:
        a = a[:, 0]
    with np.errstate(all='ignore'):
        mean = np.nanmean(a, axis=0)
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    lon = np.where(lon > 180, lon - 360, lon)
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    lat, lon, mean = lat[la_o], lon[lo_o], mean[la_o][:, lo_o]
    # 3 × 3 cells to one: 0.25°, the coast's land cells left out of the mean.
    r, c = (mean.shape[0] // 3) * 3, (mean.shape[1] // 3) * 3
    blocks = mean[:r, :c].reshape(r // 3, 3, c // 3, 3)
    with np.errstate(all='ignore'):
        coarse = np.nanmean(blocks, axis=(1, 3))
    clat = lat[:r].reshape(-1, 3).mean(axis=1)
    clon = lon[:c].reshape(-1, 3).mean(axis=1)
    step = float(round(clat[1] - clat[0], 6))
    q = np.where(np.isfinite(coarse), np.clip(np.round(coarse * 5), 0, 254), 255).astype(np.uint8)
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
            block = q[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            if (block == 255).all():
                continue
            head = b'SL' + struct.pack('<BBBfff', 1, block.shape[0], block.shape[1], float(clat[rows[0]]), float(clon[cols[0]]), step)
            with gzip.open(os.path.join(out, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.tobytes())
            count += 1
    print('salinity', year, 'tiles', count, 'grid', coarse.shape, 'step', step)


if __name__ == '__main__':
    main()
