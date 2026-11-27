"""The sea in motion for FlipSky's map, from Copernicus Marine: currents and waves, cut into tiles the app downloads.

usage: sea_motion.py <out dir> [--only cur|waves]
  Needs CMEMS_USER and CMEMS_PASS (a free Copernicus Marine account); without them it says so and does nothing.

Currents, the surface's daily mean (u east, v north), 10 days:
  the world: cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m, averaged to 0.25° -> <out>/cu/<lat>_<lon>.bin.gz, 5° tiles;
  Europe's seas at 1.5-4 km (the Baltic, the North-West shelf, Iberia-Biscay, the Black Sea, the Mediterranean)
  -> <out>/cr/<lat>_<lon>.bin.gz, 1° tiles, the first model to give a tile its sea keeps it.
Waves, the significant wave height and the direction they come from, every 6 hours:
  the world: cmems_mod_glo_wav_anfc_0.083deg_PT3H-i, averaged to 0.25°, 7 days -> <out>/wm/, 5° tiles;
  Europe's seas, their hourly models taken every 6 hours, 5 days -> <out>/wr/, 1° tiles, first model keeps a tile.

Tile format (little-endian), gzipped:
  'SM' 1 (version), int32 t0 (Unix hours of step 0), uint8 steps, uint8 hours between steps, uint8 rows, uint8 cols,
  float32 lat of row 0, float32 lon of col 0, float32 lat step, float32 lon step (rows go north, cols east),
  uint8 fields, then uint8 [fields][steps][rows][cols], 255 = land or no data:
    currents: u, v as (m/s / 0.02) + 127 (±2.54 m/s);
    waves: height as metres × 20 (0.05 m steps, to 12.7 m), direction it comes from as degrees / 360 × 254.
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

CUR_WORLD = 'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m'
CUR_SEAS = [
    'cmems_mod_bal_phy_anfc_P1D-m',
    'cmems_mod_nws_phy-cur_anfc_1.5km-3D_P1D-m',
    'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m',
    'cmems_mod_blk_phy-cur_anfc_2.5km_P1D-m',
    'cmems_mod_med_phy-cur_anfc_4.2km_P1D-m',
]
WAV_WORLD = 'cmems_mod_glo_wav_anfc_0.083deg_PT3H-i'
WAV_SEAS = [
    'cmems_mod_bal_wav_anfc_PT1H-i',
    'cmems_mod_nws_wav_anfc_1.5km_PT1H-i',
    'cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i',
    'cmems_mod_blk_wav_anfc_2.5km_PT1H-i',
    'cmems_mod_med_wav_anfc_4.2km_PT1H-i',
]


def subset(dataset, variables, start, end, depth=False):
    """One download to a temporary netCDF file; its path."""
    import copernicusmarine
    nc = os.path.join(tempfile.mkdtemp(), 'x.nc')
    kw = dict(minimum_depth=0, maximum_depth=1.5) if depth else {}
    copernicusmarine.subset(
        dataset_id=dataset, variables=variables,
        start_datetime=start.strftime('%Y-%m-%dT%H:%M:%S'), end_datetime=end.strftime('%Y-%m-%dT%H:%M:%S'),
        username=os.environ['CMEMS_USER'], password=os.environ['CMEMS_PASS'],
        output_directory=os.path.dirname(nc), output_filename=os.path.basename(nc), overwrite=True, **kw,
    )
    return nc


def read(nc, names):
    """(times as Unix hours, lat ascending, lon -180..180 ascending, {name: [time, lat, lon] float32 with NaN})."""
    import netCDF4
    ds = netCDF4.Dataset(nc)
    t = ds.variables['time']
    import calendar
    times = [calendar.timegm(x.timetuple()) // 3600 for x in netCDF4.num2date(t[:], t.units)]
    lat = np.array(ds.variables['latitude'][:], dtype=np.float64)
    lon = np.array(ds.variables['longitude'][:], dtype=np.float64)
    lon = np.where(lon > 180, lon - 360, lon)
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    out = {}
    for n in names:
        a = np.ma.filled(np.ma.masked_invalid(np.ma.array(ds.variables[n][:])), np.nan).astype(np.float32)
        if a.ndim == 4:
            a = a[:, 0]
        out[n] = a[:, la_o][:, :, lo_o]
    ds.close()
    return times, lat[la_o], lon[lo_o], out


def coarsen(a, k):
    """Mean over k×k blocks of the last two axes, ignoring NaN (all-NaN blocks stay NaN)."""
    t, h, w = a.shape
    h2, w2 = h // k * k, w // k * k
    b = a[:, :h2, :w2].reshape(t, h2 // k, k, w2 // k, k)
    with np.errstate(all='ignore'):
        return np.nanmean(b, axis=(2, 4)).astype(np.float32)


def coarsen_axis(v, k):
    return v[:len(v) // k * k].reshape(-1, k).mean(axis=1)


def write_tiles(folder, size, t0, step_h, lat, lon, fields):
    """fields: list of uint8 arrays [steps, lat, lon]; tiles of [size] degrees, an existing tile is kept."""
    os.makedirs(folder, exist_ok=True)
    lat_step = float(round(abs(lat[1] - lat[0]), 6))
    lon_step = float(round(abs(lon[1] - lon[0]), 6))
    steps = fields[0].shape[0]
    count = 0
    for la in range(int(np.floor(lat[0] / size) * size), int(np.ceil(lat[-1])), size):
        rows = np.nonzero((lat >= la) & (lat < la + size))[0]
        if len(rows) == 0:
            continue
        for lo in range(int(np.floor(lon[0] / size) * size), int(np.ceil(lon[-1])), size):
            cols = np.nonzero((lon >= lo) & (lon < lo + size))[0]
            if len(cols) == 0:
                continue
            blocks = [f[:, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1] for f in fields]
            if (blocks[0] == 255).all():
                continue
            target = os.path.join(folder, f'{la}_{lo}.bin.gz')
            if os.path.exists(target):
                continue
            assert blocks[0].shape[1] < 256 and blocks[0].shape[2] < 256, blocks[0].shape
            head = b'SM' + struct.pack('<BiBBBBffffB', 1, t0, steps, step_h, blocks[0].shape[1], blocks[0].shape[2],
                                       float(lat[rows[0]]), float(lon[cols[0]]), lat_step, lon_step, len(fields))
            with gzip.open(target, 'wb', compresslevel=9) as f:
                f.write(head + b''.join(b.tobytes() for b in blocks))
            count += 1
    return count


def q_uv(a):
    return np.where(np.isfinite(a), np.clip(np.round(a / 0.02) + 127, 0, 254), 255).astype(np.uint8)


def q_height(h):
    return np.where(np.isfinite(h), np.clip(np.round(h * 20), 0, 254), 255).astype(np.uint8)


def q_dir(d, h):
    return np.where(np.isfinite(d) & np.isfinite(h), np.clip(np.round((d % 360) / 360 * 254), 0, 254), 255).astype(np.uint8)


def currents(out, today):
    summary = []
    start = dt.datetime.combine(today, dt.time())
    end = start + dt.timedelta(days=9)
    # The world, averaged to a quarter of a degree (three cells each way).
    times, lat, lon, v = read(subset(CUR_WORLD, ['uo', 'vo'], start, end, depth=True), ['uo', 'vo'])
    u, w = coarsen(v['uo'], 3), coarsen(v['vo'], 3)
    lat, lon = coarsen_axis(lat, 3), coarsen_axis(lon, 3)
    n = write_tiles(os.path.join(out, 'cu'), 5, times[0], 24, lat, lon, [q_uv(u), q_uv(w)])
    summary.append({'dataset': CUR_WORLD, 'tiles': n, 'shape': list(u.shape)})
    print('cu', n)
    del u, w, v
    for ds in CUR_SEAS:
        try:
            times, lat, lon, v = read(subset(ds, ['uo', 'vo'], start, end, depth=True), ['uo', 'vo'])
            n = write_tiles(os.path.join(out, 'cr'), 1, times[0], 24, lat, lon, [q_uv(v['uo']), q_uv(v['vo'])])
            summary.append({'dataset': ds, 'tiles': n, 'shape': list(v['uo'].shape)})
            print('cr', ds, n)
        except Exception as e:
            print('cr', ds, 'failed:', e)
    return summary


def waves_at(ds, when, k):
    """One moment's height and direction, averaged over k×k cells (the direction as a vector, by height)."""
    times, lat, lon, v = read(subset(ds, ['VHM0', 'VMDR'], when, when), ['VHM0', 'VMDR'])
    h, d = v['VHM0'][:1], v['VMDR'][:1]
    if k > 1:
        rad = np.radians(d)
        sx, cx = coarsen(np.sin(rad) * h, k), coarsen(np.cos(rad) * h, k)
        h = coarsen(h, k)
        d = (np.degrees(np.arctan2(sx, cx)) + 360) % 360
        lat, lon = coarsen_axis(lat, k), coarsen_axis(lon, k)
    return lat, lon, h[0], d[0]


def waves(out, today):
    summary = []
    t0 = dt.datetime.combine(today, dt.time())
    for ds, folder, size, k, days in [(WAV_WORLD, 'wm', 5, 3, 7)] + [(s, 'wr', 1, 1, 5) for s in WAV_SEAS]:
        try:
            hs, ds_ = [], []
            lat = lon = None
            for i in range(days * 4):
                try:
                    lat, lon, h, d = waves_at(ds, t0 + dt.timedelta(hours=6 * i), k)
                except Exception as e:
                    # A model whose forecast ends sooner: the steps it has.
                    if not hs:
                        raise
                    print(folder, ds, 'ends after', i, 'steps:', str(e)[:120])
                    break
                hs.append(h)
                ds_.append(d)
            h, d = np.stack(hs), np.stack(ds_)
            t0h = int(t0.replace(tzinfo=dt.timezone.utc).timestamp() // 3600)
            n = write_tiles(os.path.join(out, folder), size, t0h, 6, lat, lon, [q_height(h), q_dir(d, h)])
            summary.append({'dataset': ds, 'tiles': n, 'shape': list(h.shape)})
            print(folder, ds, n)
        except Exception as e:
            print(folder, ds, 'failed:', e)
    return summary


def main():
    args = sys.argv[1:]
    out = args[0]
    only = args[args.index('--only') + 1] if '--only' in args else None
    if not os.environ.get('CMEMS_USER') or not os.environ.get('CMEMS_PASS'):
        print('sea motion: no Copernicus Marine account (CMEMS_USER, CMEMS_PASS): skipped')
        return
    today = dt.datetime.now(dt.timezone.utc).date()
    summary = {}
    if only in (None, 'cur'):
        try:
            summary['currents'] = currents(out, today)
        except Exception as e:
            print('currents failed:', e)
    if only in (None, 'waves'):
        summary['waves'] = waves(out, today)
    with open(os.path.join(out, 'motion.json'), 'w') as f:
        json.dump({'date': today.isoformat(), 'sets': summary,
                   'source': 'Generated using E.U. Copernicus Marine Service Information'}, f, indent=1)


if __name__ == '__main__':
    main()
