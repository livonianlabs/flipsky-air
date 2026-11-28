"""Sea level and sea ice for FlipSky's map, from Copernicus Marine, cut into tiles the app downloads.

usage: sea_extra.py <out dir> [--only level|ice]
  Needs CMEMS_USER and CMEMS_PASS (a free Copernicus Marine account); without them it says so and does nothing.

Sea level, every 3 hours for 5 days, Europe's seas in 1° tiles (the first model to give a tile its sea keeps it):
  the Baltic's (storm surges, seiches): cmems_mod_bal_phy_anfc_PT1H-i, sla (above the mean), about 1 nautical mile;
  the North-West shelf (tides), Iberia-Biscay, the Black Sea: their hourly zos, less the model's mean over the
  5 days, so all read from about 0 -> <out>/lb/<lat>_<lon>.bin.gz.
  The Mediterranean: its model gives sea level only as daily means (cmems_mod_med_phy-ssh_anfc_4.2km_P1D-m, zos, its
  mean over the days taken off as well): the same tiles, a step a day, where the others have no sea.
Sea ice, daily, concentration and thickness:
  the world's polar seas: cmems_mod_glo_phy_anfc_0.083deg_P1D-m, averaged to 0.25°, 7 days -> <out>/iw/, 5° tiles;
  the Baltic: cmems_mod_bal_phy_anfc_P1D-m, 5 days -> <out>/ib/, 1° tiles.

Tiles in sea_motion.py's 'SM' format:
  level: one field, metres × 25 + 127 (4 cm steps, ±5.08 m: the Channel's tides);
  ice: concentration × 250 (0..1), thickness in metres × 50 (to 5.08 m); open water and land are 255 (no tile when a
  tile has no ice at all).

Generated using E.U. Copernicus Marine Service Information.
"""
import datetime as dt
import json
import os
import sys
import tempfile

import numpy as np

from sea_motion import coarsen, coarsen_axis, read, write_tiles

LEVEL_BAL = 'cmems_mod_bal_phy_anfc_PT1H-i'
LEVEL_SEAS = [
    'cmems_mod_nws_phy-ssh_anfc_1.5km-2D_PT1H-i',
    'cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m',
    'cmems_mod_blk_phy-ssh_anfc_2.5km_PT1H-m',
]
LEVEL_MED = 'cmems_mod_med_phy-ssh_anfc_4.2km_P1D-m'
ICE_WORLD = 'cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
ICE_BAL = 'cmems_mod_bal_phy_anfc_P1D-m'


def subset(dataset, variables, start, end, south=None, north=None):
    import copernicusmarine
    nc = os.path.join(tempfile.mkdtemp(), 'x.nc')
    kw = {}
    if south is not None:
        kw.update(minimum_latitude=south, maximum_latitude=north)
    copernicusmarine.subset(
        dataset_id=dataset, variables=variables,
        start_datetime=start.strftime('%Y-%m-%dT%H:%M:%S'), end_datetime=end.strftime('%Y-%m-%dT%H:%M:%S'),
        username=os.environ['CMEMS_USER'], password=os.environ['CMEMS_PASS'],
        output_directory=os.path.dirname(nc), output_filename=os.path.basename(nc), overwrite=True, **kw,
    )
    return nc


def q_level(m):
    return np.where(np.isfinite(m), np.clip(np.round(m * 25) + 127, 0, 254), 255).astype(np.uint8)


def q_ice(c, h):
    ice = np.isfinite(c) & (c >= 0.01)
    qc = np.where(ice, np.clip(np.round(c * 250), 0, 254), 255).astype(np.uint8)
    qh = np.where(ice & np.isfinite(h), np.clip(np.round(np.nan_to_num(h) * 50), 0, 254), 255).astype(np.uint8)
    return qc, qh


def unix_h(t):
    return int(t.replace(tzinfo=dt.timezone.utc).timestamp() // 3600)


def level_of(ds, var, t0):
    """One model's sea level every 3 hours for 5 days (fewer if its forecast ends sooner)."""
    frames, lat, lon = [], None, None
    for i in range(5 * 8):
        when = t0 + dt.timedelta(hours=3 * i)
        try:
            _, lat, lon, v = read(subset(ds, [var], when, when), [var])
        except Exception as e:
            if not frames:
                raise
            print('lb', ds, 'ends after', i, 'steps:', str(e)[:120])
            break
        frames.append(v[var][0])
    return lat, lon, np.stack(frames)


def level(out, today):
    t0 = dt.datetime.combine(today, dt.time())
    summary = []
    for ds, var in [(LEVEL_BAL, 'sla')] + [(s, 'zos') for s in LEVEL_SEAS]:
        try:
            lat, lon, a = level_of(ds, var, t0)
            if var == 'zos':
                # The model's own datum: its mean over these days taken off, as the Baltic's sla is.
                a = a - np.nanmean(a)
            n = write_tiles(os.path.join(out, 'lb'), 1, unix_h(t0), 3, lat, lon, [q_level(a)])
            print('lb', ds, n, a.shape, float(np.nanmin(a)), float(np.nanmax(a)))
            summary.append({'dataset': ds, 'tiles': n, 'shape': list(a.shape)})
        except Exception as e:
            print('lb', ds, 'failed:', e)
    # The Mediterranean, daily (written last: where an hourly model has a tile, it keeps it).
    try:
        times, lat, lon, v = read(subset(LEVEL_MED, ['zos'], t0, t0 + dt.timedelta(days=5)), ['zos'])
        a = v['zos'] - np.nanmean(v['zos'])
        n = write_tiles(os.path.join(out, 'lb'), 1, times[0], 24, lat, lon, [q_level(a)])
        print('lb', LEVEL_MED, n, a.shape, float(np.nanmin(a)), float(np.nanmax(a)))
        summary.append({'dataset': LEVEL_MED, 'tiles': n, 'shape': list(a.shape)})
    except Exception as e:
        print('lb', LEVEL_MED, 'failed:', e)
    return summary


def ice(out, today):
    summary = []
    start = dt.datetime.combine(today, dt.time())
    # The world's polar seas only (the rest never freezes), each a download of its own.
    for south, north in ((40, 90), (-80, -45)):
        try:
            times, lat, lon, v = read(subset(ICE_WORLD, ['siconc', 'sithick'], start, start + dt.timedelta(days=6),
                                             south, north), ['siconc', 'sithick'])
            c, h = coarsen(v['siconc'], 3), coarsen(v['sithick'], 3)
            lat, lon = coarsen_axis(lat, 3), coarsen_axis(lon, 3)
            qc, qh = q_ice(c, h)
            n = write_tiles(os.path.join(out, 'iw'), 5, times[0], 24, lat, lon, [qc, qh])
            summary.append({'dataset': ICE_WORLD, 'lat': [south, north], 'tiles': n, 'shape': list(c.shape)})
            print('iw', south, north, n)
        except Exception as e:
            print('iw', south, north, 'failed:', e)
    try:
        times, lat, lon, v = read(subset(ICE_BAL, ['siconc', 'sithick'], start, start + dt.timedelta(days=4)),
                                  ['siconc', 'sithick'])
        qc, qh = q_ice(v['siconc'], v['sithick'])
        n = write_tiles(os.path.join(out, 'ib'), 1, times[0], 24, lat, lon, [qc, qh])
        summary.append({'dataset': ICE_BAL, 'tiles': n, 'shape': list(qc.shape)})
        print('ib', n)
    except Exception as e:
        print('ib failed:', e)
    return summary


def main():
    args = sys.argv[1:]
    out = args[0]
    only = args[args.index('--only') + 1] if '--only' in args else None
    if not os.environ.get('CMEMS_USER') or not os.environ.get('CMEMS_PASS'):
        print('sea extra: no Copernicus Marine account (CMEMS_USER, CMEMS_PASS): skipped')
        return
    today = dt.datetime.now(dt.timezone.utc).date()
    summary = {}
    if only in (None, 'level'):
        try:
            summary['level'] = level(out, today)
        except Exception as e:
            print('level failed:', e)
    if only in (None, 'ice'):
        summary['ice'] = ice(out, today)
    with open(os.path.join(out, 'extra.json'), 'w') as f:
        json.dump({'date': today.isoformat(), 'sets': summary,
                   'source': 'Generated using E.U. Copernicus Marine Service Information'}, f, indent=1)


if __name__ == '__main__':
    main()
