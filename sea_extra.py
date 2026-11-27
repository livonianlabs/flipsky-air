"""Sea level and sea ice for FlipSky's map, from Copernicus Marine, cut into tiles the app downloads.

usage: sea_extra.py <out dir> [--only level|ice]
  Needs CMEMS_USER and CMEMS_PASS (a free Copernicus Marine account); without them it says so and does nothing.

Sea level, the Baltic's (storm surges, seiches), every 3 hours for 5 days:
  cmems_mod_bal_phy_anfc_PT1H-i, sla (above the mean) -> <out>/lb/<lat>_<lon>.bin.gz, 1° tiles at about 1 nautical mile.
Sea ice, daily, concentration and thickness:
  the world's polar seas: cmems_mod_glo_phy_anfc_0.083deg_P1D-m, averaged to 0.25°, 7 days -> <out>/iw/, 5° tiles;
  the Baltic: cmems_mod_bal_phy_anfc_P1D-m, 5 days -> <out>/ib/, 1° tiles.

Tiles in sea_motion.py's 'SM' format:
  level: one field, metres × 50 + 127 (2 cm steps, ±2.54 m);
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
    return np.where(np.isfinite(m), np.clip(np.round(m * 50) + 127, 0, 254), 255).astype(np.uint8)


def q_ice(c, h):
    ice = np.isfinite(c) & (c >= 0.01)
    qc = np.where(ice, np.clip(np.round(c * 250), 0, 254), 255).astype(np.uint8)
    qh = np.where(ice & np.isfinite(h), np.clip(np.round(np.nan_to_num(h) * 50), 0, 254), 255).astype(np.uint8)
    return qc, qh


def unix_h(t):
    return int(t.replace(tzinfo=dt.timezone.utc).timestamp() // 3600)


def level(out, today):
    t0 = dt.datetime.combine(today, dt.time())
    frames, lat, lon = [], None, None
    for i in range(5 * 8):
        when = t0 + dt.timedelta(hours=3 * i)
        try:
            _, lat, lon, v = read(subset(LEVEL_BAL, ['sla'], when, when), ['sla'])
        except Exception as e:
            if not frames:
                raise
            print('lb ends after', i, 'steps:', str(e)[:120])
            break
        frames.append(v['sla'][0])
    a = np.stack(frames)
    n = write_tiles(os.path.join(out, 'lb'), 1, unix_h(t0), 3, lat, lon, [q_level(a)])
    print('lb', n, a.shape, float(np.nanmin(a)), float(np.nanmax(a)))
    return [{'dataset': LEVEL_BAL, 'tiles': n, 'shape': list(a.shape)}]


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
