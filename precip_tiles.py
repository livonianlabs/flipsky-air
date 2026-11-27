"""Precipitation for FlipSky's map, finer than MET Norway's points: how much falls each hour and the air's temperature
(rain or snow), cut into tiles the app downloads.

usage: precip_tiles.py <out dir> [--only eu|world]
  Europe: DWD's ICON-EU (opendata.dwd.de, 0.0625°, about 7 km), the newest run out, hourly for 72 hours
    -> <out>/pe/<lat>_<lon>.bin.gz, 5° tiles.
  The world: ECMWF's open data IFS (data.ecmwf.int, 0.25°, CC BY 4.0), the newest 00 UTC run, every 3 hours for 6 days
    -> <out>/pw/<lat>_<lon>.bin.gz, 5° tiles.

Tiles in sea_motion.py's 'SM' format, two fields:
  the hour's precipitation as round(sqrt(mm/h) × 40) (0.1 mm/h -> 13, 1 -> 40, 10 -> 126, up to about 40 mm/h);
  the air's temperature at 2 m as round((°C + 60) × 2).
A step holds what fell in the hours up to it (an hour for ICON-EU, three for IFS), as an hourly rate.

Sources: Deutscher Wetterdienst (DWD) open data; ECMWF open data, CC BY 4.0.
"""
import bz2
import datetime as dt
import json
import os
import sys
import urllib.request

import numpy as np

from sea_motion import write_tiles

UA = {'User-Agent': 'FlipSky data job (github.com/livonianlabs)'}
DWD = 'https://opendata.dwd.de/weather/nwp/icon-eu/grib'
ECMWF = 'https://data.ecmwf.int/forecasts'


def get(url, rng=None):
    h = dict(UA)
    if rng:
        h['Range'] = 'bytes=%d-%d' % rng
    return urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=300).read()


def decode(msg):
    """(values [lat south->north, lon west->east], lat, lon) of one GRIB message."""
    import eccodes
    h = eccodes.codes_new_from_message(msg)
    try:
        nj, ni = eccodes.codes_get(h, 'Nj'), eccodes.codes_get(h, 'Ni')
        lat0 = eccodes.codes_get(h, 'latitudeOfFirstGridPointInDegrees')
        lon0 = eccodes.codes_get(h, 'longitudeOfFirstGridPointInDegrees')
        dj = eccodes.codes_get(h, 'jDirectionIncrementInDegrees')
        di = eccodes.codes_get(h, 'iDirectionIncrementInDegrees')
        up = eccodes.codes_get(h, 'jScansPositively') == 1
        v = eccodes.codes_get_values(h).reshape(nj, ni).astype(np.float32)
        missing = eccodes.codes_get(h, 'missingValue')
    finally:
        eccodes.codes_release(h)
    v[v == missing] = np.nan
    lat = lat0 + (1 if up else -1) * np.arange(nj) * dj
    lon = (lon0 + np.arange(ni) * di + 180) % 360 - 180
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    return v[la_o][:, lo_o], lat[la_o], lon[lo_o]


def q_rate(mm_h):
    return np.where(np.isfinite(mm_h), np.clip(np.round(np.sqrt(np.clip(mm_h, 0, None)) * 40), 0, 254), 255).astype(np.uint8)


def q_temp(c):
    return np.where(np.isfinite(c), np.clip(np.round((c + 60) * 2), 0, 254), 255).astype(np.uint8)


def unix_h(t):
    return int(t.replace(tzinfo=dt.timezone.utc).timestamp() // 3600)


def europe(out):
    """ICON-EU's newest run: accumulated precipitation and 2 m temperature, hourly to 72 h."""
    now = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, minute=0, second=0, microsecond=0)
    run = None
    for back in range(3, 27, 3):
        t = now - dt.timedelta(hours=back)
        t = t.replace(hour=t.hour // 6 * 6)
        name = f'{DWD}/{t:%H}/tot_prec/icon-eu_europe_regular-lat-lon_single-level_{t:%Y%m%d%H}_072_TOT_PREC.grib2.bz2'
        try:
            urllib.request.urlopen(urllib.request.Request(name, headers=UA, method='HEAD'), timeout=60)
            run = t
            break
        except Exception:
            continue
    if run is None:
        raise SystemExit('no ICON-EU run found')

    def field(var, step):
        url = f'{DWD}/{run:%H}/{var.lower()}/icon-eu_europe_regular-lat-lon_single-level_{run:%Y%m%d%H}_{step:03d}_{var}.grib2.bz2'
        return decode(bz2.decompress(get(url)))

    prev, lat, lon = field('TOT_PREC', 0)
    rates, temps = [], []
    for s in range(1, 73):
        acc, lat, lon = field('TOT_PREC', s)
        rates.append(np.clip(acc - prev, 0, None))
        prev = acc
        t2, _, _ = field('T_2M', s)
        temps.append(t2 - 273.15)
    r, t = np.stack(rates), np.stack(temps)
    n = write_tiles(os.path.join(out, 'pe'), 5, unix_h(run) + 1, 1, lat, lon, [q_rate(r), q_temp(t)])
    print('pe', run, n)
    return {'model': 'ICON-EU', 'run': run.isoformat(), 'tiles': n, 'shape': list(r.shape)}


def world(out):
    """IFS's newest 00 UTC run: precipitation and 2 m temperature every 3 hours to 144 h (by byte ranges)."""
    today = dt.datetime.now(dt.timezone.utc).date()
    for d in (today, today - dt.timedelta(days=1)):
        s = d.strftime('%Y%m%d')
        try:
            get(f'{ECMWF}/{s}/00z/ifs/0p25/oper/{s}000000-144h-oper-fc.index')
            break
        except Exception as e:
            print('no IFS run', s, e)
    else:
        raise SystemExit('no IFS run found')
    run = dt.datetime(d.year, d.month, d.day)

    def fields(step):
        base = f'{ECMWF}/{s}/00z/ifs/0p25/oper/{s}000000-{step}h-oper-fc'
        rows = [json.loads(l) for l in get(base + '.index').decode().splitlines() if l.strip()]
        out = {}
        for r in rows:
            if r.get('param') in ('tp', '2t') and r.get('param') not in out:
                out[r['param']] = decode(get(base + '.grib2', (r['_offset'], r['_offset'] + r['_length'] - 1)))
        return out

    prev = fields(0)['tp'][0]
    rates, temps = [], []
    lat = lon = None
    for step in range(3, 145, 3):
        f = fields(step)
        (acc, lat, lon), (t2, _, _) = f['tp'], f['2t']
        rates.append(np.clip(acc - prev, 0, None) * 1000 / 3)
        prev = acc
        temps.append(t2 - 273.15)
    r, t = np.stack(rates), np.stack(temps)
    n = write_tiles(os.path.join(out, 'pw'), 5, unix_h(run) + 3, 3, lat, lon, [q_rate(r), q_temp(t)])
    print('pw', run, n)
    return {'model': 'ECMWF IFS', 'run': run.isoformat(), 'tiles': n, 'shape': list(r.shape)}


def main():
    args = sys.argv[1:]
    out = args[0]
    only = args[args.index('--only') + 1] if '--only' in args else None
    summary = {}
    for name, fn in (('eu', europe), ('world', world)):
        if only not in (None, name):
            continue
        try:
            summary[name] = fn(out)
        except BaseException as e:
            print(name, 'failed:', e)
    with open(os.path.join(out, 'precip.json'), 'w') as f:
        json.dump({'sets': summary, 'source': 'DWD ICON-EU open data; ECMWF open data (CC BY 4.0)'}, f, indent=1)


if __name__ == '__main__':
    main()
