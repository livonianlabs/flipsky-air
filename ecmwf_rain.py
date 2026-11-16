"""The chance of rain by day for the whole world, from ECMWF's open ensemble probabilities (CC BY 4.0).

usage: ecmwf_rain.py <out dir>
  Takes the latest 00 UTC ensemble's "ep" file of data.ecmwf.int (0.25°): tpg1, the chance of at least 1 mm of
  precipitation in 24 hours, for the windows 0-24, 12-36 … 216-240 h, fetched by byte ranges from its index (not the
  whole file), and writes <out>/rp/<lat>_<lon>.bin.gz, 5-degree tiles. MET Norway gives a chance of rain only around
  Scandinavia; the app uses these elsewhere.

Tile format (little-endian), gzipped:
  'RP' 1, int32 base (Unix hours of the run), uint8 windows (window i is base+12i … base+12i+24 h),
  uint8 rows, uint8 cols, float32 lat of row 0, float32 lon of col 0, float32 step (rows go north, cols east),
  then uint8 [windows][rows][cols]: percent, 255 = no data.

Source: ECMWF open data (www.ecmwf.int), CC BY 4.0.
"""
import datetime as dt
import gzip
import json
import os
import struct
import sys
import urllib.request

import numpy as np

HOST = 'https://data.ecmwf.int/forecasts'
UA = {'User-Agent': 'FlipSky data job (github.com/livonianlabs)'}


def get(url, rng=None):
    h = dict(UA)
    if rng:
        h['Range'] = 'bytes=%d-%d' % rng
    return urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=300).read()


def latest():
    """The newest 00 UTC run whose ep file is out: today's, else yesterday's."""
    today = dt.datetime.now(dt.timezone.utc).date()
    for d in (today, today - dt.timedelta(days=1), today - dt.timedelta(days=2)):
        s = d.strftime('%Y%m%d')
        base = f'{HOST}/{s}/00z/ifs/0p25/enfo/{s}000000-240h-enfo-ep'
        try:
            return d, base, get(base + '.index').decode()
        except Exception as e:
            print('no run', s, e)
    raise SystemExit('no ECMWF ep file found')


def main():
    out = os.path.join(sys.argv[1], 'rp')
    import eccodes
    day, base, index = latest()
    rows = [json.loads(l) for l in index.splitlines() if l.strip()]
    wins = sorted((r for r in rows if r['param'] == 'tpg1' and '-' in r['step']), key=lambda r: int(r['step'].split('-')[0]))
    fields = []
    lat0 = lon0 = step = None
    for r in wins:
        start = int(r['step'].split('-')[0])
        if start % 12 or start // 12 != len(fields):
            continue  # keep the regular 12-hour spaced windows only
        msg = get(base + '.grib2', (r['_offset'], r['_offset'] + r['_length'] - 1))
        h = eccodes.codes_new_from_message(msg)
        nj, ni = eccodes.codes_get(h, 'Nj'), eccodes.codes_get(h, 'Ni')
        lat0 = eccodes.codes_get(h, 'latitudeOfFirstGridPointInDegrees')
        lon0 = eccodes.codes_get(h, 'longitudeOfFirstGridPointInDegrees')
        step = eccodes.codes_get(h, 'iDirectionIncrementInDegrees')
        fields.append(eccodes.codes_get_values(h).reshape(nj, ni))
        eccodes.codes_release(h)
    if not fields:
        raise SystemExit('no tpg1 windows')
    a = np.stack(fields)  # windows, lat (north to south), lon (from lon0 east)
    lat = lat0 - np.arange(a.shape[1]) * step
    lon = (lon0 + np.arange(a.shape[2]) * step + 180) % 360 - 180
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    lat, lon, a = lat[la_o], lon[lo_o], a[:, la_o][:, :, lo_o]
    q = np.where(np.isfinite(a), np.clip(np.round(a), 0, 100), 255).astype(np.uint8)
    base_h = int(dt.datetime(day.year, day.month, day.day, tzinfo=dt.timezone.utc).timestamp() // 3600)
    os.makedirs(out, exist_ok=True)
    count = 0
    for la in range(-90, 90, 5):
        rr = np.nonzero((lat >= la) & (lat < la + 5))[0]
        if len(rr) == 0:
            continue
        for lo in range(-180, 180, 5):
            cc = np.nonzero((lon >= lo) & (lon < lo + 5))[0]
            if len(cc) == 0:
                continue
            block = q[:, rr[0]:rr[-1] + 1, cc[0]:cc[-1] + 1]
            head = b'RP' + struct.pack('<BiBBBfff', 1, base_h, block.shape[0], block.shape[1], block.shape[2],
                                       float(lat[rr[0]]), float(lon[cc[0]]), float(step))
            with gzip.open(os.path.join(out, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.tobytes())
            count += 1
    print('rain chance', day, 'windows', len(fields), 'tiles', count)


if __name__ == '__main__':
    main()
