"""Wind gusts for the whole world, from ECMWF's open high-resolution forecast (CC BY 4.0).

usage: ecmwf_gust.py <out dir>
  Takes the latest 00 UTC IFS run of data.ecmwf.int (0.25°): 10fg, the strongest 10 m gust in the 3 hours before each
  step, for steps 3 … 72 h, fetched by byte ranges from each step's index, and writes <out>/gu/<lat>_<lon>.bin.gz,
  5-degree tiles. MET Norway gives gusts only around Scandinavia; the app uses these elsewhere.

Tile format (little-endian), gzipped:
  'GU' 1, int32 base (Unix hours of the run), uint8 steps, uint8 every (hours; step i ends at base + every*(i+1)),
  uint8 rows, uint8 cols, float32 lat of row 0, float32 lon of col 0, float32 step (rows go north, cols east),
  then uint8 [steps][rows][cols]: gust in 0.5 m/s, 255 = no data.

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
EVERY = 3
STEPS = 24  # 3-72 h


def get(url, rng=None):
    h = dict(UA)
    if rng:
        h['Range'] = 'bytes=%d-%d' % rng
    return urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=300).read()


def run_day():
    today = dt.datetime.now(dt.timezone.utc).date()
    for d in (today, today - dt.timedelta(days=1), today - dt.timedelta(days=2)):
        s = d.strftime('%Y%m%d')
        try:
            get(f'{HOST}/{s}/00z/ifs/0p25/oper/{s}000000-72h-oper-fc.index')
            return d
        except Exception as e:
            print('no run', s, e)
    raise SystemExit('no ECMWF run found')


def main():
    out = os.path.join(sys.argv[1], 'gu')
    import eccodes
    day = run_day()
    s = day.strftime('%Y%m%d')
    fields = []
    lat0 = lon0 = step = None
    for i in range(STEPS):
        h = EVERY * (i + 1)
        base = f'{HOST}/{s}/00z/ifs/0p25/oper/{s}000000-{h}h-oper-fc'
        rows = [json.loads(l) for l in get(base + '.index').decode().splitlines() if l.strip()]
        r = next((x for x in rows if x['param'] == '10fg'), None)
        if r is None:
            fields.append(None)
            continue
        g = eccodes.codes_new_from_message(get(base + '.grib2', (r['_offset'], r['_offset'] + r['_length'] - 1)))
        nj, ni = eccodes.codes_get(g, 'Nj'), eccodes.codes_get(g, 'Ni')
        lat0 = eccodes.codes_get(g, 'latitudeOfFirstGridPointInDegrees')
        lon0 = eccodes.codes_get(g, 'longitudeOfFirstGridPointInDegrees')
        step = eccodes.codes_get(g, 'iDirectionIncrementInDegrees')
        fields.append(eccodes.codes_get_values(g).reshape(nj, ni))
        eccodes.codes_release(g)
    shape = next(f.shape for f in fields if f is not None)
    a = np.stack([f if f is not None else np.full(shape, np.nan) for f in fields])
    lat = lat0 - np.arange(a.shape[1]) * step
    lon = (lon0 + np.arange(a.shape[2]) * step + 180) % 360 - 180
    la_o, lo_o = np.argsort(lat), np.argsort(lon)
    lat, lon, a = lat[la_o], lon[lo_o], a[:, la_o][:, :, lo_o]
    q = np.where(np.isfinite(a), np.clip(np.round(a * 2), 0, 254), 255).astype(np.uint8)
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
            head = b'GU' + struct.pack('<BiBBBBfff', 1, base_h, block.shape[0], EVERY, block.shape[1], block.shape[2],
                                       float(lat[rr[0]]), float(lon[cc[0]]), float(step))
            with gzip.open(os.path.join(out, f'{la}_{lo}.bin.gz'), 'wb', compresslevel=9) as f:
                f.write(head + block.tobytes())
            count += 1
    print('gusts', day, 'steps', len(fields), 'tiles', count)


if __name__ == '__main__':
    main()
