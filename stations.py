"""Weather stations' normals for FlipSky's trips: what the days and nights and the rain are usually like, half-month by
half-month, at the real stations of NOAA's GHCN-Daily (public domain), so islands and coasts get their own climate
(Tenerife's south is not its north).

usage: stations.py <out dir> [--from 2015] [--to 2024]
  Reads NCEI's by-year files (all stations, TMAX/TMIN/PRCP), keeps the values without a quality flag, and for each
  station and half-month (1-15, 16-end) averages the day's high and low and counts the days with 1 mm of rain or more.
  A half-month counts with at least 5 years of it. Writes <out>/<lat>_<lon>.json, 5-degree tiles:
  [[name, lat, lon, elevation m, [24 highs ×10], [24 lows ×10], [24 rainy-day percents]], ...] (null where too little).
Run by hand once a year (the stations workflow); committed to flipsky-air's stations/.
"""
import gzip
import io
import json
import math
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

BASE = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily'
UA = {'User-Agent': 'FlipSky data job (github.com/livonianlabs)'}


def fetch(path):
    return urllib.request.urlopen(urllib.request.Request(f'{BASE}/{path}', headers=UA), timeout=600)


def main():
    args = sys.argv[1:]
    out = args[0]
    y0 = int(args[args.index('--from') + 1]) if '--from' in args else 2015
    y1 = int(args[args.index('--to') + 1]) if '--to' in args else 2024
    # Per chunk: sums and counts per (station, element, half-month, year); rain as wet days of all days.
    parts = []
    for year in range(y0, y1 + 1):
        raw = fetch(f'by_year/{year}.csv.gz').read()
        print('year', year, len(raw) // 1_000_000, 'MB', flush=True)
        for chunk in pd.read_csv(io.BytesIO(raw), compression='gzip', header=None, chunksize=5_000_000,
                                 names=['id', 'date', 'el', 'v', 'm', 'q', 's', 't'], usecols=['id', 'date', 'el', 'v', 'q'],
                                 dtype={'id': str, 'date': np.int32, 'el': str, 'v': np.int32, 'q': str}):
            c = chunk[chunk['el'].isin(['TMAX', 'TMIN', 'PRCP']) & chunk['q'].isna()]
            month = (c['date'] // 100) % 100
            day = c['date'] % 100
            v = np.where(c['el'] == 'PRCP', (c['v'] >= 10).astype(np.int32), c['v'])
            c = pd.DataFrame({'id': c['id'], 'el': c['el'], 'b': ((month - 1) * 2 + (day > 15)).astype(np.int8), 'y': np.int16(year), 'v': v})
            parts.append(c.groupby(['id', 'el', 'b', 'y'])['v'].agg(['sum', 'count']).reset_index())
        del raw
    allp = pd.concat(parts).groupby(['id', 'el', 'b', 'y'])[['sum', 'count']].sum().reset_index()
    agg = allp.groupby(['id', 'el', 'b']).agg(total=('sum', 'sum'), n=('count', 'sum'), years=('y', 'nunique')).reset_index()
    agg = agg[(agg['years'] >= 5) & (agg['n'] > 0)]
    per = {}
    for sid, el, b, total, n in zip(agg['id'], agg['el'], agg['b'], agg['total'], agg['n']):
        d = per.setdefault(sid, {'TMAX': [None] * 24, 'TMIN': [None] * 24, 'PRCP': [None] * 24})
        # Temperatures in tenths of °C already; rain as the percent of days with 1 mm or more.
        d[el][int(b)] = int(round(total / n)) if el != 'PRCP' else int(round(100 * total / n))
    stations = {}
    for line in fetch('ghcnd-stations.txt').read().decode('latin-1').splitlines():
        sid = line[0:11]
        stations[sid] = (float(line[12:20]), float(line[21:30]), float(line[31:37]), line[41:71].strip())
    tiles = {}
    kept = 0
    for sid, d in per.items():
        if sid not in stations or sum(v is not None for v in d['TMAX']) < 12:
            continue
        lat, lon, elev, name = stations[sid]
        key = f'{math.floor(lat / 5) * 5}_{math.floor(lon / 5) * 5}'
        tiles.setdefault(key, []).append([name.title(), round(lat, 4), round(lon, 4), round(elev), d['TMAX'], d['TMIN'], d['PRCP']])
        kept += 1
    os.makedirs(out, exist_ok=True)
    for k, v in tiles.items():
        with open(os.path.join(out, k + '.json'), 'w', encoding='utf-8') as f:
            json.dump(v, f, ensure_ascii=False, separators=(',', ':'))
    print('stations', kept, 'tiles', len(tiles))


if __name__ == '__main__':
    main()
