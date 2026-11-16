"""Bathing places for FlipSky's water card: the EU bathing water register (EEA, via EMODnet Human Activities, CC BY 4.0)
cut into 1-degree tiles the app downloads.

usage: beaches.py <out dir>
  Writes <out>/<lat>_<lon>.json (lower-left corner, whole degrees): [[name, lat, lon, kind, quality], ...]
  kind: c coastal, t transitional (river mouths, gulfs), l lake, r river; quality: 1 excellent, 2 good, 3 sufficient,
  4 poor, 0 not classified. Each place once, from its latest reported season.
Run by hand once a year (the register changes with each bathing season) and commit the tiles to flipsky-air's beaches/.
"""
import json
import os
import re
import sys
import urllib.parse
import urllib.request

WFS = 'https://ows.emodnet-humanactivities.eu/wfs'
KIND = {'Coastal Bathing Water': 'c', 'Transitional Bathing Water': 't', 'Lake Bathing Water': 'l', 'River Bathing Water': 'r'}


def fetch(year):
    """All places reported for [year], in pages."""
    out, start = [], 0
    while True:
        q = urllib.parse.urlencode({
            'service': 'WFS', 'version': '2.0.0', 'request': 'GetFeature', 'typeNames': 'emodnet:bathingwaters',
            'outputFormat': 'application/json', 'count': 10000, 'startIndex': start, 'CQL_FILTER': f"year='{year}'",
            'propertyName': 'bwid,lat,lon,bwname,bwatercat,class,year',
        })
        req = urllib.request.Request(WFS + '?' + q, headers={'User-Agent': 'FlipSky data job (github.com/livonianlabs)'})
        page = json.load(urllib.request.urlopen(req, timeout=300))['features']
        out += [f['properties'] for f in page]
        if len(page) < 10000:
            return out
        start += 10000


def nice(name):
    """'BALTIJAS JŪRA, PELDVIETA "VĀRZAS"' -> 'Vārzas'; 'VECĀĶI' -> 'Vecāķi'."""
    quoted = re.findall(r'"([^"]+)"', name)
    n = quoted[-1] if quoted else name
    # "Baltijas jūra, Ainaži": the sea's name before the place's.
    n = n.split(',')[-1].strip() or n.strip()
    if n.isupper():
        n = ' '.join(w[:1] + w[1:].lower() for w in n.split(' '))
        n = re.sub(r'-(\w)', lambda m: '-' + m.group(1).upper(), n)
    return n


def main():
    out = sys.argv[1]
    latest = {}
    for year in ('2021', '2022', '2023', '2024', '2025'):
        for p in fetch(year):
            if p.get('lat') is None or p.get('lon') is None:
                continue
            latest[p['bwid']] = p
        print(year, len(latest))
    tiles = {}
    for p in latest.values():
        kind = KIND.get(p.get('bwatercat'), 'c')
        q = (p.get('class') or '0')[:1]
        q = int(q) if q.isdigit() else 0
        la, lo = int(p['lat'] // 1), int(p['lon'] // 1)
        tiles.setdefault(f'{la}_{lo}', []).append([nice(p['bwname'] or p['bwid']), round(p['lat'], 5), round(p['lon'], 5), kind, q])
    os.makedirs(out, exist_ok=True)
    for k, v in tiles.items():
        with open(os.path.join(out, k + '.json'), 'w', encoding='utf-8') as f:
            json.dump(v, f, ensure_ascii=False, separators=(',', ':'))
    print('places', len(latest), 'tiles', len(tiles))


if __name__ == '__main__':
    main()
