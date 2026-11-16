"""The world's named beaches from OpenStreetMap (via Overture Maps' base theme, land class "beach"), merged into the
bathing places' tiles that beaches.py writes from the EU register.

usage: beaches_world.py <beaches dir> [--csv world_beaches.csv] [--release 2026-10-21.0]
  Without --csv, queries Overture's GeoParquet on S3 with DuckDB (only the columns needed, about a minute).
  Each beach becomes [name, lat, lon, "o", 0, english name or ""] in <dir>/<lat>_<lon>.json; one within 1 km of a
  register place, or of a beach of the same name already kept within 2 km, is left out.
  Run after beaches.py, once a year.

Beaches © OpenStreetMap contributors (ODbL), via Overture Maps.
"""
import csv
import json
import math
import os
import re
import sys


def query(release, out):
    import duckdb
    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs; SET s3_region='us-west-2'; SET enable_progress_bar=false;")
    u = f's3://overturemaps-us-west-2/release/{release}/theme=base/type=land/*.parquet'
    con.execute(f"""COPY (SELECT names.primary AS name, names.common['en'] AS en,
        round((bbox.ymin+bbox.ymax)/2, 5) AS lat, round((bbox.xmin+bbox.xmax)/2, 5) AS lon
        FROM read_parquet('{u}') WHERE class='beach' AND names.primary IS NOT NULL) TO '{out}' (HEADER, DELIMITER ',')""")


def km(a, b, c, d):
    return math.hypot((a - c) * 111, (b - d) * 111 * math.cos(math.radians(a)))


def main():
    args = sys.argv[1:]
    folder = args[0]
    path = args[args.index('--csv') + 1] if '--csv' in args else os.path.join(folder, '..', 'world_beaches.csv')
    if '--csv' not in args:
        query(args[args.index('--release') + 1] if '--release' in args else '2026-10-21.0', path)
    tiles = {}
    for f in os.listdir(folder):
        if f.endswith('.json'):
            # Only the register's own (a rerun starts again from them).
            tiles[f[:-5]] = [r for r in json.load(open(os.path.join(folder, f), encoding='utf-8')) if r[3] != 'o']
    added = 0
    for row in csv.DictReader(open(path, encoding='utf-8')):
        try:
            lat, lon = float(row['lat']), float(row['lon'])
        except (TypeError, ValueError):
            continue
        name = (row['name'] or '').strip()
        en = (row.get('en') or '').strip()
        if not name or len(name) > 60:
            continue
        key = f'{math.floor(lat)}_{math.floor(lon)}'
        here = tiles.setdefault(key, [])
        near = False
        for r in here:
            d = km(lat, lon, r[1], r[2])
            # A register name the source mangled ("Playa Fa?ab?"): OpenStreetMap's name for the same place (within 1 km).
            if r[3] != 'o' and '?' in r[0] and d < 1:
                r[0] = name
            if (r[3] != 'o' and d < 1) or (r[3] == 'o' and r[0] == name and d < 2):
                near = True
                break
        if near:
            continue
        here.append([name, round(lat, 5), round(lon, 5), 'o', 0, en if en != name else ''])
        added += 1
    # Names still mangled: a part in brackets dropped ("Playa Troya Ii (am?ricas Ii)"), else the place left out.
    for k, v in tiles.items():
        kept = []
        for r in v:
            if '?' in r[0]:
                r[0] = re.sub(r'\s*\([^)]*\?[^)]*\)', '', r[0]).strip()
                if '?' in r[0] or not r[0]:
                    continue
            kept.append(r)
        tiles[k] = kept
    for k, v in tiles.items():
        if v:
            with open(os.path.join(folder, k + '.json'), 'w', encoding='utf-8') as f:
                json.dump(v, f, ensure_ascii=False, separators=(',', ':'))
    print('beaches added', added, 'tiles', sum(1 for v in tiles.values() if v))


if __name__ == '__main__':
    main()
