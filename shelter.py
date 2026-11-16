"""How open each beach is to the sea, direction by direction: a bay is closed on most sides, a cape open on many.

usage: shelter.py <beaches dir> [--land land-polygons-split-4326.zip]
  For every place in <dir>/*.json (beaches.py, beaches_world.py), casts 16 rays (N, NNE, … every 22.5°) over
  OpenStreetMap's land polygons (osmdata.openstreetmap.de, ODbL) and measures the open water in each: from the first
  water within 1 km of the beach to the next land, up to 20 km. Appends to each place a 7th field: 16 integers,
  hundreds of metres of water (0: land that way; 200: open sea, 20 km or more).
  The app compares them with where the wind comes from: wind over a long stretch of water reaches the beach with
  its waves; wind from the land leaves the water calm.

Land polygons © OpenStreetMap contributors (ODbL).
"""
import glob
import json
import math
import os
import sys
import tempfile
import urllib.request
import zipfile

import numpy as np

LAND = 'https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip'
DIRS = 16
STEP_KM = 0.1
REACH_KM = 20.0
START_KM = 1.0


def main():
    args = sys.argv[1:]
    folder = args[0]
    path = args[args.index('--land') + 1] if '--land' in args else None
    if path is None:
        path = os.path.join(tempfile.mkdtemp(), 'land.zip')
        req = urllib.request.Request(LAND, headers={'User-Agent': 'FlipSky data job (github.com/livonianlabs)'})
        with urllib.request.urlopen(req, timeout=1800) as r, open(path, 'wb') as f:
            while True:
                b = r.read(1 << 22)
                if not b:
                    break
                f.write(b)
    import pyogrio
    import shapely
    d = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as z:
        z.extractall(d)
    shp = glob.glob(os.path.join(d, '**', '*.shp'), recursive=True)[0]
    # Raw reading: the geometries as WKB, no geopandas needed.
    geoms = np.asarray(shapely.from_wkb(pyogrio.raw.read(shp, columns=[])[2]))
    tree = shapely.STRtree(geoms)
    print('land polygons', len(geoms), flush=True)
    steps = np.arange(STEP_KM, REACH_KM + STEP_KM / 2, STEP_KM)
    angles = np.radians(np.arange(DIRS) * 360 / DIRS)
    done = 0
    for fp in sorted(glob.glob(os.path.join(folder, '*.json'))):
        rows = json.load(open(fp, encoding='utf-8'))
        if not rows:
            continue
        # The land around all the tile's beaches, once: their rays reach 20 km out.
        lats = [r[1] for r in rows]; lons = [r[2] for r in rows]
        kk = max(math.cos(math.radians(max(abs(min(lats)), abs(max(lats))))), .05)
        pad = REACH_KM / 111.0
        near = tree.query(shapely.box(min(lons) - pad / kk, min(lats) - pad, max(lons) + pad / kk, max(lats) + pad))
        land = shapely.union_all(geoms[near]) if len(near) else None
        if land is not None:
            shapely.prepare(land)
        for r in rows:
            lat, lon = r[1], r[2]
            k = math.cos(math.radians(lat))
            # All the ray points of the beach at once: DIRS × steps.
            dy = np.outer(np.cos(angles), steps) / 111.0
            dx = np.outer(np.sin(angles), steps) / (111.0 * max(k, .05))
            ys, xs = lat + dy, lon + dx
            if land is None:
                fetch = [int(REACH_KM * 10)] * DIRS
            else:
                inside = shapely.contains_xy(land, xs.ravel(), ys.ravel()).reshape(DIRS, len(steps))
                fetch = []
                for i in range(DIRS):
                    row = inside[i]
                    head = ~row[: int(START_KM / STEP_KM)]
                    first = int(np.argmax(head)) if head.any() else -1
                    if first < 0:
                        fetch.append(0)
                        continue
                    rest = row[first:]
                    hit = int(np.argmax(rest)) if rest.any() else len(rest)
                    fetch.append(min(int(round(hit * STEP_KM * 10)), int(REACH_KM * 10)))
            if len(r) < 6:
                r += [''] * (6 - len(r))
            r[6:] = [fetch]
            done += 1
        with open(fp, 'w', encoding='utf-8') as f:
            json.dump(rows, f, ensure_ascii=False, separators=(',', ':'))
        if done % 5000 < len(rows):
            print('beaches', done, flush=True)
    print('shelter done', done)


if __name__ == '__main__':
    main()
