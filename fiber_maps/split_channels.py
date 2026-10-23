#!/usr/bin/env python3
"""Split every 2-channel .tif under a folder into <name>_nucleus.tif and <name>_fibronectin.tif,
written next to the original. Usage: python fiber_maps/split_channels.py <folder>
"""
import sys, os
from PIL import Image

root = sys.argv[1]
done = skip = err = 0
for dirpath, dirs, files in os.walk(root):
    for fn in sorted(files):
        if not fn.lower().endswith('.tif'):
            continue
        if fn.endswith('_nucleus.tif') or fn.endswith('_fibronectin.tif'):
            continue  # don't re-split outputs
        src = os.path.join(dirpath, fn)
        base = fn[:-4]
        out_nuc = os.path.join(dirpath, base + '_nucleus.tif')
        out_fib = os.path.join(dirpath, base + '_fibronectin.tif')
        if os.path.exists(out_nuc) and os.path.exists(out_fib):
            skip += 1
            continue
        try:
            im = Image.open(src)
            nf = getattr(im, 'n_frames', 1)
            if nf != 2:
                print(f"WARN not 2 frames ({nf}): {src}", flush=True)
                err += 1
                continue
            im.seek(0); im.save(out_nuc)
            im.seek(1); im.save(out_fib)
            done += 1
            if done % 200 == 0:
                print(f"  progress: {done} split, {skip} skipped", flush=True)
        except Exception as e:
            print(f"ERROR {src}: {e}", flush=True)
            err += 1
print(f"DONE split={done} skipped={skip} errors={err}", flush=True)
