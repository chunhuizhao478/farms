#!/usr/bin/env python3
"""
make_gif.py -- assemble numbered PNG frames into a looping GIF with one shared palette.

WHAT    reads frame_000.png, frame_001.png, ... from a folder, builds ONE 256-colour palette from a montage of
        evenly spaced frames (so colours do not flicker between frames), quantizes every frame with it (no
        dithering), and writes a looping GIF: FRAME_MS per frame, the last frame held HOLD_LAST_MS.
USAGE   python make_gif.py FRAME_DIR OUT.gif
"""
import glob
import os
import sys

from PIL import Image

FRAME_MS = 100          # 10 frames per second
HOLD_LAST_MS = 2500     # pause on the final state before looping
PALETTE_SAMPLES = 8     # frames combined to build the shared palette


def main(frame_dir, out):
    paths = sorted(glob.glob(os.path.join(frame_dir, "frame_*.png")))
    frames = [Image.open(p).convert("RGB") for p in paths]
    w, h = frames[0].size
    picks = frames[:: max(1, len(frames) // PALETTE_SAMPLES)][:PALETTE_SAMPLES] + [frames[-1]]
    montage = Image.new("RGB", (w, h * len(picks)))
    for k, im in enumerate(picks):
        montage.paste(im, (0, k * h))
    palette = montage.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    quantized = [im.quantize(palette=palette, dither=Image.Dither.NONE) for im in frames]
    timing = os.path.join(frame_dir, "durations.txt")          # per-frame ms written by the renderer, if any
    if os.path.exists(timing):
        durations = [int(line) for line in open(timing) if line.strip()]
        assert len(durations) == len(quantized), "durations.txt has %d lines for %d frames" % (len(durations), len(quantized))
    else:
        durations = [FRAME_MS] * (len(quantized) - 1) + [HOLD_LAST_MS]
    quantized[0].save(out, save_all=True, append_images=quantized[1:], duration=durations, loop=0, optimize=False)
    print("wrote %s: %d frames, %dx%d, %.1f MB" % (out, len(quantized), w, h, os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
