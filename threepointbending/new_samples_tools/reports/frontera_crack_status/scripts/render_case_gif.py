"""
render_case_gif.py -- frames of the damage field inside a translucent 3D specimen for ANY 3PB case (pvpython).

WHAT    the approved C1 format (render_c1_gif_volume.py) for any case: translucent light-grey specimen with dark
        feature edges, the damage d volume-rendered on a 0..1 colour bar, title, key, per-frame label with the
        simulation time t and the load-line displacement (0.1 mm per unit t); frames from --t-from to the last output,
        then a half-circle orbit about the vertical axis at the final state; durations.txt and range.txt for the GIF.
        The resampling slab for the volume is sized from the case's own final damage: nodes with d >= SLAB_D above
        the supports (y >= SLAB_Y_MIN_MM) inside the central window SLAB_X_WINDOW_MM, padded by SLAB_PAD_MM.
STEPS   1 reader + merged grid + shell + edges   2 slab from the final damage   3 styles + labels   4 camera
        5 time frames   6 orbit frames
USAGE   pvpython render_case_gif.py EXODUS OUT_DIR --title "..." --t-from T [--test] [--orbit-only]
"""
import argparse
import math
import os

import numpy as np
from paraview import servermanager
from paraview.simple import (ColorBy, CreateView, ExtractSurface, FeatureEdges, GetColorTransferFunction,
                             GetOpacityTransferFunction, GetScalarBar, MergeBlocks, OpenDataFile, ResampleToImage,
                             SaveScreenshot, Show, Text)
from vtk.util.numpy_support import vtk_to_numpy

IMAGE = [1200, 760]
SHELL_OPACITY = 0.18
SHELL_RGB = [0.72, 0.74, 0.78]
DAMAGE_RGB = [0.0, 0.86, 0.86, 0.86,  0.35, 0.98, 0.72, 0.55,  0.7, 0.89, 0.29, 0.20,  1.0, 0.55, 0.0, 0.05]
DAMAGE_OPACITY = [0.0, 0.0, 0.5, 0.0,  0.08, 0.0, 0.5, 0.0,  0.4, 0.25, 0.5, 0.0,  0.8, 0.75, 0.5, 0.0,  1.0, 1.0, 0.5, 0.0]
OPACITY_UNIT_M = 2.0e-4
VOXEL_M = 5.0e-5                 # 0.05 mm resampling grid
SLAB_D = 0.3                     # damage that sizes the slab
SLAB_Y_MIN_MM = 1.0              # ignore support damage at the bottom
SLAB_X_WINDOW_MM = (7.0, 21.0)   # cracks stay in the middle of the 28 mm beam
SLAB_PAD_MM = 0.75
KEY = "colour: damage d inside the specimen\ngrey: specimen, translucent"
CAMERA_DIR = (0.80, 0.42, 0.85)
CAMERA_DIST_M = 0.047
LOAD_MM_PER_T = 0.1
ORBIT_DEG, ORBIT_FRAMES = 180.0, 60
TIME_MS, HOLD_FINAL_TIME_MS, ORBIT_MS, HOLD_END_MS = 100, 1200, 80, 1500


def rotate_about_y(direction, degrees):
    a = math.radians(degrees)
    x, y, z = direction
    return (x * math.cos(a) + z * math.sin(a), y, -x * math.sin(a) + z * math.cos(a))


def build(path):
    reader = OpenDataFile(path)
    for prop in ("NodeBlockFields", "PointVariables"):
        if reader.GetProperty(prop) is not None:
            setattr(reader, prop, ["d"])
            break
    merged = MergeBlocks(Input=reader)
    shell = ExtractSurface(Input=merged)
    edges = FeatureEdges(Input=shell)
    edges.FeatureAngle = 75.0
    edges.BoundaryEdges = 0
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0
    return reader, merged, shell, edges


def slab_bounds(merged, t_last, box):
    """Resampling box: the final damage's x-extent (padded), full height and depth of the specimen (metres)."""
    merged.UpdatePipeline(t_last)
    grid = servermanager.Fetch(merged)
    xyz = vtk_to_numpy(grid.GetPoints().GetData())
    d = vtk_to_numpy(grid.GetPointData().GetArray("d"))
    s = 1e-3 if box[1] < 1.0 else 1.0
    x, y = xyz[:, 0] / s, xyz[:, 1] / s
    sel = (d >= SLAB_D) & (y >= SLAB_Y_MIN_MM) & (x >= SLAB_X_WINDOW_MM[0]) & (x <= SLAB_X_WINDOW_MM[1])
    x0, x1 = (x[sel].min(), x[sel].max()) if sel.any() else (13.0, 15.0)
    x0, x1 = max(box[0] / s, x0 - SLAB_PAD_MM), min(box[1] / s, x1 + SLAB_PAD_MM)
    bounds = [float(v) for v in (x0 * s, x1 * s, box[2], box[3], box[4], box[5])]   # plain floats for ParaView
    dims = [int(round((bounds[2 * k + 1] - bounds[2 * k]) / VOXEL_M)) + 1 for k in range(3)]
    return bounds, dims


def solid(view, source, rgb, opacity=1.0):
    shown = Show(source, view)
    shown.ColorArrayName = ["POINTS", ""]
    shown.DiffuseColor = rgb
    shown.AmbientColor = rgb
    shown.Opacity = opacity
    shown.Specular = 0.0
    return shown


def style(view, damage, shell, edges):
    solid(view, shell, SHELL_RGB, SHELL_OPACITY)
    lines = solid(view, edges, [0.2, 0.2, 0.2])
    lines.LineWidth = 1.5
    volume = Show(damage, view)
    volume.Representation = "Volume"
    ColorBy(volume, ("POINTS", "d"))
    lut = GetColorTransferFunction("d")
    lut.RGBPoints = DAMAGE_RGB
    lut.ColorSpace = "Lab"
    lut.RescaleTransferFunction(0.0, 1.0)
    lut.AutomaticRescaleRangeMode = "Never"
    pwf = GetOpacityTransferFunction("d")
    pwf.Points = DAMAGE_OPACITY
    pwf.RescaleTransferFunction(0.0, 1.0)
    volume.ScalarOpacityUnitDistance = OPACITY_UNIT_M
    volume.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    for name, value in {"Title": "damage d", "ComponentTitle": "", "TitleColor": [0.1, 0.1, 0.1],
                        "LabelColor": [0.1, 0.1, 0.1], "TitleFontSize": 18, "LabelFontSize": 16,
                        "AutomaticLabelFormat": 0, "LabelFormat": "{:.2f}", "RangeLabelFormat": "{:.0f}",
                        "UseCustomLabels": 1, "CustomLabels": [0.25, 0.5, 0.75], "AddRangeLabels": 1,
                        "WindowLocation": "Any Location", "Position": [0.90, 0.16], "ScalarBarLength": 0.5}.items():
        if bar.GetProperty(name) is not None:
            setattr(bar, name, value)


def labels(view, title):
    for text, size, where in ((title, 17, "Upper Left Corner"), (KEY, 15, "Lower Right Corner")):
        shown = Show(Text(Text=text), view)
        shown.FontSize = size
        shown.Color = [0.1, 0.1, 0.1]
        shown.WindowLocation = where
    clock_src = Text(Text="")
    clock = Show(clock_src, view)
    clock.FontSize = 22
    clock.Color = [0.1, 0.1, 0.1]
    clock.WindowLocation = "Lower Left Corner"
    return clock_src


def camera(view, bounds, direction):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    centre = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)]
    n = sum(c * c for c in direction) ** 0.5
    view.CameraFocalPoint = centre
    view.CameraPosition = [centre[k] + CAMERA_DIST_M * direction[k] / n for k in range(3)]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraViewAngle = 30.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exodus")
    ap.add_argument("out_dir")
    ap.add_argument("--title", required=True)
    ap.add_argument("--t-from", type=float, required=True)
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--orbit-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    reader, merged, shell, edges = build(a.exodus)
    times = [t for t in reader.TimestepValues if t >= a.t_from - 1e-9]
    reader.UpdatePipeline(times[-1])                   # bounds are a +-1e299 sentinel until data is loaded
    box = reader.GetDataInformation().GetBounds()
    bounds, dims = slab_bounds(merged, times[-1], box)
    damage = ResampleToImage(Input=merged)
    damage.UseInputBounds = 0
    damage.SamplingBounds = bounds
    damage.SamplingDimensions = dims
    print("slab x %.2f..%.2f mm, grid %s" % (bounds[0] * 1e3, bounds[1] * 1e3, dims))
    view = CreateView("RenderView")
    view.ViewSize = IMAGE
    view.OrientationAxesVisibility = 1
    for name, value in (("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"),
                        ("Background", [1.0, 1.0, 1.0]), ("UseDepthPeeling", 1), ("MaximumNumberOfPeels", 8),
                        ("UseDepthPeelingForVolumes", 1)):
        if view.GetProperty(name) is not None:
            setattr(view, name, value)
    style(view, damage, shell, edges)
    clock_src = labels(view, a.title)
    final = "t = %.1f (final)    load-line displacement %.2f mm" % (times[-1], LOAD_MM_PER_T * times[-1])
    if a.test:
        view.ViewTime = times[-1]
        clock_src.Text = final
        camera(view, box, CAMERA_DIR)
        SaveScreenshot(os.path.join(a.out_dir, "test_main.png"), view, ImageResolution=IMAGE)
        print("wrote test view")
        return
    camera(view, box, CAMERA_DIR)
    if not a.orbit_only:
        for k, t in enumerate(times):
            view.ViewTime = t
            clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (t, LOAD_MM_PER_T * t)
            SaveScreenshot(os.path.join(a.out_dir, "frame_%03d.png" % k), view, ImageResolution=IMAGE)
    view.ViewTime = times[-1]
    clock_src.Text = final
    for j in range(1, ORBIT_FRAMES + 1):
        camera(view, box, rotate_about_y(CAMERA_DIR, ORBIT_DEG * j / ORBIT_FRAMES))
        SaveScreenshot(os.path.join(a.out_dir, "frame_%03d.png" % (len(times) - 1 + j)), view, ImageResolution=IMAGE)
    durations = [TIME_MS] * (len(times) - 1) + [HOLD_FINAL_TIME_MS] + [ORBIT_MS] * (ORBIT_FRAMES - 1) + [HOLD_END_MS]
    with open(os.path.join(a.out_dir, "durations.txt"), "w") as fh:
        fh.write("\n".join(str(ms) for ms in durations) + "\n")
    with open(os.path.join(a.out_dir, "range.txt"), "w") as fh:
        fh.write("t%g-%g\n" % (round(times[0], 1), round(times[-1], 1)))     # e.g. t18-26.6, as for C1
    print("wrote %d time frames (t = %.1f .. %.1f) + %d orbit frames" % (len(times), times[0], times[-1], ORBIT_FRAMES))


if __name__ == "__main__":
    main()
