"""
render_c1_gif_broken.py -- frames of the BROKEN specimen of case C1: the fully damaged band removed, for a GIF.

WHAT    the phase-field crack is a band of material with d ~ 1, so its d = 0.95 contour is two parallel sheets with
        nothing drawn between them. Here the band's cells (any node with d > D_BROKEN) are removed instead, so the
        specimen's own boundary runs along the crack faces and the specimen reads as split, joined only by the
        uncracked ligament. Optional warp by the computed displacement (WARP_SCALE, 1 = true scale) opens the gap.
        The specimen surface is coloured by d (grey 0 -> red 1, crack faces red), translucent (SHELL_OPACITY),
        with dark feature edges. Title, colour bar, per-frame t and load-line displacement; after the last time a
        half-circle orbit about the vertical axis (durations.txt for make_gif.py).
STEPS   1 reader (d, disp) -> merged grid -> threshold d <= D_BROKEN -> optional warp -> surface + edges
        2 style + labels  3 camera  4 time frames  5 orbit frames
USAGE   pvpython render_c1_gif_broken.py EXODUS OUT_DIR [--test] [--warp SCALE] [--opacity A]
        (--test: final state only, from the main camera and from the right end, to choose the settings)
"""
import math
import os
import sys

from paraview.simple import (Contour, CreateView, ExtractSurface, FeatureEdges, MergeBlocks, OpenDataFile,
                             SaveScreenshot, Show, Text, Threshold, WarpByVector)

T_FROM = 18.0
IMAGE = [1200, 760]
D_BROKEN = 0.95                  # cells touching a node above this are removed (the open crack)
D_FACES = 0.98                   # crack faces drawn just inside the removed band, clear of the jagged cut faces
TITLE = "C1   (c) two holes, not connected: xL 3.00 / xR 2.00 mm, r 1.0 mm   |   broken material (d > 0.95) removed"
CAMERA_DIR = (0.80, 0.42, 0.85)
CAMERA_DIST_M = 0.047
LOAD_MM_PER_T = 0.1
ORBIT_DEG, ORBIT_FRAMES = 180.0, 60
TIME_MS, HOLD_FINAL_TIME_MS, ORBIT_MS, HOLD_END_MS = 100, 1200, 80, 1500


def option(name, default):
    return float(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def build(path, warp_scale):
    reader = OpenDataFile(path)
    for prop in ("NodeBlockFields", "PointVariables"):
        if reader.GetProperty(prop) is not None:
            setattr(reader, prop, list(getattr(reader, prop).Available))   # d and the displacement
            break
    merged = MergeBlocks(Input=reader)
    intact = Threshold(Input=merged, Scalars=["POINTS", "d"])
    for name, value in (("LowerThreshold", -1.0), ("UpperThreshold", D_BROKEN), ("ThresholdMethod", "Between"),
                        ("AllScalars", 1)):
        if intact.GetProperty(name) is not None:
            setattr(intact, name, value)
    body = intact
    faces = Contour(Input=merged, ContourBy=["POINTS", "d"], Isosurfaces=[D_FACES])    # smooth crack faces
    if warp_scale > 0:
        body = WarpByVector(Input=intact, Vectors=["POINTS", "disp"], ScaleFactor=warp_scale)
        faces = WarpByVector(Input=faces, Vectors=["POINTS", "disp"], ScaleFactor=warp_scale)
    shell = ExtractSurface(Input=body)
    # outline from the UNCUT specimen surface (box, holes, notch), moved with the same displacement: edges taken
    # from the cut shell would trace every jagged facet of the crack faces as dark specks
    edges = FeatureEdges(Input=ExtractSurface(Input=merged))
    edges.FeatureAngle = 75.0
    edges.BoundaryEdges = 0
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0
    if warp_scale > 0:
        edges = WarpByVector(Input=edges, Vectors=["POINTS", "disp"], ScaleFactor=warp_scale)
    return reader, shell, edges, faces


def solid(view, source, rgb, opacity=1.0, ambient=0.0, diffuse=1.0):
    shown = Show(source, view)
    shown.ColorArrayName = ["POINTS", ""]
    shown.DiffuseColor = rgb
    shown.AmbientColor = rgb
    shown.Opacity = opacity
    shown.Specular = 0.0
    shown.Ambient = ambient
    shown.Diffuse = diffuse
    return shown


def style(view, shell, edges, faces, opacity):
    solid(view, shell, [0.72, 0.74, 0.78], opacity, ambient=0.6, diffuse=0.4)   # flat light: jagged cut facets read even
    solid(view, faces, [0.78, 0.06, 0.08])                   # the crack faces, opaque red
    lines = Show(edges, view)
    lines.ColorArrayName = ["POINTS", ""]
    lines.DiffuseColor = [0.2, 0.2, 0.2]
    lines.AmbientColor = [0.2, 0.2, 0.2]
    lines.LineWidth = 1.5


def labels(view, warp_scale):
    title = Show(Text(Text=TITLE), view)
    title.FontSize = 17
    title.Color = [0.1, 0.1, 0.1]
    title.WindowLocation = "Upper Left Corner"
    note = ("red: crack faces (broken band, d > 0.95, removed)      grey: specimen, translucent\n"
            + ("shape: computed displacement x %g" % warp_scale if warp_scale > 0 else "shape: undeformed"))
    key = Show(Text(Text=note), view)
    key.FontSize = 15
    key.Color = [0.1, 0.1, 0.1]
    key.WindowLocation = "Lower Right Corner"
    clock_src = Text(Text="")
    clock = Show(clock_src, view)
    clock.FontSize = 22
    clock.Color = [0.1, 0.1, 0.1]
    clock.WindowLocation = "Lower Left Corner"
    return clock_src


def rotate_about_y(direction, degrees):
    a = math.radians(degrees)
    x, y, z = direction
    return (x * math.cos(a) + z * math.sin(a), y, -x * math.sin(a) + z * math.cos(a))


def camera(view, bounds, direction):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    centre = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)]
    n = sum(c * c for c in direction) ** 0.5
    view.CameraFocalPoint = centre
    view.CameraPosition = [centre[k] + CAMERA_DIST_M * direction[k] / n for k in range(3)]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraViewAngle = 30.0


def main(path, out_dir):
    test = "--test" in sys.argv
    warp_scale = option("--warp", 0.0)
    opacity = option("--opacity", 0.35)
    os.makedirs(out_dir, exist_ok=True)
    reader, shell, edges, faces = build(path, warp_scale)
    times = [t for t in reader.TimestepValues if t >= T_FROM - 1e-9]
    view = CreateView("RenderView")
    view.ViewSize = IMAGE
    view.OrientationAxesVisibility = 1
    for name, value in (("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"),
                        ("Background", [1.0, 1.0, 1.0]), ("UseDepthPeeling", 1), ("MaximumNumberOfPeels", 8)):
        if view.GetProperty(name) is not None:
            setattr(view, name, value)
    style(view, shell, edges, faces, opacity)
    clock_src = labels(view, warp_scale)
    reader.UpdatePipeline(times[0])
    bounds = reader.GetDataInformation().GetBounds()      # undeformed box: the camera does not follow the warp
    final = "t = %.1f (final)    load-line displacement %.2f mm" % (times[-1], LOAD_MM_PER_T * times[-1])
    if test:
        view.ViewTime = times[-1]
        clock_src.Text = final
        tag = "w%g_o%g" % (warp_scale, opacity)
        for name, direction in (("main", CAMERA_DIR), ("rightend", rotate_about_y(CAMERA_DIR, 45.0))):
            camera(view, bounds, direction)
            SaveScreenshot(os.path.join(out_dir, "test_%s_%s.png" % (tag, name)), view, ImageResolution=IMAGE)
        print("wrote test views", tag)
        return
    camera(view, bounds, CAMERA_DIR)
    for k, t in enumerate(times):
        view.ViewTime = t
        clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (t, LOAD_MM_PER_T * t)
        SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % k), view, ImageResolution=IMAGE)
    clock_src.Text = final
    for j in range(1, ORBIT_FRAMES + 1):
        camera(view, bounds, rotate_about_y(CAMERA_DIR, ORBIT_DEG * j / ORBIT_FRAMES))
        SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % (len(times) - 1 + j)), view, ImageResolution=IMAGE)
    durations = [TIME_MS] * (len(times) - 1) + [HOLD_FINAL_TIME_MS] + [ORBIT_MS] * (ORBIT_FRAMES - 1) + [HOLD_END_MS]
    with open(os.path.join(out_dir, "durations.txt"), "w") as fh:
        fh.write("\n".join(str(ms) for ms in durations) + "\n")
    print("wrote %d time frames + %d orbit frames" % (len(times), ORBIT_FRAMES))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
