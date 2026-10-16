"""
render_c1_gif.py -- frames of the damage field on the 3D specimen of case C1, for a GIF (pvpython).

WHAT    outside surface of the beam coloured by damage d (light grey 0 -> dark red 1), dark feature edges for the
        shape, an oblique camera that shows the front face (crack), the top and one end, a title, and a per-frame
        label with the simulation time t and the load-line displacement (0.1 mm per unit t). One PNG per output
        time from T_FROM to the end.
STEPS   1 reader (d only) + surface + feature edges   2 colour map + scalar bar + labels   3 camera   4 frames
USAGE   pvpython render_c1_gif.py EXODUS OUT_DIR [--test]     (--test: only the last time, to check the camera)
"""
import os
import sys

from paraview.simple import (CreateView, ExtractSurface, FeatureEdges, GetColorTransferFunction, GetScalarBar,
                             MergeBlocks, OpenDataFile, SaveScreenshot, Show, Text, ColorBy)

T_FROM = 18.0                    # first frame: damage builds at the notch before the crack jumps at t ~ 21.6
IMAGE = [1200, 760]
TITLE = "C1   (c) two holes, not connected: xL 3.00 / xR 2.00 mm, r 1.0 mm   |   damage d on the specimen surface"
CAMERA_DIR = (0.42, 0.34, 1.0)   # from the specimen centre toward the camera: front face, top and right end
CAMERA_DIST_M = 0.049
LOAD_MM_PER_T = 0.1              # func_loading = -1e-4 * t (m)


def build(path):
    reader = OpenDataFile(path)
    for prop in ("NodeBlockFields", "PointVariables"):
        if reader.GetProperty(prop) is not None:
            setattr(reader, prop, ["d"])
            break
    merged = MergeBlocks(Input=reader)          # one grid: block interfaces (support/load patches) are not surfaces
    surface = ExtractSurface(Input=merged)
    edges = FeatureEdges(Input=surface)
    edges.FeatureAngle = 75.0        # true edges are ~90 deg; lower angles catch mesh folds at the BC lines
    edges.BoundaryEdges = 0          # open patch boundaries (support and load lines) only add squiggles
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0
    return reader, surface, edges


def style(view, surface, edges):
    shown = Show(surface, view)
    shown.Ambient = 0.25             # lift the faces turned away from the light
    shown.Diffuse = 0.85
    ColorBy(shown, ("POINTS", "d"))
    lut = GetColorTransferFunction("d")
    lut.RGBPoints = [0.0, 0.86, 0.86, 0.86,  0.35, 0.98, 0.72, 0.55,  0.7, 0.89, 0.29, 0.20,  1.0, 0.55, 0.0, 0.05]
    lut.ColorSpace = "Lab"
    lut.RescaleTransferFunction(0.0, 1.0)
    lut.AutomaticRescaleRangeMode = "Never"
    shown.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    for name, value in {"Title": "damage d", "ComponentTitle": "", "TitleColor": [0.1, 0.1, 0.1],
                        "LabelColor": [0.1, 0.1, 0.1], "TitleFontSize": 18, "LabelFontSize": 16,
                        "AutomaticLabelFormat": 0, "LabelFormat": "{:.2f}", "RangeLabelFormat": "{:.0f}",
                        "UseCustomLabels": 1, "CustomLabels": [0.25, 0.5, 0.75], "AddRangeLabels": 1,
                        "WindowLocation": "Any Location", "Position": [0.90, 0.12], "ScalarBarLength": 0.5}.items():
        if bar.GetProperty(name) is not None:
            setattr(bar, name, value)
    lines = Show(edges, view)
    lines.DiffuseColor = [0.15, 0.15, 0.15]
    lines.AmbientColor = [0.15, 0.15, 0.15]
    lines.LineWidth = 1.5
    ColorBy(lines, None)


def labels(view):
    title = Show(Text(Text=TITLE), view)
    title.FontSize = 17
    title.Color = [0.1, 0.1, 0.1]
    title.WindowLocation = "Upper Left Corner"
    clock_src = Text(Text="")
    clock = Show(clock_src, view)
    clock.FontSize = 22
    clock.Color = [0.1, 0.1, 0.1]
    clock.WindowLocation = "Lower Left Corner"
    return clock_src


def camera(view, bounds):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    centre = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)]
    n = sum(c * c for c in CAMERA_DIR) ** 0.5
    view.CameraFocalPoint = centre
    view.CameraPosition = [centre[k] + CAMERA_DIST_M * CAMERA_DIR[k] / n for k in range(3)]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraViewAngle = 30.0


def main(path, out_dir, test):
    os.makedirs(out_dir, exist_ok=True)
    reader, surface, edges = build(path)
    times = [t for t in reader.TimestepValues if t >= T_FROM - 1e-9]
    if test:
        times = times[-1:]
    view = CreateView("RenderView")
    view.ViewSize = IMAGE
    view.OrientationAxesVisibility = 1
    for name, value in (("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"),
                        ("Background", [1.0, 1.0, 1.0])):
        if view.GetProperty(name) is not None:
            setattr(view, name, value)
    style(view, surface, edges)
    clock_src = labels(view)
    reader.UpdatePipeline(times[0])
    camera(view, reader.GetDataInformation().GetBounds())
    for k, t in enumerate(times):
        view.ViewTime = t
        clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (t, LOAD_MM_PER_T * t)
        SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % k), view, ImageResolution=IMAGE)
    print("wrote %d frames, t = %.1f .. %.1f" % (len(times), times[0], times[-1]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--test" in sys.argv)
