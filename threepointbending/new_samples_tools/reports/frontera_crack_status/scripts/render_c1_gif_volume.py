"""
render_c1_gif_volume.py -- frames of the damage field inside a translucent 3D specimen, case C1, for a GIF (pvpython).

WHAT    the X-ray view of render_c1_gif_xray.py with the damage itself filling the crack: the beam's outer surface
        translucent light grey with dark feature edges, and the damage d VOLUME-rendered inside it on a 0..1 colour
        bar (light grey -> orange -> dark red; opacity rises with d, so intact material stays clear and the broken
        band reads solid). A title, a key, the colour bar and a per-frame label with the simulation time t and the
        load-line displacement (0.1 mm per unit t). Frames from T_FROM to the last output, then a half-circle orbit
        about the vertical axis at the final state; durations.txt for make_gif.py.
STEPS   1 reader (d) -> merged grid -> shell surface + feature edges; d resampled on a grid around the crack
        2 shell, edges, damage volume (ray cast) + colour bar
        3 labels   4 camera   5 time frames   6 orbit frames
USAGE   pvpython render_c1_gif_volume.py EXODUS OUT_DIR [--test] [--orbit-only]
        (--test: final state from the main camera and from the right end only)
"""
import math
import os
import sys

from paraview.simple import (ColorBy, CreateView, ExtractSurface, FeatureEdges, GetColorTransferFunction,
                             GetOpacityTransferFunction, GetScalarBar, MergeBlocks, OpenDataFile, ResampleToImage,
                             SaveScreenshot, Show, Text)

T_FROM = 18.0
IMAGE = [1200, 760]
SHELL_OPACITY = 0.18
SHELL_RGB = [0.72, 0.74, 0.78]
# damage colour (d, r, g, b) and opacity (d, opacity, midpoint, sharpness); OPACITY_UNIT_M = path length (m) over
# which a point reaches the listed opacity, set near the broken band's width so the band is solid seen face-on
DAMAGE_RGB = [0.0, 0.86, 0.86, 0.86,  0.35, 0.98, 0.72, 0.55,  0.7, 0.89, 0.29, 0.20,  1.0, 0.55, 0.0, 0.05]
DAMAGE_OPACITY = [0.0, 0.0, 0.5, 0.0,  0.08, 0.0, 0.5, 0.0,  0.4, 0.25, 0.5, 0.0,  0.8, 0.75, 0.5, 0.0,  1.0, 1.0, 0.5, 0.0]
OPACITY_UNIT_M = 2.0e-4
# the volume is sampled on a regular grid around the crack line (the ray-cast mapper mixes with the translucent
# shell; the unstructured-grid mapper cannot): x 12-16 mm, full height and depth, 0.05 mm spacing
SAMPLE_BOUNDS_M = [0.012, 0.016, 0.0, 0.008, 0.0, 0.008]
SAMPLE_DIMS = [81, 161, 161]
TITLE = "C1   (c) two holes, not connected: xL 3.00 / xR 2.00 mm, r 1.0 mm"
KEY = "colour: damage d inside the specimen\ngrey: specimen, translucent"
CAMERA_DIR = (0.80, 0.42, 0.85)
CAMERA_DIST_M = 0.047
LOAD_MM_PER_T = 0.1
ORBIT_DEG, ORBIT_FRAMES = 180.0, 60
TIME_MS, HOLD_FINAL_TIME_MS, ORBIT_MS, HOLD_END_MS = 100, 1200, 80, 1500


def rotate_about_y(direction, degrees):
    """Camera direction turned about the vertical (y) axis; positive angles swing it from +z toward +x first."""
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
    damage = ResampleToImage(Input=merged)
    damage.UseInputBounds = 0
    damage.SamplingBounds = SAMPLE_BOUNDS_M
    damage.SamplingDimensions = SAMPLE_DIMS
    return reader, damage, shell, edges


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


def labels(view):
    for text, size, where in ((TITLE, 17, "Upper Left Corner"), (KEY, 15, "Lower Right Corner")):
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


def main(path, out_dir, test, orbit_only):
    os.makedirs(out_dir, exist_ok=True)
    reader, damage, shell, edges = build(path)
    times = [t for t in reader.TimestepValues if t >= T_FROM - 1e-9]
    view = CreateView("RenderView")
    view.ViewSize = IMAGE
    view.OrientationAxesVisibility = 1
    for name, value in (("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"),
                        ("Background", [1.0, 1.0, 1.0]), ("UseDepthPeeling", 1), ("MaximumNumberOfPeels", 8),
                        ("UseDepthPeelingForVolumes", 1)):
        if view.GetProperty(name) is not None:
            setattr(view, name, value)
    style(view, damage, shell, edges)
    clock_src = labels(view)
    reader.UpdatePipeline(times[-1])
    bounds = reader.GetDataInformation().GetBounds()
    final = "t = %.1f (final)    load-line displacement %.2f mm" % (times[-1], LOAD_MM_PER_T * times[-1])
    if test:
        view.ViewTime = times[-1]
        clock_src.Text = final
        for name, direction in (("main", CAMERA_DIR), ("rightend", rotate_about_y(CAMERA_DIR, 45.0))):
            camera(view, bounds, direction)
            SaveScreenshot(os.path.join(out_dir, "test_%s.png" % name), view, ImageResolution=IMAGE)
        print("wrote test views")
        return
    camera(view, bounds, CAMERA_DIR)
    if not orbit_only:
        for k, t in enumerate(times):
            view.ViewTime = t
            clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (t, LOAD_MM_PER_T * t)
            SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % k), view, ImageResolution=IMAGE)
        print("wrote %d time frames, t = %.1f .. %.1f" % (len(times), times[0], times[-1]))
    view.ViewTime = times[-1]
    clock_src.Text = final
    for j in range(1, ORBIT_FRAMES + 1):
        camera(view, bounds, rotate_about_y(CAMERA_DIR, ORBIT_DEG * j / ORBIT_FRAMES))
        SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % (len(times) - 1 + j)), view, ImageResolution=IMAGE)
    print("wrote %d orbit frames (%.0f degrees)" % (ORBIT_FRAMES, ORBIT_DEG))
    durations = [TIME_MS] * (len(times) - 1) + [HOLD_FINAL_TIME_MS] + [ORBIT_MS] * (ORBIT_FRAMES - 1) + [HOLD_END_MS]
    with open(os.path.join(out_dir, "durations.txt"), "w") as fh:
        fh.write("\n".join(str(ms) for ms in durations) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--test" in sys.argv, "--orbit-only" in sys.argv)
