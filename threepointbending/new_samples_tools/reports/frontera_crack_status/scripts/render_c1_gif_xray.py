"""
render_c1_gif_xray.py -- frames of the crack plane inside a translucent 3D specimen, case C1, for a GIF (pvpython).

WHAT    the beam's outer surface drawn translucent (light grey, SHELL_OPACITY) with dark feature edges, the crack as
        the opaque iso-surface d = D_CRACK (dark red) and the damage zone as a faint iso-surface d = D_ZONE (orange),
        so the crack plane is seen through the specimen, back faces included. A title, a colour key and a per-frame
        label with the simulation time t and the load-line displacement (0.1 mm per unit t). One PNG per output time
        from T_FROM to the end.
STEPS   1 reader (d) -> merged grid -> shell surface, feature edges, two contours   2 styles + labels   3 camera
        4 frames
USAGE   pvpython render_c1_gif_xray.py EXODUS OUT_DIR [--test]   (--test: last time only, one PNG per camera in
        TEST_CAMERAS, to choose the angle)
"""
import os
import sys

from paraview.simple import (Contour, CreateView, ExtractSurface, FeatureEdges, MergeBlocks, OpenDataFile,
                             SaveScreenshot, Show, Text, ColorBy)

T_FROM = 18.0
IMAGE = [1200, 760]
D_CRACK = 0.95                   # crack surface (fully broken)
D_ZONE = 0.5                     # damage zone
SHELL_OPACITY = 0.18
SHELL_RGB = [0.72, 0.74, 0.78]
CRACK_RGB = [0.78, 0.06, 0.08]
ZONE_RGB = [0.98, 0.55, 0.20]
ZONE_OPACITY = 0.18              # faint, so the crack under it stays red
TITLE = "C1   (c) two holes, not connected: xL 3.00 / xR 2.00 mm, r 1.0 mm"
KEY = "red: crack surface, d = 0.95\norange: damage zone, d = 0.5\ngrey: specimen, translucent"
CAMERA_DIR = (0.80, 0.42, 0.85)  # chosen from TEST_CAMERAS: front face, right end and top, crack plane well open
CAMERA_DIST_M = 0.047
TEST_CAMERAS = {"A": (0.42, 0.34, 1.0), "B": (0.80, 0.42, 0.85), "C": (1.0, 0.55, 0.45)}
LOAD_MM_PER_T = 0.1


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
    crack = Contour(Input=merged, ContourBy=["POINTS", "d"], Isosurfaces=[D_CRACK])
    zone = Contour(Input=merged, ContourBy=["POINTS", "d"], Isosurfaces=[D_ZONE])
    return reader, shell, edges, crack, zone


def solid(view, source, rgb, opacity=1.0, specular=0.0):
    shown = Show(source, view)
    shown.ColorArrayName = ["POINTS", ""]      # solid colour (no array)
    shown.DiffuseColor = rgb
    shown.AmbientColor = rgb
    shown.Opacity = opacity
    shown.Specular = specular
    return shown


def style(view, shell, edges, crack, zone):
    solid(view, shell, SHELL_RGB, SHELL_OPACITY)
    lines = solid(view, edges, [0.2, 0.2, 0.2])
    lines.LineWidth = 1.5
    solid(view, zone, ZONE_RGB, ZONE_OPACITY)
    solid(view, crack, CRACK_RGB, 1.0)


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


def main(path, out_dir, test):
    os.makedirs(out_dir, exist_ok=True)
    reader, shell, edges, crack, zone = build(path)
    times = [t for t in reader.TimestepValues if t >= T_FROM - 1e-9]
    view = CreateView("RenderView")
    view.ViewSize = IMAGE
    view.OrientationAxesVisibility = 1
    for name, value in (("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"),
                        ("Background", [1.0, 1.0, 1.0]), ("EnableRenderOnInteraction", 0),
                        ("UseDepthPeeling", 1), ("MaximumNumberOfPeels", 8)):
        if view.GetProperty(name) is not None:
            setattr(view, name, value)
    style(view, shell, edges, crack, zone)
    clock_src = labels(view)
    reader.UpdatePipeline(times[-1])
    bounds = reader.GetDataInformation().GetBounds()
    if test:
        view.ViewTime = times[-1]
        clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (times[-1], LOAD_MM_PER_T * times[-1])
        for name, direction in TEST_CAMERAS.items():
            camera(view, bounds, direction)
            SaveScreenshot(os.path.join(out_dir, "test_camera_%s.png" % name), view, ImageResolution=IMAGE)
        print("wrote %d test views at t = %.1f" % (len(TEST_CAMERAS), times[-1]))
        return
    camera(view, bounds, CAMERA_DIR)
    for k, t in enumerate(times):
        view.ViewTime = t
        clock_src.Text = "t = %.1f    load-line displacement %.2f mm" % (t, LOAD_MM_PER_T * t)
        SaveScreenshot(os.path.join(out_dir, "frame_%03d.png" % k), view, ImageResolution=IMAGE)
    print("wrote %d frames, t = %.1f .. %.1f" % (len(times), times[0], times[-1]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--test" in sys.argv)
