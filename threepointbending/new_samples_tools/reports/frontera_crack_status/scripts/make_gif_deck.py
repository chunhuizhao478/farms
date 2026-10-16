"""
make_gif_deck.py -- PowerPoint deck of the crack GIFs: one 16:9 slide per sample, named, with its animated GIF.

WHAT    for every case whose folder under ROOT holds <ID>_crack_plane_xray_damage_t*.gif (campaign order C1-C3, D1-D9,
        E1-E2, F1-F2), one slide: the title names the sample (case ID, sample class, variant), a line below gives its
        dimensions and the sample (folder) name, and the GIF sits centred under them at the largest size that fits.
        GIFs are embedded byte-for-byte, so they stay animated in PowerPoint. The slide background takes the GIF's own
        background colour (the shared GIF palette stores white as 252,252,252), so the GIF shows no box edge. Text in
        Arial, which every viewer has (Office's Calibri is missing outside Office, e.g. in Quick Look). Cases without a
        GIF yet get no slide and are listed on the console.
STEPS   1 find each case's GIF   2 one slide per case: background, title, sub-lines, picture   3 save   4 checks
CHECKS  reopens the saved deck: one GIF picture per slide, content type image/gif, embedded bytes identical (md5) to the
        source GIF, picture inside the slide.
USAGE   python make_gif_deck.py [ROOT] [OUT.pptx]          (needs python-pptx)
"""
import glob
import hashlib
import os
import sys

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, MSO_AUTO_SIZE, PP_ALIGN
from pptx.util import Emu, Inches, Pt

ROOT = "/Volumes/Extreme SSD/threepointbending/new_sample"
OUT = os.path.join(ROOT, "new_samples_crack_damage.pptx")
GIF_PATTERN = "%s_crack_plane_xray_damage_t*.gif"

# (ID, title naming the sample, dimensions line, sample folder) in campaign order
CASES = [
    ("C1", "(c) two holes, not connected", "xL 3.00 mm / xR 2.00 mm, r 1.0 mm", "case_2holes_asym_3d_xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm"),
    ("C2", "(c) two holes, not connected", "xL 3.00 mm / xR 2.50 mm, r 1.0 mm", "case_2holes_asym_3d_xL3.00mm_xR2.50mm_y3.6mm_radi1.0mm"),
    ("C3", "(c) two holes, not connected", "xL 3.50 mm / xR 2.50 mm, r 1.0 mm", "case_2holes_asym_3d_xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm"),
    ("D1", "(d) two connected holes, overlap", "DL 1.5 mm / DR 2.0 mm, x 1.00 mm", "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.00mm_y3.6mm_overlap"),
    ("D2", "(d) two connected holes, downline slot", "DL 1.5 mm / DR 2.0 mm, x 1.50 mm", "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline"),
    ("D3", "(d) two connected holes, upline slot", "DL 1.5 mm / DR 2.0 mm, x 1.50 mm", "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_upline"),
    ("D4", "(d) two connected holes, downline slot", "DL 1.5 mm / DR 2.0 mm, x 2.00 mm", "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_downline"),
    ("D5", "(d) two connected holes, midline slot", "DL 1.5 mm / DR 2.0 mm, x 2.00 mm", "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_midline"),
    ("D6", "(d) two connected holes, overlap", "DL 2.0 mm / DR 2.0 mm, x 1.00 mm", "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x1.00mm_y3.6mm_overlap"),
    ("D7", "(d) two connected holes, downline slot", "DL 2.0 mm / DR 2.0 mm, x 1.50 mm", "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x1.50mm_y3.6mm_downline"),
    ("D8", "(d) two connected holes, downline slot", "DL 2.0 mm / DR 2.0 mm, x 2.00 mm", "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_downline"),
    ("D9", "(d) two connected holes, upline slot", "DL 2.0 mm / DR 2.0 mm, x 2.00 mm", "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline"),
    ("E1", "(e) sine-wave void, centred", "D 1.0 mm, amplitude 2.0 mm", "case_sinewave_3d_D1.0mm_y1.6mm_amp2.0mm_centered"),
    ("E2", "(e) sine-wave void, offset one amplitude", "D 1.0 mm, amplitude 2.0 mm", "case_sinewave_3d_D1.0mm_y1.6mm_amp2.0mm_1ampoffset"),
    ("F1", "(f) pillar structures, welded mismatch 2onN", "", "case_pillars_3d_welded_mismatch_2onN"),
    ("F2", "(f) pillar structures, welded mismatch 3onN", "", "case_pillars_3d_welded_mismatch_3onN"),
]

SLIDE_W, SLIDE_H = Emu(12192000), Emu(6858000)          # 16:9 widescreen, 13.333 x 7.5 in
MARGIN_X = Inches(0.45)
TITLE_TOP, TITLE_H, TITLE_PT = Inches(0.18), Inches(0.62), 28
SUB_TOP, SUB_H, SUB_PT, SAMPLE_PT = Inches(0.78), Inches(0.5), 16, 12
GIF_TOP, GIF_BOTTOM_MARGIN = Inches(1.3), Inches(0.1)
DARK, MID, LIGHT = RGBColor(0x1A, 0x1A, 0x1A), RGBColor(0x40, 0x40, 0x40), RGBColor(0x80, 0x80, 0x80)
FONT = "Arial"


def find_gifs(root):
    """(case, gif path) for every case with a GIF, in campaign order; the IDs without one."""
    found, missing = [], []
    for case in CASES:
        hits = sorted(glob.glob(os.path.join(root, case[3], GIF_PATTERN % case[0])))
        if hits:
            found.append((case, hits[-1]))
        else:
            missing.append(case[0])
    return found, missing


def gif_size(path):
    """Pixel width and height from the GIF logical-screen descriptor (bytes 6-9, little-endian)."""
    with open(path, "rb") as fh:
        head = fh.read(10)
    assert head[:6] in (b"GIF87a", b"GIF89a"), "%s is not a GIF" % path
    return int.from_bytes(head[6:8], "little"), int.from_bytes(head[8:10], "little")


def gif_background(path):
    """The GIF's background colour: its top-left pixel (the renders have a plain white background)."""
    return Image.open(path).convert("RGB").getpixel((2, 2))


def style_runs(paragraph, text, size_pt, colour, bold=False):
    run = paragraph.add_run()
    run.text = text
    run.font.name = FONT
    run.font.size = Pt(size_pt)
    run.font.bold = bold
    run.font.color.rgb = colour


def add_case_slide(prs, case, gif):
    cid, name, dims, folder = case
    slide = prs.slides.add_slide(prs.slide_layouts[5])                  # "Title Only": a real slide title
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = RGBColor(*gif_background(gif))
    title = slide.shapes.title
    title.left, title.top, title.width, title.height = MARGIN_X, TITLE_TOP, SLIDE_W - 2 * MARGIN_X, TITLE_H
    tf = title.text_frame
    tf.clear()
    tf.word_wrap = False
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.vertical_anchor = MSO_ANCHOR.TOP
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.paragraphs[0].alignment = PP_ALIGN.LEFT
    style_runs(tf.paragraphs[0], "%s   %s" % (cid, name), TITLE_PT, DARK, bold=True)

    sub = slide.shapes.add_textbox(MARGIN_X, SUB_TOP, SLIDE_W - 2 * MARGIN_X, SUB_H).text_frame
    sub.word_wrap = False
    sub.margin_left = sub.margin_right = sub.margin_top = sub.margin_bottom = 0
    first = sub.paragraphs[0]
    if dims:
        style_runs(first, dims, SUB_PT, MID)
        first = sub.add_paragraph()
    style_runs(first, "sample: " + folder, SAMPLE_PT, LIGHT)

    w_px, h_px = gif_size(gif)
    height = SLIDE_H - GIF_TOP - GIF_BOTTOM_MARGIN
    width = int(height * w_px / h_px)
    if width > SLIDE_W - 2 * MARGIN_X:                                   # width-limited (not with 1200 x 760 GIFs)
        width = SLIDE_W - 2 * MARGIN_X
        height = int(width * h_px / w_px)
    slide.shapes.add_picture(gif, int((SLIDE_W - width) / 2), GIF_TOP, width=width, height=height)


def check_deck(out, found):
    """Each slide holds exactly one picture: a GIF byte-identical to its source and inside the slide."""
    prs = Presentation(out)
    assert len(prs.slides) == len(found), "%d slides for %d GIFs" % (len(prs.slides), len(found))
    for slide, ((cid, *_), gif) in zip(prs.slides, found):
        pics = [s for s in slide.shapes if s.shape_type == 13]          # MSO_SHAPE_TYPE.PICTURE
        assert len(pics) == 1, "%s: %d pictures" % (cid, len(pics))
        img = pics[0].image
        assert img.content_type == "image/gif", "%s: %s" % (cid, img.content_type)
        assert hashlib.md5(img.blob).hexdigest() == hashlib.md5(open(gif, "rb").read()).hexdigest(), "%s: bytes differ" % cid
        p = pics[0]
        assert p.left >= 0 and p.top >= 0 and p.left + p.width <= SLIDE_W and p.top + p.height <= SLIDE_H, "%s: off slide" % cid
        assert slide.shapes.title.text_frame.text.startswith(cid), "%s: title" % cid
    return len(prs.slides)


def main(root, out):
    found, missing = find_gifs(root)
    prs = Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
    for case, gif in found:
        add_case_slide(prs, case, gif)
    partial = out + ".part.pptx"
    prs.save(partial)
    os.replace(partial, out)                   # swap in whole: a copy open in PowerPoint is never half-written
    n = check_deck(out, found)
    print("wrote %s: %d slides (%s), %.1f MB; checks passed (animated GIFs embedded byte-identical)"
          % (out, n, " ".join(c[0] for c, _ in found), os.path.getsize(out) / 1e6))
    if missing:
        print("no GIF yet (no slide): " + " ".join(missing))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT, sys.argv[2] if len(sys.argv) > 2 else OUT)
