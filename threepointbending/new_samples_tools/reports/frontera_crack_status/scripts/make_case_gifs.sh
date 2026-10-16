#!/bin/bash
# make_case_gifs.sh -- the approved crack GIF (render_case_gif.py + make_gif.py) for the listed case IDs, two at a time.
#   Each GIF goes next to its case's solution file:  $DST/<case>/<ID>_crack_plane_xray_damage_t<from>-<to>.gif
#   Start time: about 3.5 before the case's crack onset (first step with >= 5 fixed-point iterations in the Frontera
#   log), rounded down to 0.5, so every GIF shows the build-up and the jump. Frames + per-case logs go to $FRAMES_ROOT.
# USAGE   bash make_case_gifs.sh C2 D1 D2 D3
#         (run under a memory kill-switch on the 36 GB laptop, e.g. run_guarded.sh <logfile> bash make_case_gifs.sh ...;
#          8 cases took 4 min with a peak of 2.9 GB for both lanes)
HERE=$(cd "$(dirname "$0")" && pwd)                    # render_case_gif.py and make_gif.py sit next to this script
DST=${DST:-"/Volumes/Extreme SSD/threepointbending/new_sample"}
FRAMES_ROOT=${FRAMES_ROOT:-${TMPDIR:-/tmp}/tpb_gif_frames}
PV=${PV:-/Applications/ParaView-6.2.0-RC1.app/Contents/bin/pvpython}
PY=${PY:-/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python}
mkdir -p "$FRAMES_ROOT"

spec() {   # ID -> "case folder|title|t_from"
  case $1 in
    C1) echo "case_2holes_asym_3d_xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm|C1   (c) two holes, not connected: xL 3.00 / xR 2.00 mm, r 1.0 mm|18.0" ;;
    C2) echo "case_2holes_asym_3d_xL3.00mm_xR2.50mm_y3.6mm_radi1.0mm|C2   (c) two holes, not connected: xL 3.00 / xR 2.50 mm, r 1.0 mm|18.5" ;;
    C3) echo "case_2holes_asym_3d_xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm|C3   (c) two holes, not connected: xL 3.50 / xR 2.50 mm, r 1.0 mm|19.5" ;;
    D1) echo "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.00mm_y3.6mm_overlap|D1   (d) two connected holes: DL 1.5 / DR 2.0 mm, x 1.00 mm, overlap|20.5" ;;
    D2) echo "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline|D2   (d) two connected holes: DL 1.5 / DR 2.0 mm, x 1.50 mm, downline slot|24.0" ;;
    D3) echo "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_upline|D3   (d) two connected holes: DL 1.5 / DR 2.0 mm, x 1.50 mm, upline slot|23.0" ;;
    D4) echo "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_downline|D4   (d) two connected holes: DL 1.5 / DR 2.0 mm, x 2.00 mm, downline slot|24.0" ;;
    D5) echo "case_2holes_connected_3d_DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_midline|D5   (d) two connected holes: DL 1.5 / DR 2.0 mm, x 2.00 mm, midline slot|25.0" ;;
    D6) echo "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x1.00mm_y3.6mm_overlap|D6   (d) two connected holes: DL 2.0 / DR 2.0 mm, x 1.00 mm, overlap|24.0" ;;
    D7) echo "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x1.50mm_y3.6mm_downline|D7   (d) two connected holes: DL 2.0 / DR 2.0 mm, x 1.50 mm, downline slot|22.0" ;;
    D8) echo "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_downline|D8   (d) two connected holes: DL 2.0 / DR 2.0 mm, x 2.00 mm, downline slot|24.0" ;;
    D9) echo "case_2holes_connected_3d_DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline|D9   (d) two connected holes: DL 2.0 / DR 2.0 mm, x 2.00 mm, upline slot|26.0" ;;
    E1) echo "case_sinewave_3d_D1.0mm_y1.6mm_amp2.0mm_centered|E1   (e) sine-wave void, D 1.0 mm, amplitude 2.0 mm, centred|19.0" ;;
    E2) echo "case_sinewave_3d_D1.0mm_y1.6mm_amp2.0mm_1ampoffset|E2   (e) sine-wave void, D 1.0 mm, amplitude 2.0 mm, offset one amplitude|19.0" ;;
    F1) echo "case_pillars_3d_welded_mismatch_2onN|F1   (f) pillar structures, welded mismatch 2onN|20.5" ;;
    F2) echo "case_pillars_3d_welded_mismatch_3onN|F2   (f) pillar structures, welded mismatch 3onN|21.0" ;;
  esac
}

run_case() {   # ID
  local s=$(spec $1) folder title t_from frames log
  [ -n "$s" ] || { echo "$1 unknown case ID"; return; }
  folder=${s%%|*}; s=${s#*|}; title=${s%%|*}; t_from=${s#*|}
  frames=$FRAMES_ROOT/gif_frames_$1
  log=$FRAMES_ROOT/gif_$1.log
  rm -rf "$frames"
  "$PV" --force-offscreen-rendering "$HERE/render_case_gif.py" "$DST/$folder/elasticity_exodus.e" "$frames" \
        --title "$title" --t-from $t_from > "$log" 2>&1 || { echo "$1 render FAILED (see $log)"; return; }
  "$PY" "$HERE/make_gif.py" "$frames" "$DST/$folder/$1_crack_plane_xray_damage_$(cat "$frames/range.txt").gif" >> "$log" 2>&1 \
     && echo "$1 $(tail -n 1 "$log")" || echo "$1 GIF FAILED (see $log)"
}

ids=("$@")
lane() { for id in "$@"; do run_case $id; done; }
half=$(( (${#ids[@]} + 1) / 2 ))
lane "${ids[@]:0:$half}" &
lane "${ids[@]:$half}" &
wait
echo GIFS_DONE
