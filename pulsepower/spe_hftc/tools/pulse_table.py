"""Per-pulse extents and energies of a Vista case CSV, with --recover resets handled.
usage: python3 pulse_table.py [--json=out.json] label=path [label=path ...]
A recover appends rows that restart at the checkpoint time: rows of the old segment after that
time are orphaned and dropped; dissipated_energy_total restarts from zero at each recover, so the
old-segment value at the checkpoint time is added as an offset (external_work and damping_work
continue across a recover and need no offset)."""
import csv, sys, json
GC = 2721.0  # J/m^2
PAPER = {5: 7.97, 11: 13.90, 15: 15.66, 20: 20.45}
def load(path):
    R = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(path))]
    segs = [[R[0]]]
    for a, b in zip(R, R[1:]):
        if b["time"] <= a["time"]: segs.append([b])
        else: segs[-1].append(b)
    clean, offset, resets = [], 0.0, []
    for i, seg in enumerate(segs):
        if i + 1 < len(segs):
            t_cp = segs[i + 1][0]["time"]
            keep = [r for r in seg if r["time"] <= t_cp + 1e-12]
            pre = keep[-1]["dissipated_energy_total"] if keep else 0.0
            resets.append((t_cp * 1e6, pre + offset, seg[-1]["time"] * 1e6))
        else:
            keep = seg
        for r in keep:
            r = dict(r); r["crackE"] = r["dissipated_energy_total"] + offset; clean.append(r)
        if i + 1 < len(segs): offset += pre
    return clean, resets
def windows(rows, key="fracture_extent"):
    """(pulse, value, complete, t_us) at the end of each 60 us window; the last entry may be the window in progress."""
    out = []
    for k in range(1, 21):
        win = [r for r in rows if (k - 1) * 60.0 < r["time"] * 1e6 <= k * 60.0]
        if not win: break
        complete = win[-1]["time"] * 1e6 >= k * 60.0 - 1.0
        out.append((k, win[-1][key], complete, win[-1]["time"] * 1e6))
        if not complete: break
    return out
def main(args):
    res, jpath = {}, None
    for arg in args:
        if arg.startswith("--json="): jpath = arg[7:]; continue
        label, path = arg.split("=", 1)
        rows, resets = load(path); last = rows[-1]
        ext = windows(rows); ce = windows(rows, "crackE")
        print(f"== {label}: {len(rows)} rows after cleaning, last t = {last['time']*1e6:.1f} us")
        for t_cp, pre, t_old in resets:
            print(f"   recover at {t_cp:.1f} us (old segment ran to {t_old:.1f} us, {t_old - t_cp:.1f} us orphaned); crack energy before the reset {pre:.1f} J/m")
        damp = last.get("damping_work", 0.0)   # no dashpots (free faces): no column
        print(f"   crack energy {last['crackE']:.1f} J/m (eq. crack length {last['crackE']/GC*1e3:.0f} mm), pressure work {-last['external_work']:.0f} J/m, "
              f"dashpot {damp:.0f} J/m, crack/pressure {100*last['crackE']/(-last['external_work']):.1f} %")
        print("   extent by pulse (mm):", ", ".join(f"{k}: {e*1e3:.2f}" + ("" if c else f" (partial, {t:.0f} us)") for k, e, c, t in ext))
        incs = [(ext[i][0], (ext[i][1] - ext[i-1][1]) * 1e3) for i in range(1, len(ext)) if ext[i][2]]
        print("   increment per pulse (mm):", ", ".join(f"{k}: {d:.2f}" for k, d in incs))
        print("   crack energy at window ends (J/m):", ", ".join(f"{k}: {e:.0f}" for k, e, c, t in ce))
        print("   vs paper (model, measured):", {k: (round(e * 1e3, 2), PAPER[k]) for k, e, c, t in ext if k in PAPER and c})
        res[label] = {"last_t_us": last["time"] * 1e6, "crackE": last["crackE"], "ext_work": -last["external_work"], "damp": last.get("damping_work", 0.0),
                      "extent_mm": [(k, round(e * 1e3, 3), c) for k, e, c, t in ext], "resets": resets}
    if jpath: json.dump(res, open(jpath, "w"), indent=1)
if __name__ == "__main__": main(sys.argv[1:])
