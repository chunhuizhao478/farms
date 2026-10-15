#!/bin/bash
# Remote, read-only: for every Frontera restart log, the last simulation time and the hard steps (time steps whose
# MultiApp fixed-point loop took >= 5 iterations, the sign of a crack onset or jump): count, first and last t, and
# the most iterations of any step.  Main-app "Time Step" lines only (the fracture0 sub-app lines are skipped).
# Run as: ssh HOST 'bash -s' < frontera_hard_steps.sh
P=/scratch2/10024/zhaochun/projects/farms_cdms_04192026/threepointbending
echo __B__
for f in $P/jobs_*_recover/3dptrbend_*_recover.o7990*; do
  [ -f "$f" ] || continue
  awk -v id="${f##*.o}" '
    /Time Step/ && !/fracture0/ { if (t != "" && kmax >= 4) { n++; if (first == "") first = t; last = t }
                                  if (kmax + 1 > most) most = kmax + 1
                                  s = $0; sub(/.*time = /, "", s); t = s + 0; kmax = -1; next }
    /Beginning fixed point iteration/ { s = $0; sub(/.*Beginning fixed point iteration /, "", s); k = s + 0; if (k > kmax) kmax = k }
    END { printf "H %s t=%s hard=%d first=%s last=%s most_its=%d\n", id, t, n, (first == "" ? "-" : first), (last == "" ? "-" : last), most }
  ' "$f"
done
echo __E__
