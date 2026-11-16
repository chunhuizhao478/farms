#!/bin/bash
# Login-node driver for the Vista 3D cube run: submits submit_vista_3d_eta0p45.sbatch, watches run.log, and resubmits with
# RECOVER=1 until farms finishes. A job whose run.log is silent for STALL_S seconds while running is cancelled (the
# hypre dead-lock of 2026-11-13/15 shows up as a running job with a silent log) and the run continues from the newest
# checkpoint. Stops after MAXSEG segments or after two segments in a row that end in under 15 min.
#   nohup ./run_chain_vista.sh [RECOVER_FIRST=1|0] > chain.log 2>&1 &
# FIRST_FROM=<dir>/<step> makes the first segment recover from that explicit checkpoint. SBATCH_OPTS adds sbatch options.
HERE=$(cd $(dirname $0) && pwd); cd $HERE
FIRST_FROM=${FIRST_FROM:-}
REC=${1:-1}; MAXSEG=${MAXSEG:-8}; STALL_S=${STALL_S:-1500}; EXTRA=${EXTRA:-dt_quiet=5e-7}
L=$SCRATCH/spe3d_vista/cube_homog_eta0p45_absorb/run.log; short=0
for seg in $(seq 1 $MAXSEG); do
  mkdir -p $SCRATCH/spe3d_vista/slurm_logs
  J=$(sbatch -o $SCRATCH/spe3d_vista/slurm_logs/%x.o%j -e $SCRATCH/spe3d_vista/slurm_logs/%x.e%j ${SBATCH_OPTS:-} --export=ALL,RECOVER=$REC,RECOVER_FROM=$FIRST_FROM,EXTRA="$EXTRA" submit_vista_3d_eta0p45.sbatch 2>&1 | grep -o "Submitted batch job [0-9]*" | grep -o "[0-9]*$")
  [ -n "$J" ] || { echo "[$(date)] sbatch failed, stop"; exit 1; }
  echo "[$(date)] segment $seg: job $J (RECOVER=$REC FROM=$FIRST_FROM)"; REC=1; FIRST_FROM=""
  while :; do
    sleep 300
    st=$(squeue -h -j $J -o %t 2>/dev/null); rc=$?
    [ $rc -ne 0 ] && [ -z "$st" ] && sacct -n -X -j $J -o State > /dev/null 2>&1 || true
    if [ -z "$st" ]; then sleep 20; [ -z "$(squeue -h -j $J -o %t 2>/dev/null)" ] && break; continue; fi
    if [ "$st" = R ] && [ -f $L ]; then
      age=$(( $(date +%s) - $(stat -c %Y $L) )); run=$(squeue -h -j $J -o %M)
      echo "[$(date)] $J running $run, $(grep -a '^Time Step' $L | tail -1), log age $age s"
      if [ $age -gt $STALL_S ] && [[ "$run" == *:*:* ]]; then echo "[$(date)] STALL: cancelling $J"; scancel $J; fi
    fi
  done
  s=$(sacct -n -X -j $J -o State,Elapsed | tr -s ' '); echo "[$(date)] job $J ended:$s"
  case "$s" in *COMPLETED*) echo "[$(date)] farms finished, chain done"; exit 0;; esac
  el=$(sacct -n -X -j $J -o ElapsedRaw | tr -d ' '); if [ "${el:-0}" -lt 900 ]; then short=$((short+1)); else short=0; fi
  [ $short -ge 2 ] && { echo "[$(date)] two short segments in a row, stop"; exit 1; }
done
echo "[$(date)] MAXSEG reached, stop"
