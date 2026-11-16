#!/bin/bash
# Login-node driver for the Vista 3D cube run: keeps submit_vista_3d_eta0p45.sbatch going until farms finishes.
# - A job whose run.log is silent for STALL_S seconds while running is cancelled (the hypre dead-lock of 2026-11-13/15
#   shows up as a running job with a silent log) and the run continues from the newest checkpoint.
# - As soon as a segment starts running, the next one is queued with --dependency=afterany, so it gains queue priority
#   while the current one runs (queue waits were 4 to 10 h); it is cancelled if farms finishes.
# - Stops after MAXSEG segments or after two segments in a row that end in under 15 min.
#   nohup ./run_chain_vista.sh [RECOVER_FIRST=1|0] > chain.log 2>&1 &
# Env: ADOPT=<jobid> takes over an already submitted job as the first segment; SBATCH_OPTS adds sbatch options (e.g. -t);
# EXTRA are the input overrides; FIRST_FROM=<dir>/<step> makes the first submitted segment recover from that checkpoint.
HERE=$(cd $(dirname $0) && pwd); cd $HERE

OUTROOT=$SCRATCH/spe3d_vista
# one driver at a time: a lock folder on scratch (removed on exit)
LOCK=$OUTROOT/chain.lockdir; mkdir -p $OUTROOT/slurm_logs; mkdir $LOCK 2>/dev/null || { echo "[$(date)] another chain driver holds $LOCK, exit"; exit 0; }
trap "rmdir $LOCK" EXIT; trap "exit 1" TERM INT
REC=${1:-1}; MAXSEG=${MAXSEG:-12}; STALL_S=${STALL_S:-1500}; EXTRA=${EXTRA:-dt_quiet=5e-7}
L=$OUTROOT/cube_homog_eta0p45_absorb/run.log; short=0; J=${ADOPT:-}; NEXT=""
submit() { # $1 = RECOVER flag, $2 = explicit checkpoint or empty, $3 = extra sbatch options
  sbatch -o $OUTROOT/slurm_logs/%x.o%j -e $OUTROOT/slurm_logs/%x.e%j ${SBATCH_OPTS:-} $3 --export=ALL,RECOVER=$1,RECOVER_FROM=$2,EXTRA="$EXTRA" submit_vista_3d_eta0p45.sbatch 2>&1 | grep -o "Submitted batch job [0-9]*" | grep -o "[0-9]*$"
}
for seg in $(seq 1 $MAXSEG); do
  if [ -z "$J" ]; then J=$(submit $REC "${FIRST_FROM:-}" ""); [ -n "$J" ] || { echo "[$(date)] sbatch failed, stop"; exit 1; }; fi
  echo "[$(date)] segment $seg: job $J"; REC=1; FIRST_FROM=""
  while :; do
    sleep 300
    st=$(squeue -h -j $J -o %t 2>/dev/null)
    if [ -z "$st" ]; then sleep 30; [ -z "$(squeue -h -j $J -o %t 2>/dev/null)" ] && break; continue; fi
    if [ "$st" = R ]; then
      if [ -z "$NEXT" ]; then NEXT=$(submit 1 "" "--dependency=afterany:$J"); echo "[$(date)] next segment queued behind $J: job ${NEXT:-FAILED}"; fi
      if [ -f $L ]; then
        age=$(( $(date +%s) - $(stat -c %Y $L) )); run=$(squeue -h -j $J -o %M)
        echo "[$(date)] $J running $run, $(grep -a '^Time Step' $L | tail -1), log age $age s"
        if [ $age -gt $STALL_S ] && [[ "$run" == *:*:* ]]; then echo "[$(date)] STALL: cancelling $J"; scancel $J; fi
      fi
    fi
  done
  s=$(sacct -n -X -j $J -o State,Elapsed | tr -s ' '); echo "[$(date)] job $J ended:$s"
  case "$s" in *COMPLETED*) [ -n "$NEXT" ] && scancel $NEXT; echo "[$(date)] farms finished, chain done"; exit 0;; esac
  el=$(sacct -n -X -j $J -o ElapsedRaw | tr -d ' '); if [ "${el:-0}" -lt 900 ]; then short=$((short+1)); else short=0; fi
  if [ $short -ge 2 ]; then [ -n "$NEXT" ] && scancel $NEXT; echo "[$(date)] two short segments in a row, stop"; exit 1; fi
  J=$NEXT; NEXT=""
done
echo "[$(date)] MAXSEG reached, stop"
