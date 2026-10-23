#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Vista: timing and validation run of the homogeneous-Gc, eta = 0.45, absorbing-face case on the
# production mesh (spe_cube_h0p3.msh, 14.6 M tets, 2.4 M nodes). Splits the mesh for the rank count of this job if the split
# does not exist yet, then runs the first 3 us of pulse 1 (about 50 steps: 0.3 us rise at dt 3e-8, then dt 1e-7) and writes
# the Exodus frames at INITIAL and FINAL, so that the step cost, the memory per node and the time of one full-mesh Exodus write
# are measured before the 20-pulse chain is launched. Same physics as the 2D att_homog_eta0p45_absorb case: AT1, spectral
# split, Gc = 8 l ft^2/(3E) = 2721 J/m^2 (l = 0.6 mm), P_peak = 214 MPa at the gap, zeta_pulse = 1e-8 s, quench after 25 us.
# Submit from this folder:  sbatch submit_vista_3d_timing.sbatch      (8 nodes x 144 = 1152 ranks, 1.5 h, about 12 SUs)
#   other sizes:  sbatch -N 4 -n 576 submit_vista_3d_timing.sbatch ;  sbatch -N 16 -n 2304 submit_vista_3d_timing.sbatch
#   (the split, the output folder timing_<N>n and the perf graph at step 25 and at the end follow the job size)
# =============================================================================
#SBATCH -J spe3d_timing
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 1152
#SBATCH -t 01:30:00
#SBATCH -A EAR20006
#SBATCH -o spe3d_timing.o%j
#SBATCH -e spe3d_timing.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_3d
cd $CODE || exit 1
NR=$SLURM_NTASKS
MESH=../mesh/spe_cube_h0p3.msh
SPLIT=../mesh/spe_cube_h0p3_split${NR}.cpr       # one split per rank count (--use-split needs the matching count)
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20"
OUT=cases_vista/timing_${SLURM_NNODES}n
mkdir -p $OUT
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  ibrun -n 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > $OUT/split.log 2>&1
  echo "[$(date)] split exit $? ; $(du -sh $SPLIT 2>/dev/null)"
fi
echo "[$(date)] timing run on $SLURM_NNODES nodes, $NR ranks"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT \
      Executioner/end_time=3e-6 Outputs/exodus/time_step_interval=1000 Outputs/perf/time_step_interval=25 \
      Outputs/file_base=$OUT/out > $OUT/run.log 2>&1
echo "[$(date)] timing run exit $?"
grep -E "Time Step|Solve Converged|Finished Executing" $OUT/run.log | tail -5
