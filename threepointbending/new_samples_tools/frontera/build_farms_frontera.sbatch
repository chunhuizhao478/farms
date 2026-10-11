#!/bin/bash
#SBATCH -J build_farms_cdms        # Job name
#SBATCH -o build_farms_cdms.o%j    # Name of stdout output file
#SBATCH -e build_farms_cdms.e%j    # Name of stderr error file
#SBATCH -p development             # Queue (partition) name
#SBATCH -N 1                       # Total # of nodes
#SBATCH -n 1                       # Total # of mpi tasks
#SBATCH -t 02:00:00                # Run time (hh:mm:ss)
#SBATCH --mail-type=all            # Send email at begin and end of job
#SBATCH -A EAR20006                # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

# Rebuild ../../../farms-opt (the executable every jobs_* file calls as ../../farms-opt) after the
# scratch purge of 2026-10 removed every file of farms_cdms_04192026.  Builds against the MOOSE
# that is still on $WORK, with that build's own compiler and MPI (its petscvariables: gcc 9.1.0 +
# mvapich2-x 2.3).  NOT SUBMITTED YET -- needs the user's go (README.md, Phase 6).
# Submit from this folder after the branch is checked out into PROJECT_BASE.

#=== USER PARAMETERS - MODIFY THIS ===
PROJECT_BASE="/scratch2/10024/zhaochun/projects/farms_cdms_04192026"
export MOOSE_DIR="$WORK/projects/moose94b6528"   # fork commit 94b6528, = the Mac build's MOOSE
#=====================================

# Modules first: Lmod scripts do not survive `set -u`
module swap intel gcc/9.1.0
module swap impi mvapich2-x/2.3
module list
set -eo pipefail

export METHOD=opt
cd "$PROJECT_BASE"
echo "building $(git rev-parse --short HEAD 2>/dev/null || echo unknown) against $MOOSE_DIR"
make -j 56
./farms-opt --version
