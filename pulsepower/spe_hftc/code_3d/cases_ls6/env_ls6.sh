# Environment for building/running MOOSE + farms on TACC Lonestar6 (AMD EPYC 7763, 128 cores per node).
# gcc 13.2 + Intel MPI 2021.18 (no OpenMPI under gcc on LS6). MOOSE at the commit of the Vista build (6e0969c3),
# so that the Vista checkpoints of the 3D cube run can be recovered here.
module purge > /dev/null 2>&1
module load gcc/13.2.0 impi/21.18 cmake/4.1.1 > /dev/null 2>&1
export SPE3D_ROOT=$WORK/spe3d
export MOOSE_DIR=$SPE3D_ROOT/moose-6e0969c3
export FARMS_DIR=$SPE3D_ROOT/farms_cdms
export METHODS=opt
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77
# The MOOSE build needs Python >= 3.9 (python3 and python3-config); the system python3 is 3.6 and the python modules do
# not load under gcc/13.2.0, so the 3.12 installation is put on PATH directly. pyyaml/packaging/jinja2 for the MOOSE
# scripts are in $SPE3D_ROOT/pydeps (pip install --target, done once on the login node).
export PATH=/opt/apps/python/3.12.11/bin:$PATH
export PYTHONPATH=$SPE3D_ROOT/pydeps${PYTHONPATH:+:$PYTHONPATH}
