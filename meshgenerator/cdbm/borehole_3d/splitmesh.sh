mpirun -np 8 ./farms-opt -i meshgenerator/cdbm/borehole_3d/meshcheck.i --split-mesh 8 --split-file meshgenerator/cdbm/borehole_3d/foo.cpr
mpirun -np 8 ./farms-opt -i meshgenerator/cdbm/borehole_3d/meshcheck.i --use-split --split-file meshgenerator/cdbm/borehole_3d/foo.cpr