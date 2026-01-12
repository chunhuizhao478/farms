# TPV2052D Test

## CZM Formulation

This guide will help you run the TPV2052D test using the CZM formulation.

### Steps to Run the Code:

1. **Run the simulation**:  
   Use the `mpirun` command to execute the test with 8 processors:
   ```bash
   mpirun -np 8 ./farms-opt -i examples/benchmark_tpv2052D/czm/tpv2052D_quad.i
2. **Post-process the results**
   Run the first post-processing script to process the output:
   (node_per_elem = 4 for quad4, node_per_elem = 3 for tria3)
   ```bash
   python3 postprocess/tpv2052d/postmain.py 
3. **Compare the data** 
   Execute the data comparison script to compare the results:
   ```bash
   python3 postprocess/tpv2052d/datacompare.py

### The benchmark results as saved in the folder:

1. postprocess/tpv2052d/outputs_quad_200m_saved/
2. postprocess/tpv2052d/outputs_tria_200m_saved/