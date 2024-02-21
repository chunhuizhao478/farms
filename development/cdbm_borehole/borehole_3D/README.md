3D Borehole Test 

We setup the case as in the paper:

Borehole Breakouts Induced in Arkosic Sandstones and a Discrete Element Analysis,  H. Lee1 • T. Moon2 • B. C. Haimson3

Test Feb 20
MOOSE Version:           git commit 2bbd8e54e1 on 2024-02-13
LibMesh Version:         
PETSc Version:           3.20.3
SLEPc Version:           3.20.1
Current Time:            Tue Feb 20 19:42:57 2024
Executable Timestamp:    Tue Feb 20 16:56:47 2024

Parallelism:
  Num Processors:          8
  Num Threads:             1

Mesh: 
  Parallel Type:           distributed (pre-split)
  Mesh Dimension:          3
  Spatial Dimension:       3
  Nodes:                   
    Total:                 34314
    Local:                 4741
    Min/Max/Avg:           3859/4741/4289
  Elems:                   
    Total:                 185195
    Local:                 23150
    Min/Max/Avg:           23149/23150/23149
  Num Subdomains:          1
  Num Partitions:          1
  Partitioner:             parmetis

Nonlinear System:
  Num DOFs:                102942
  Num Local DOFs:          14223
  Variables:               { "disp_x" "disp_y" "disp_z" } 
  Finite Element Types:    "LAGRANGE" 
  Approximation Orders:    "FIRST" 

Auxiliary System:
  Num DOFs:                13457671
  Num Local DOFs:          1686319
  Variables:               { "disp_cdbm_x" "disp_cdbm_y" "disp_cdbm_z" "vel_cdbm_x" "vel_cdbm_y" "vel_cdbm_z" 
                             } { "mu_s" "mu_d" } "ini_shear_stress" { "tangent_jump_rate" "tria_area_aux" 
                             } { "nodal_area" "B_old" } { "xi_old" "I2_old" "mu_old" "lambda_old" "gamma_old" 
                             "alpha_in" "B_in" } { "alpha_in_dummy" "B_in_dummy" } { "alpha_grad_x" "alpha_grad_y" 
                             "alpha_grad_z" "stress_xx" "stress_yy" ... "sts_total_11" "sts_total_12" 
                             "sts_total_20" "sts_total_21" "sts_total_22" } 
  Finite Element Types:    "LAGRANGE" "MONOMIAL" "LAGRANGE" "MONOMIAL" "LAGRANGE" "MONOMIAL" "MONOMIAL" 
                             "MONOMIAL" 
  Approximation Orders:    "FIRST" "CONSTANT" "FIRST" "CONSTANT" "FIRST" "CONSTANT" "FIRST" "CONSTANT" 
                             

Execution Information:
  Executioner:             Transient
  TimeStepper:             ConstantDT
  TimeIntegrator:          CentralDifference
  Solver Mode:             Linear

LEGACY MODES ENABLED:
 This application uses the legacy material output option: material properties are output only on TIMESTEP_END, not INITIAL. To remove this message, set 'use_legacy_material_output' to false in this application. If there are gold output files that contain material property output for which output occurs on INITIAL, then these will generate diffs due to zero values being stored, and these tests should be re-golded.

Currently Executing
  Performing Initial Setup.....
    Setting Up Materials
      Finished Computing Initial Material Values                                         [  5.69 s] [  713 MB]
    Finished Setting Up Materials                                                        [  5.69 s] [  713 MB]
sub_app0: Parallelism:
sub_app0:   Num Processors:          8
sub_app0:   Num Threads:             1
sub_app0: 
sub_app0: Mesh: 
sub_app0:   Parallel Type:           distributed (pre-split)
sub_app0:   Mesh Dimension:          3
sub_app0:   Spatial Dimension:       3
sub_app0:   Nodes:                   
sub_app0:     Total:                 34314
sub_app0:     Local:                 4741
sub_app0:     Min/Max/Avg:           3859/4741/4289
sub_app0:   Elems:                   
sub_app0:     Total:                 185195
sub_app0:     Local:                 23150
sub_app0:     Min/Max/Avg:           23149/23150/23149
sub_app0:   Num Subdomains:          1
sub_app0:   Num Partitions:          1
sub_app0:   Partitioner:             parmetis
sub_app0: 
sub_app0: Nonlinear System:
sub_app0:   Num DOFs:                1851950
sub_app0:   Num Local DOFs:          231500
sub_app0:   Variables:               { "alpha_sub" "B_sub" } { "alpha_sub_dummy" "B_sub_dummy" } 
sub_app0:   Finite Element Types:    "MONOMIAL" "MONOMIAL" 
sub_app0:   Approximation Orders:    "CONSTANT" "FIRST" 
sub_app0: 
sub_app0: Auxiliary System:
sub_app0:   Num DOFs:                5555850
sub_app0:   Num Local DOFs:          694500
sub_app0:   Variables:               { "alpha_old" "B_old" } { "alpha_old_dummy" "B_old_dummy" } { "xi_old" "I2_old" 
sub_app0:                              "mu_old" "lambda_old" "gamma_old" "alpha_checked" "B_checked" } { "alpha_checked_dummy" 
sub_app0:                              "B_checked_dummy" } { "alpha_grad_x_sub" "alpha_grad_y_sub" "alpha_grad_z_sub" 
sub_app0:                              "mechanical_strain_rate_sub" "track_Cd" } 
sub_app0:   Finite Element Types:    "MONOMIAL" "MONOMIAL" "MONOMIAL" "MONOMIAL" "MONOMIAL" 
sub_app0:   Approximation Orders:    "CONSTANT" "FIRST" "CONSTANT" "FIRST" "CONSTANT" 
sub_app0: 
sub_app0: Execution Information:
sub_app0:   Executioner:             Transient
sub_app0:   TimeStepper:             ConstantDT
sub_app0:   TimeIntegrator:          ExplicitSSPRungeKutta
sub_app0:   Solver Mode:             Linear
sub_app0: 
sub_app0: LEGACY MODES ENABLED:
sub_app0:  This application uses the legacy material output option: material properties are output only on TIMESTEP_END, not INITIAL. To remove this message, set 'use_legacy_material_output' to false in this application. If there are gold output files that contain material property output for which output occurs on INITIAL, then these will generate diffs due to zero values being stored, and these tests should be re-golded.
sub_app0: 
sub_app0: 
sub_app0: Time Step 0, time = 0
  Still Performing Initial Setup
    Computing Initial Material Properties                                                [  6.64 s] [  649 MB]
  Finished Performing Initial Setup                                                      [ 53.53 s] [  609 MB]
