(moose) andyz@wirelessprv-10-193-217-157 farms % bash development/cdbm_borehole/borehole_3D_testspeed/splitmesh.sh


*** Warning, This code is deprecated and will be removed in future versions:
Please update your main.C to adapt new main function in MOOSE framework, see'test/src/main.C in MOOSE as an example of moose::main()'. 



*** Info ***
TensorMechanics Action: selecting 'total small strain' formulation. Use `incremental = true` to select 'incremental small strain' instead.
Splitting 8 ways...
    - writing 1 files per process...
Mesh meta data written into "/Users/andyz/projects/farms/development/cdbm_borehole/borehole_3D/foo_main.cpr/meta_data_mesh.rd".
Finished Setting Up                                                                      [  3.94 s] [  351 MB]
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There are 3 unused database options. They are:
Option left: name:--split-file value: development/cdbm_borehole/borehole_3D/foo_main.cpr source: command line
Option left: name:--split-mesh value: 8 source: command line
Option left: name:-i value: development/cdbm_borehole/borehole_3D_testspeed/test_borehole_main.i source: command line


*** Warning, This code is deprecated and will be removed in future versions:
Please update your main.C to adapt new main function in MOOSE framework, see'test/src/main.C in MOOSE as an example of moose::main()'. 



*** Info ***
TensorMechanics Action: selecting 'total small strain' formulation. Use `incremental = true` to select 'incremental small strain' instead.
sub_app0: Initializingsub_app0: 
sub_app0:   Finished Initializing Equation Systems                                       [  0.34 s] [  233 MB]
sub_app0: Finished Initializing                                                          [  0.36 s] [  233 MB]
Currently Setting Up
  Finished Instantiating Sub-Apps                                                        [  1.34 s] [  233 MB]
The following total 143 aux variables:
  B
  I1
  I2
  Jacobian_mult_0000
  Jacobian_mult_0001
  Jacobian_mult_0002
  Jacobian_mult_0010
  Jacobian_mult_0011
  Jacobian_mult_0012
  Jacobian_mult_0020
  Jacobian_mult_0021
  Jacobian_mult_0022
  Jacobian_mult_0100
  Jacobian_mult_0101
  Jacobian_mult_0102
  Jacobian_mult_0110
  Jacobian_mult_0111
  Jacobian_mult_0112
  Jacobian_mult_0120
  Jacobian_mult_0121
  Jacobian_mult_0122
  Jacobian_mult_0200
  Jacobian_mult_0201
  Jacobian_mult_0202
  Jacobian_mult_0210
  Jacobian_mult_0211
  Jacobian_mult_0212
  Jacobian_mult_0220
  Jacobian_mult_0221
  Jacobian_mult_0222
  Jacobian_mult_1000
  Jacobian_mult_1001
  Jacobian_mult_1002
  Jacobian_mult_1010
  Jacobian_mult_1011
  Jacobian_mult_1012
  Jacobian_mult_1020
  Jacobian_mult_1021
  Jacobian_mult_1022
  Jacobian_mult_1100
  Jacobian_mult_1101
  Jacobian_mult_1102
  Jacobian_mult_1110
  Jacobian_mult_1111
  Jacobian_mult_1112
  Jacobian_mult_1120
  Jacobian_mult_1121
  Jacobian_mult_1122
  Jacobian_mult_1200
  Jacobian_mult_1201
  Jacobian_mult_1202
  Jacobian_mult_1210
  Jacobian_mult_1211
  Jacobian_mult_1212
  Jacobian_mult_1220
  Jacobian_mult_1221
  Jacobian_mult_1222
  Jacobian_mult_2000
  Jacobian_mult_2001
  Jacobian_mult_2002
  Jacobian_mult_2010
  Jacobian_mult_2011
  Jacobian_mult_2012
  Jacobian_mult_2020
  Jacobian_mult_2021
  Jacobian_mult_2022
  Jacobian_mult_2100
  Jacobian_mult_2101
  Jacobian_mult_2102
  Jacobian_mult_2110
  Jacobian_mult_2111
  Jacobian_mult_2112
  Jacobian_mult_2120
  Jacobian_mult_2121
  Jacobian_mult_2122
  Jacobian_mult_2200
  Jacobian_mult_2201
  Jacobian_mult_2202
  Jacobian_mult_2210
  Jacobian_mult_2211
  Jacobian_mult_2212
  Jacobian_mult_2220
  Jacobian_mult_2221
  Jacobian_mult_2222
  alpha_damagedvar
  elastic_strain_00
  elastic_strain_01
  elastic_strain_02
  elastic_strain_10
  elastic_strain_11
  elastic_strain_12
  elastic_strain_20
  elastic_strain_21
  elastic_strain_22
  eps_e_00
  eps_e_01
  eps_e_02
  eps_e_10
  eps_e_11
  eps_e_12
  eps_e_20
  eps_e_21
  eps_e_22
  eps_p_00
  eps_p_01
  eps_p_02
  eps_p_10
  eps_p_11
  eps_p_12
  eps_p_20
  eps_p_21
  eps_p_22
  eps_total_00
  eps_total_01
  eps_total_02
  eps_total_10
  eps_total_11
  eps_total_12
  eps_total_20
  eps_total_21
  eps_total_22
  gamma_damaged
  lambda
  shear_modulus
  stress_00
  stress_01
  stress_02
  stress_10
  stress_11
  stress_12
  stress_20
  stress_21
  stress_22
  sts_total_00
  sts_total_01
  sts_total_02
  sts_total_10
  sts_total_11
  sts_total_12
  sts_total_20
  sts_total_21
  sts_total_22
  xi
are added for automatic output by MaterialOutputAction.
  Initializing
    Finished Initializing Equation Systems                                               [  0.54 s] [  778 MB]
    Finished Initializing Displaced Equation System                                      [  0.58 s] [  882 MB]
  Finished Initializing                                                                  [  1.18 s] [  875 MB]
Finished Setting Up                                                                      [  4.38 s] [  877 MB]
Framework Information:
MOOSE Version:           git commit 2bbd8e54e1 on 2024-02-13
LibMesh Version:         
PETSc Version:           3.20.3
SLEPc Version:           3.20.1
Current Time:            Wed Feb 21 18:19:19 2024
Executable Timestamp:    Wed Feb 21 17:20:28 2024

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
  Num DOFs:                29850709
  Num Local DOFs:          3731891
  Variables:               "B_old" { "xi_old" "I2_old" "mu_old" "lambda_old" "gamma_old" "alpha_in" "B_in" 
                             } { "alpha_in_dummy" "B_in_dummy" } { "alpha_grad_x" "alpha_grad_y" "alpha_grad_z" 
                             "B" "I1" ... "sts_total_12" "sts_total_20" "sts_total_21" "sts_total_22" 
                             "xi" } 
  Finite Element Types:    "LAGRANGE" "MONOMIAL" "MONOMIAL" "MONOMIAL" 
  Approximation Orders:    "FIRST" "CONSTANT" "FIRST" "CONSTANT" 

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
      Finished Computing Initial Material Values                                         [  6.02 s] [  730 MB]
    Finished Setting Up Materials                                                        [  6.02 s] [  730 MB]
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
sub_app0:   Num DOFs:                5185460
sub_app0:   Num Local DOFs:          648200
sub_app0:   Variables:               { "alpha_old" "B_old" } { "alpha_old_dummy" "B_old_dummy" } { "xi_old" "I2_old" 
sub_app0:                              "mu_old" "lambda_old" "gamma_old" "alpha_checked" "B_checked" } { "alpha_checked_dummy" 
sub_app0:                              "B_checked_dummy" } { "alpha_grad_x_sub" "alpha_grad_y_sub" "alpha_grad_z_sub" 
sub_app0:                              } 
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
    Finished Computing Initial Material Properties                                       [  6.47 s] [  659 MB]
  Finished Performing Initial Setup                                                      [ 53.71 s] [  714 MB]

Time Step 0, time = 0

Time Step 1, time = 5.9e-08, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 1, time = 5.9e-08, dt = 5.9e-08
sub_app0:  Solve Converged!
Still Executing
    Computing Residual                                                                   [  7.14 s] [  672 MB]
    Currently Computing Jacobian                                                         [  7.84 s] [  663 MB]
 Solve Converged!
  Finished Solving                                                                       [ 15.17 s] [  624 MB]
Still Executing
Time Step 2, time = 1.18e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 2, time = 1.18e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
.
    Computing Residual                                                                   [  8.78 s] [  630 MB]
 Solve Converged!
  Finished Solving                                                                       [  9.04 s] [  635 MB]
Still Executing
Time Step 3, time = 1.77e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 3, time = 1.77e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
.
    Computing Residual                                                                   [ 10.02 s] [  589 MB]
 Solve Converged!
  Finished Solving                                                                       [ 10.28 s] [  513 MB]
Still Executing.
Time Step 4, time = 2.36e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 4, time = 2.36e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
.
    Computing Residual                                                                   [  9.45 s] [  528 MB]
 Solve Converged!
  Finished Solving                                                                       [  9.68 s] [  453 MB]
Still Executing..
Time Step 5, time = 2.95e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 5, time = 2.95e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
..
    Computing Residual.                                                                  [ 12.02 s] [  433 MB]
 Solve Converged!
  Finished Solving                                                                       [ 12.22 s] [  358 MB]
Still Executing.
Finished Executing                                                                       [228.06 s] [  402 MB]