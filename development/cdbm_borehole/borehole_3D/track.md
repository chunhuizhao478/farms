(moose) andyz@wirelessprv-10-193-217-157 farms % bash development/cdbm_borehole/borehole_3D/splitmesh.sh        


*** Warning, This code is deprecated and will be removed in future versions:
Please update your main.C to adapt new main function in MOOSE framework, see'test/src/main.C in MOOSE as an example of moose::main()'. 



*** Info ***
TensorMechanics Action: selecting 'total small strain' formulation. Use `incremental = true` to select 'incremental small strain' instead.
Splitting 8 ways...
    - writing 1 files per process...
Mesh meta data written into "/Users/andyz/projects/farms/development/cdbm_borehole/borehole_3D/foo_main.cpr/meta_data_mesh.rd".
Finished Setting Up                                                                      [  3.99 s] [  326 MB]
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There are 3 unused database options. They are:
Option left: name:--split-file value: development/cdbm_borehole/borehole_3D/foo_main.cpr source: command line
Option left: name:--split-mesh value: 8 source: command line
Option left: name:-i value: development/cdbm_borehole/borehole_3D/test_borehole_main.i source: command line


*** Warning, This code is deprecated and will be removed in future versions:
Please update your main.C to adapt new main function in MOOSE framework, see'test/src/main.C in MOOSE as an example of moose::main()'. 



*** Info ***
TensorMechanics Action: selecting 'total small strain' formulation. Use `incremental = true` to select 'incremental small strain' instead.
sub_app0: Initializingsub_app0: 
sub_app0:   Finished Initializing Equation Systems                                       [  0.39 s] [  241 MB]
sub_app0: Finished Initializing                                                          [  0.41 s] [  241 MB]
Currently Setting Up
  Finished Instantiating Sub-Apps                                                        [  1.38 s] [  241 MB]
The following total 37 aux variables:
  I1
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
  sts_total_00
  sts_total_01
  sts_total_02
  sts_total_10
  sts_total_11
  sts_total_12
  sts_total_20
  sts_total_21
  sts_total_22
are added for automatic output by MaterialOutputAction.
  Initializing
    Finished Initializing Equation Systems                                               [  2.55 s] [  831 MB]
  Finished Initializing                                                                  [  6.56 s] [  772 MB]
Finished Setting Up                                                                      [  9.79 s] [  772 MB]
Framework Information:
MOOSE Version:           git commit 2bbd8e54e1 on 2024-02-13
LibMesh Version:         
PETSc Version:           3.20.3
SLEPc Version:           3.20.1
Current Time:            Wed Feb 21 18:25:00 2024
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
      Finished Computing Initial Material Values                                         [  7.05 s] [  641 MB]
    Finished Setting Up Materials                                                        [  7.07 s] [  641 MB]
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
    Finished Computing Initial Material Properties                                       [  8.05 s] [  600 MB]
  Finished Performing Initial Setup                                                      [ 56.91 s] [  607 MB]

Time Step 0, time = 0

Time Step 1, time = 5.9e-08, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 1, time = 5.9e-08, dt = 5.9e-08
sub_app0:  Solve Converged!
Still Executing
    Computing Residual                                                                   [  7.85 s] [  623 MB]
    Currently Computing Jacobian                                                         [  8.32 s] [  647 MB]
 Solve Converged!
  Finished Solving                                                                       [ 16.32 s] [  644 MB]
Still Executing.
Time Step 2, time = 1.18e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 2, time = 1.18e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
..
    Computing Residual.                                                                  [ 14.59 s] [  576 MB]
    Currently Computing Jacobian.                                                        [ 14.33 s] [  564 MB]
 Solve Converged!
  Finished Solving                                                                       [ 29.02 s] [  581 MB]
Still Executing.
Time Step 3, time = 1.77e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 3, time = 1.77e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
....
    Computing Residual...                                                                [ 24.64 s] [  586 MB]
    Currently Computing Jacobian...                                                      [ 21.97 s] [  594 MB]
 Solve Converged!
  Finished Solving                                                                       [ 46.74 s] [  601 MB]
Still Executing..
Time Step 4, time = 2.36e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 4, time = 2.36e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
....
    Computing Residual...                                                                [ 21.28 s] [  569 MB]
    Currently Computing Jacobian...                                                      [ 22.86 s] [  584 MB]
 Solve Converged!
  Finished Solving                                                                       [ 44.22 s] [  611 MB]
Still Executing...
Time Step 5, time = 2.95e-07, dt = 5.9e-08
sub_app0: 
sub_app0: Time Step 5, time = 2.95e-07, dt = 5.9e-08
sub_app0:  Solve Converged!
....
    Computing Residual..                                                                 [ 20.01 s] [  581 MB]
    Currently Computing Jacobian...                                                      [ 24.87 s] [  571 MB]
 Solve Converged!
  Finished Solving                                                                       [ 44.95 s] [  582 MB]
Still Executing..
Finished Executing                                                                       [422.19 s] [  454 MB]