

*** Info ***
TensorMechanics Action: selecting 'total small strain' formulation. Use `incremental = true` to select 'incremental small strain' instead.
Currently Setting Up....sub_app0: Initializingsub_app0: 
sub_app0:   Initializing Equation Systems.....                                           [ 38.98 s] [  515 MB]
sub_app0: Finished Initializing                                                          [ 42.17 s] [  550 MB]
  Finished Instantiating Sub-Apps                                                        [ 79.26 s] [  552 MB]
The following total 28 aux variables:
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
are added for automatic output by MaterialOutputAction.
  Initializing
    Updating Because Mesh Changed
      Finished Updating Mesh                                                             [  3.25 s] [  492 MB]
    Finished Updating Because Mesh Changed                                               [  3.50 s] [  519 MB]
    Initializing Equation Systems..............                                          [ 84.48 s] [  587 MB]
    Initializing Displaced Equation System...............                                [ 86.00 s] [  451 MB]
  Finished Initializing                                                                  [182.24 s] [  521 MB]
Finished Setting Up                                                                      [307.43 s] [  497 MB]
Framework Information:
MOOSE Version:           git commit 2810766e21 on 2023-11-15
LibMesh Version:         
PETSc Version:           3.20.0
SLEPc Version:           3.19.2
Current Time:            Mon Nov 20 19:39:39 2023
Executable Timestamp:    Mon Nov 20 19:33:29 2023

Parallelism:
  Num Processors:          8
  Num Threads:             1

Mesh: 
  Parallel Type:           replicated
  Mesh Dimension:          2
  Spatial Dimension:       2
  Nodes:                   
    Total:                 849514
    Local:                 106549
    Min/Max/Avg:           105463/106549/106189
  Elems:                   
    Total:                 1694224
    Local:                 211778
    Min/Max/Avg:           211778/211778/211778
  Num Subdomains:          2
  Num Partitions:          8
  Partitioner:             metis

Nonlinear System:
  Num DOFs:                1699028
  Num Local DOFs:          213098
  Variables:               { "disp_x" "disp_y" } 
  Finite Element Types:    "LAGRANGE" 
  Approximation Orders:    "FIRST" 

Auxiliary System:
  Num DOFs:                116085570
  Num Local DOFs:          14515373
  Variables:               { "resid_x" "resid_y" "resid_slipweakening_x" "resid_slipweakening_y" "disp_slipweakening_x" 
                             "disp_slipweakening_y" "vel_slipweakening_x" "vel_slipweakening_y" "accel_slipweakening_x" 
                             "accel_slipweakening_y" } { "mu_s" "mu_d" } "ini_shear_stress" { "tangent_jump_rate" 
                             "tria_area_aux" } { "nodal_area" "B_old" } { "xi_old" "I2_old" "mu_old" 
                             "lambda_old" "gamma_old" "alpha_in" "B_in" } { "alpha_in_dummy" "B_in_dummy" 
                             } { "alpha_grad_x" "alpha_grad_y" "mechanical_strain_rate" "track_Cd" "eqv_plastic_strain" 
                             ... "eps_total_11" "eps_total_12" "eps_total_20" "eps_total_21" "eps_total_22" 
                             } 
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
  Performing Initial Setup
    Finished Setting Up Outputs                                                          [  0.93 s] [  464 MB]
    Finished Projecting Initial Solutions                                                [  6.62 s] [  471 MB]
    Currently Setting Up Materials
      Computing Initial Material Values........                                          [ 50.17 s] [  404 MB]
    Finished Setting Up Materials                                                        [ 50.21 s] [  404 MB]
sub_app0: Parallelism:
sub_app0:   Num Processors:          8
sub_app0:   Num Threads:             1
sub_app0: 
sub_app0: Mesh: 
sub_app0:   Parallel Type:           replicated
sub_app0:   Mesh Dimension:          2
sub_app0:   Spatial Dimension:       2
sub_app0:   Nodes:                   
sub_app0:     Total:                 849514
sub_app0:     Local:                 106549
sub_app0:     Min/Max/Avg:           105463/106549/106189
sub_app0:   Elems:                   
sub_app0:     Total:                 1694224
sub_app0:     Local:                 211778
sub_app0:     Min/Max/Avg:           211778/211778/211778
sub_app0:   Num Subdomains:          2
sub_app0:   Num Partitions:          8
sub_app0:   Partitioner:             metis
sub_app0: 
sub_app0: Nonlinear System:
sub_app0:   Num DOFs:                13553792
sub_app0:   Num Local DOFs:          1694224
sub_app0:   Variables:               { "alpha_sub" "B_sub" } { "alpha_sub_dummy" "B_sub_dummy" } 
sub_app0:   Finite Element Types:    "MONOMIAL" "MONOMIAL" 
sub_app0:   Approximation Orders:    "CONSTANT" "FIRST" 
sub_app0: 
sub_app0: Auxiliary System:
sub_app0:   Num DOFs:                42355600
sub_app0:   Num Local DOFs:          5294450
sub_app0:   Variables:               { "alpha_old" "B_old" } { "alpha_old_dummy" "B_old_dummy" } { "xi_old" "I2_old" 
sub_app0:                              "mu_old" "lambda_old" "gamma_old" "alpha_checked" "B_checked" } { "alpha_checked_dummy" 
sub_app0:                              "B_checked_dummy" } { "alpha_grad_x_sub" "alpha_grad_y_sub" "mechanical_strain_rate_sub" 
sub_app0:                              "track_Cd" } 
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
sub_app0: Currently Performing Initial Setupsub_app0: 
sub_app0:   Finished Building SemiLocalElemMap                                           [  7.85 s] [  398 MB]
sub_app0: Finished Performing Initial Setup                                              [ 15.71 s] [  400 MB]
sub_app0: 