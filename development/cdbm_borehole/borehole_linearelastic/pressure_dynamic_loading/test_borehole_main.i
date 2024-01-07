##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4 [Check!]
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677

# under the stress state, 
# -135MPa 70MPa; 70MPa -120e6
# the principal stress is 198e6 

# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.5
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.4) * 198e6) = 233.6m

# Use mesh size = 25m, resolved by more than 7 elements [Check!]

# Use 1% overstress

# Use dt = 2e-4
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../../meshgenerator/cdbm/borehole_wofaults_small/borehole_wofaults.msh'
  []
  [./subdomain_id] 
    input = msh
    type = ElementSubdomainIDGenerator 
    element_ids = '
    21 265 499 955 995 1227 1543 2001 2397 2413 2747 2827 3071 3855 4137 4305 4315 4497 4521 4575 4693 4787 5059 5633 5867 6125 6219 6403 6487 6831 6929 7319 7365 7719 7833 7849 8029 8043 8071 8299 8305 8595 8596 9345 9375 9423 9603 9667 9939 9943 10231 10535 10715 10729 10905 11069 11221 11289 11401 11413 11523 11619 11853 12893 12967 13047 13063 13113 13199 13327 13353 13541 13643 13711 13721 13737 13769 13839 13939 13955 20 264 498 954 994 1226 1542 2000 2396 2412 2746 2826 3070 3854 4136 4304 4314 4496 4520 4574 4692 4786 5058 5632 5866 6124 6218 6402 6486 6830 6928 7318 7364 7718 7832 7848 8028 8042 8070 8298 8304 8594 8597 9344 9374 9422 9602 9666 9938 9942 10230 10534 10714 10728 10904 11068 11220 11288 11400 11412 11522 11618 11852 12892 12966 13046 13062 13112 13198 13326 13352 13540 13642 13710 13720 13736 13768 13838 13938 13954 23 267 501 953 993 1225 1541 2003 2395 2411 2745 2825 3073 3853 4135 4307 4317 4495 4519 4577 4691 4785 5057 5635 5869 6127 6221 6401 6489 6829 6927 7321 7363 7717 7831 7851 8027 8041 8073 8297 8303 8593 8599 9343 9377 9425 9601 9665 9937 9945 10229 10533 10713 10731 10903 11071 11223 11287 11403 11415 11525 11617 11855 12891 12969 13049 13065 13115 13201 13325 13355 13539 13641 13713 13719 13739 13767 13837 13941 13953 22 266 500 952 992 1224 1540 2002 2394 2410 2744 2824 3072 3852 4134 4306 4316 4494 4518 4576 4690 4784 5056 5634 5868 6126 6220 6400 6488 6828 6926 7320 7362 7716 7830 7850 8026 8040 8072 8296 8302 8592 8598 9342 9376 9424 9600 9664 9936 9944 10228 10532 10712 10730 10902 11070 11222 11286 11402 11414 11524 11616 11854 12890 12968 13048 13064 13114 13200 13324 13354 13538 13640 13712 13718 13738 13766 13836 13940 13952 812 813 1414 1415 1676 1677 1684 1685 1720 1721 1838 1839 2072 2073 2086 2087 2196 2197 2240 2241 2594 2595 2774 2775 3724 3725 3818 3819 3892 3893 4002 4003 4412 4413 4616 4617 4840 4841 4888 4889 4930 4931 5088 5089 5148 5149 5328 5329 5636 5637 5670 5671 5772 5773 5984 5985 6228 6229 6358 6359 6392 6393 6560 6561 6590 6591 6738 6739 6770 6771 7098 7099 7130 7131 7434 7435 7548 7549 8010 8011 8022 8023 8162 8163 8406 8407 8412 8413 8538 8539 8576 8577 8938 8939 9036 9037 9354 9355 9468 9469 9694 9695 10090 10091 10686 10687 10734 10735 11074 11075 11580 11581 11606 11607 12616 12617 13120 13121 13280 13281 13302 13303 13422 13423 13486 13487 13840 13841 14124 14125 14258 14259 14276 14277 14370 14371 14388 14389 15070 15071 
    '
    subdomain_ids = '
    1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 
    '
  []
  [./sidesets_leftbc]
    input = subdomain_id
    type = SideSetsAroundSubdomainGenerator
    block = 1000
    new_boundary = 'left'
    external_only = true
  []
  [./sidesets_rightbc]
    input = sidesets_leftbc
    type = SideSetsAroundSubdomainGenerator
    block = 1001
    new_boundary = 'right'
    external_only = true
  []
  [./sidesets_topbc]
    input = sidesets_rightbc
    type = SideSetsAroundSubdomainGenerator
    block = 1002
    new_boundary = 'top'
    external_only = true
  []
  [./sidesets_bottombc]
    input = sidesets_topbc
    type = SideSetsAroundSubdomainGenerator
    block = 1003
    new_boundary = 'bottom'
    external_only = true
  []
  [./sidesets_boreholebc]
    input = sidesets_bottombc
    type = SideSetsAroundSubdomainGenerator
    block = 1004
    new_boundary = 'borehole'
    external_only = true
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0  -1000  0'
    new_boundary = corner_ptr1
    input = sidesets_boreholebc
  []
  [./extranodeset2]
      type = ExtraNodesetGenerator
      coord = '1000  0  0'
      new_boundary = corner_ptr2
      input = extranodeset1
  []
  [./extranodeset3]
      type = ExtraNodesetGenerator
      coord = '0   1000  0'
      new_boundary = corner_ptr3
      input = extranodeset2
  []
  [./extranodeset4]
      type = ExtraNodesetGenerator
      coord = '-1000  0  0'
      new_boundary = corner_ptr4
      input = extranodeset3
  []
  use_displaced_mesh = true
  displacements = 'disp_x disp_y'
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  #damping ratio
  q = 0.2
  
[]

[AuxVariables]
  [./arctanyx]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
      [../]
    [../]
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_y
  []
  [./Reactionx]
    type = StiffPropDamping
    variable = 'disp_x'
    component = '0'
  []
  [./Reactiony]
    type = StiffPropDamping
    variable = 'disp_y'
    component = '1'
  []
[]

[Materials]
  [stress_medium]
      type = ComputeLinearElasticStress
  []
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2670
  []
[]

[AuxKernels]
  [getarctanyx]
    type = CompArctanyx
    variable = arctanyx
  []
[]

[Functions]
  [func_pressure_x]
    type = InitialStressXYPressureBorehole_x_fast
  []
  [func_pressure_y]
    type = InitialStressXYPressureBorehole_y_fast
  []
[]

[Executioner]
  type = Transient
  dt = 2e-4
  end_time = 40.0
  # num_steps = 5
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 200
  [sample_snapshots]
    type = Exodus
    interval = 2000
  []
  [snapshots]
    type = Exodus
    interval = 2000
    overwrite = true
  []
  [checkpoints]
    type = Checkpoint
    interval = 2000
    num_files = 2
  []
[]

[BCs]
  [Pressure_borehole_x]
    type = Pressure
    boundary = borehole
    function = func_pressure_x
    variable = disp_x
  []
  [Pressure_borehole_y]
    type = Pressure
    boundary = borehole
    function = func_pressure_y
    variable = disp_y
  []
  [./dashpot_top_x]
    type = NonReflectDashpotBC
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    p_wave_speed = 6000
    shear_wave_speed = 3464
    boundary = top
  []
  [./dashpot_top_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = top
  []
  [./dashpot_bottom_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_bottom_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_left_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_left_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_right_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
  [./dashpot_right_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
[]