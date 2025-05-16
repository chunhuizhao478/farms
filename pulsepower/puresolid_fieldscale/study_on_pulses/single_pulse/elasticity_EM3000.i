E = 3.2e10
nu = 0.2
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

Gc_const = 57
l = 4e-3 # N * h, N: number of elements, h: element size

[Adaptivity]
  max_h_level = 3
  marker = 'combo'
  cycles_per_step = 2
  [Markers]
      [./combo]
          type = ComboMarker
          markers = 'damage_marker strain_energy_marker'
      [../]
      [damage_marker]
        type = ValueThresholdMarker
        variable = d
        refine = 0.5
      []
      [strain_energy_marker]
        type = ValueThresholdMarker
        variable = psie_active
        refine = '${fparse 1.0*0.5*Gc_const/l}'
      []      
  []
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc_const=${Gc_const};l=${l}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = psie_active
    source_variable = psie_active
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshfile/fieldscale_small.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [fy]
  []
  [d]
    family = LAGRANGE
    order = FIRST
  []
  #err measurement of active strain energy
  [eng_err]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [vel_z]
    family = LAGRANGE
    order = FIRST
  []
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  [accel_z]
    family = LAGRANGE
    order = FIRST
  []
  #
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [./error_measure]
    type = ErrorPsiMeasure
    variable = eng_err
  [../]
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.25
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.5
    execute_on = 'TIMESTEP_END'
  []
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.25
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.5
    execute_on = 'TIMESTEP_END'
  []
  [accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = 0.25
    execute_on = 'TIMESTEP_END'
  []
  [vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = 0.5
    execute_on = 'TIMESTEP_END'
  []
  #get pulse load aux
  [get_pulse_load_aux]
    type = FunctionAux 
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
[]

[Physics/SolidMechanics/Dynamic]
  [all]
    add_variables = true
    hht_alpha = 0.11
    newmark_beta = 0.25
    newmark_gamma = 0.5
    use_automatic_differentiation = true
    # mass_damping_coefficient = 0.1
    # stiffness_damping_coefficient = 0.1
    density = 2450
  []
[]

[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 3.630e5
    shape_param_beta = 6.677e5
    rise_time = 2e-6
    single_pulse_duration = 2e-5
    EM = 3000
    gap = 0.05
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    discharge_center = '0.5 0.5 0.025'
    number_of_pulses = 1
  []
[]

[BCs]
  #confinement
  [./Pressure]
    #assign pressure on inner surface
    [pressure_inner]
      boundary = Cylinder
      function = func_tri_pulse
      displacements = 'disp_x disp_y disp_z'
      use_displaced_mesh = true
    []             
  []    
  # fix all displacements in top
  [./fix_top_z]
    type = DirichletBC
    variable = disp_z
    boundary = Top
    value = 0
  []
  #fix all displacements in bottom
  [./fix_bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = Bottom
    value = 0
  []
  # fix ptr
  [./fix_cptr1_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr
    value = 0
  []
  [./fix_cptr2_y]
    type = DirichletBC
    variable = disp_y
    boundary = corner_ptr
    value = 0
  []
  [./fix_cptr3_z]
    type = DirichletBC
    variable = disp_z
    boundary = corner_ptr
    value = 0
  []
[]

[BCs]
  #add dampers
  [damp_left_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = Left
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_left_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = Left
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_left_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = Left
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  #
  [damp_right_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = Right
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_right_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = Right
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_right_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = Right
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  #
  [damp_front_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = Front
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_front_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = Front
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_front_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = Front
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  [] 
  # 
  [damp_back_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = Back
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_back_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = Back
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  []
  [damp_back_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = Back
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = 2220
    p_wave_speed = 3630
    density = 2450
  [] 
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  # [density]
  #   type = ADGenericConstantMaterial
  #   prop_names = 'density'
  #   prop_values = '2450'
  # []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 0.5e-8
  end_time = 2e-5

  fixed_point_max_its = 5
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10

  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
[]

[Outputs]
  exodus = true
  time_step_interval = 100
  print_linear_residuals = false
  csv = true
[]

