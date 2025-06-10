E = 50e9
nu = 0.373
# ft = 25.5e6
Gc_const = 57
density = 2600
dx_min = 2.5e-5

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  1e-4 
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#fieldscale small: dx = 1e-3 < l = 1.64e-3, 3x adaptivity levels

[Adaptivity]
  max_h_level = 5
  marker = 'combo'
  cycles_per_step = 1
  [Markers]
      [./combo]
        type = FarmsComboMarker
        markers = 'damage_marker strain_energy_marker'
        meshsize_marker = 'meshsize_marker'
        block = 5
      [../]
      [damage_marker]
        type = ValueThresholdMarker
        variable = d
        refine = 0.5
        block = 5
      []
      [strain_energy_marker]
        type = ValueThresholdMarker
        variable = psie_active
        refine = '${fparse 1.0*3/8*Gc_const/l}'
        block = 5
      []   
      # if mesh_size > dxmin, refine
      # if mesh_size < dxmin/100, coarsen (which never happens)
      # otherwise, do nothing
      [meshsize_marker]
        type = ValueThresholdMarker
        variable = mesh_size
        refine = '${dx_min}'
        coarsen = '${fparse dx_min/100}'
        third_state = DO_NOTHING
        block = 5
      [] 
  []
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc_const=${Gc_const};l=${l};dx_min=${dx_min}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
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
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active mesh_size'
  []
[]

[GlobalParams]
  # displacements = 'disp_x disp_y disp_z'
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    # file =  '../../../meshfile/fieldscale_test1.msh'
    file =  '../../../2dmeshfile/fieldscale_test2_2d.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  # displacements = 'disp_x disp_y disp_z'
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  # [disp_z]
  # []
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
  # [vel_z]
  #   family = LAGRANGE
  #   order = FIRST
  # []
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  # [accel_z]
  #   family = LAGRANGE
  #   order = FIRST
  # []
  #
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  #
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
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
  # [accel_z]
  #   type = NewmarkAccelAux
  #   variable = accel_z
  #   displacement = disp_z
  #   velocity = vel_z
  #   beta = 0.25
  #   execute_on = 'TIMESTEP_END'
  # []
  # [vel_z]
  #   type = NewmarkVelAux
  #   variable = vel_z
  #   acceleration = accel_z
  #   gamma = 0.5
  #   execute_on = 'TIMESTEP_END'
  # []
  #get pulse load aux
  [get_pulse_load_aux]
    type = FunctionAux 
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  #mesh size aux
  [./max]
    type = ElementLengthAux
    variable = mesh_size
    method = max
    execute_on = TIMESTEP_BEGIN
  [../]
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
    density = ${density}
  []
[]

[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 4.658e5
    shape_param_beta = 4.661e5
    rise_time = 3e-6
    single_pulse_duration = 4e-5
    EM = 0.03
    gap = 0.02
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    discharge_center = '0 0 0'
    number_of_pulses = 10
    peak_pressure = 150e6
  []
[]

[BCs]
  #confinement
  [./Pressure]
    #assign pressure on inner surface
    [pressure_inner]
      boundary = 3
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = true
    []             
  []
  # fix all displacements in top
  # [./fix_top_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = 3
  #   value = 0
  # []
  # #fix all displacements in bottom
  # [./fix_bottom_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = 6
  #   value = 0
  # []    
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
  # [./fix_cptr3_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = corner_ptr
  #   value = 0
  # []
[]

[BCs]
  #add dampers
  [damp_outer_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    # displacements = 'disp_x disp_y disp_z'
    # velocities = 'vel_x vel_y vel_z'
    # accelerations = 'accel_x accel_y accel_z'
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = 1
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
  []
  [damp_outer_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    # displacements = 'disp_x disp_y disp_z'
    # velocities = 'vel_x vel_y vel_z'
    # accelerations = 'accel_x accel_y accel_z'
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = 1
    beta = 0.25
    gamma = 0.5
    alpha = 0.11
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
  []
  # [damp_outer_z]
  #   type = FarmsNonReflectDashpotBC
  #   variable = disp_z
  #   displacements = 'disp_x disp_y disp_z'
  #   velocities = 'vel_x vel_y vel_z'
  #   accelerations = 'accel_x accel_y accel_z'
  #   component = 2
  #   boundary = 4
  #   beta = 0.25
  #   gamma = 0.5
  #   alpha = 0.11
  #   shear_wave_speed = ${Cs}
  #   p_wave_speed = ${Cp}
  #   density = ${density}
  # [] 
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
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  automatic_scaling = true

  l_max_its = 100
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15

  # dt = 0.5e-7
  end_time = 20e-5

  fixed_point_max_its = 5
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 1e-8
    cutback_factor_at_failure = 0.5
    optimal_iterations = 5
    growth_factor = 1.25
    max_time_step_bound = 1e-7
  []
  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
[]

[Outputs]
  exodus = true
  time_step_interval = 10
  print_linear_residuals = false
  csv = true
  [checkpoint]
      type = Checkpoint
      time_step_interval = 20
      num_files = 2
  []
[]