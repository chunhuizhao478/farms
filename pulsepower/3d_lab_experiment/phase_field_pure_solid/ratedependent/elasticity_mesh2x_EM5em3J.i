E = 50e9
nu = 0.373
# ft = 25.5e6
Gc_const = 40
density = 2600

# Rate-dependent phase field viscosity (Hofacker & Miehe 2012)
# Typical range: 1e-6 to 1e-3 Pa*s for dynamic fracture
# Higher eta = slower damage evolution, prevents spurious erosion
eta_viscosity = 1e-5

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  2e-4 #l = 5e-4, ft = 61 MPa
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
# HHT-alpha damping: 0 = no damping, -0.05 to -0.1 = conservative damping
# Enable HHT damping to reduce high-frequency numerical oscillations
# that may cause spurious damage erosion from reflected waves
hht_alpha = 0 # Try -0.05 if damage erosion persists with rate-dependent formulation
#----------------------------------------------------#

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture_mesh2x_EM5em3J.i
    cli_args = 'Gc_const=${Gc_const};l=${l};eta_viscosity=${eta_viscosity}'
    execute_on = 'INITIAL TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = 'd dissipated_energy_density'
    source_variable = 'd dissipated_energy_density'
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active mesh_size'
  []
  [pp_transfer_dissipated_energy_total]
    type = MultiAppPostprocessorTransfer
    from_multi_app = 'fracture'
    from_postprocessor = 'dissipated_energy_dynamic'
    to_postprocessor = 'dissipated_energy_dynamic'
    reduction_type = 'sum' #this should not have effect in a single app
  []
  [pp_transfer_dissipated_energy_first_step]
    type = MultiAppPostprocessorTransfer
    from_multi_app = 'fracture'
    from_postprocessor = 'dissipated_energy_first_step'
    to_postprocessor = 'dissipated_energy_first_step'
    reduction_type = 'sum' #this should not have effect in a single app
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

#initial damage box 1
bottom_left1 = '-0.002 -4e-4 0'
top_right1 = '0.002 4e-4 0.06'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../3dmeshfile/cylinder_sample_refined_cross.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  [./subdomain_id]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left1}
    top_right = ${top_right1}
    location = INSIDE
    block_id = 2
    input = extranodeset1
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = FIRST
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
  []
  [disp_z]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxVariables]
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
  #
  [vel_z]
    family = LAGRANGE
    order = FIRST
  []
  #
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
  #
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  #
  [dissipated_energy_density]
    family = MONOMIAL
    order = CONSTANT
  []
  #reaction force
  [fx]
  []
  [fy]
  []
  [fz]
  []
  [fdampx]
  []
  [fdampy]
  []
  [fdampz]
  []
  [fdampx_top]
  []
  [fdampy_top]
  []
  [fdampz_top]
  []
  [fdampx_bottom]
  []
  [fdampy_bottom]
  []
  [fdampz_bottom]
  []
[]

[AuxKernels]
  [./error_measure]
    type = ErrorPsiMeasure
    variable = eng_err
  [../]
  #
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  #
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  #
  [accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  #get pulse load aux
  [get_pulse_load_aux]
    type = FunctionAux
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  #mesh size aux
  [mesh_size_aux]
    type = MeshSize
    variable = mesh_size
    #execute_on = 'TIMESTEP_END'
    execute_on = 'INITIAL'
  []
[]

[Kernels]
  [solid_x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y disp_z'
  []
  [solid_y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y disp_z'
  []
  [solid_z]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y disp_z'
  []
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_x
    acceleration = accel_x
    density = ${density}
    alpha = ${hht_alpha}
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_y
    acceleration = accel_y
    density = ${density}
    alpha = ${hht_alpha}
  []
  [inertia_z]
    type = ADInertialForce
    variable = disp_z
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_z
    acceleration = accel_z
    density = ${density}
    alpha = ${hht_alpha}
  []
[]

[Functions]
  # [func_tri_pulse]
  #   type = ElkPulseLoadExperimentWu2022Paper #need to adopt function
  #   shape_param_alpha = 4.658e5
  #   shape_param_beta = 4.661e5
  #   rise_time = 3e-6
  #   single_pulse_duration = 1e-5
  #   Pmax_coefficients = '-8.48306 105.794 -451.486 696.525 -146.249'
  #   discharge_center = '0 0 0.03'
  #   number_of_pulses = 10
  #   r_max_mm = 4.5
  # []
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 4.658e5
    shape_param_beta = 4.661e5
    rise_time = 3e-6
    single_pulse_duration = 1e-5
    EM = 0.005
    gap = 0.008
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    fitting_param_exponent = 0.25
    discharge_center = '0 0 0.03'
    number_of_pulses = 10
    base_factor = 8000
    # peak_pressure = 200e6 #if peak pressure is specified, the depth variation is ignored
  []
[]

[BCs]
  #confinement
  [./Pressure]
    #assign pressure on inner surface
    [pressure_inner]
      boundary = 3 #confirm the boundary id
      function = func_tri_pulse
      displacements = 'disp_x disp_y disp_z'
      use_displaced_mesh = false
      save_in_disp_x = fx
      save_in_disp_y = fy
      save_in_disp_z = fz
    []
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
  [./fix_cptr2_z]
    type = DirichletBC
    variable = disp_z
    boundary = corner_ptr
    value = 0
  []
  #add dampers
  #outer #4
  [damp_outer_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = 4
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampx
  []
  [damp_outer_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = 4
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampy
  []
  [damp_outer_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = 4
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampz
  []
  #top boundary #5
  [damp_top_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = 5
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampx_top
  []
  [damp_top_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = 5
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampy_top
  []
  [damp_top_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = 5
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampz_top
  []
  #bottom #2
  [damp_bottom_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = 2
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampx_bottom
  []
  [damp_bottom_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = 2
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampy_bottom
  []
  [damp_bottom_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = 2
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampz_bottom
  []
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G density'
    prop_values = '${K} ${G} ${density}'
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
    #output_properties = 'elastic_strain psie_active psie'
    output_properties = 'psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    #output_properties = 'stress'
    #outputs = exodus
  []
  [strain]
    type = ADComputeSmallStrain
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

  #scalable to large problems
  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
  # petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -snes_type -ksp_gmres_restart'
  petsc_options_value = 'gmres hypre boomeramg 0.7 4 5 0.3 vinewtonrsls 100'

  # automatic_scaling = true
  line_search = 'basic'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  # Add more iterations before failure
  nl_max_its = 30

  # dt = 0.5e-7
  end_time = 10e-5

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 1e-8
    iteration_window = 0 #the adaptive time stepping happens at number of iterations <-> 'optimal_iterations plus/minus iteration_window'
    cutback_factor_at_failure = 0.5
    optimal_iterations = 20
    growth_factor = 1.25
    max_time_step_bound = 1e-7
  []
  [./TimeIntegrator]
    type = NewmarkBeta
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    time_step_interval = 40
    show = 'd pulse_load_aux vel_x vel_y vel_z'
  [../]
  [checkpoint]
      type = Checkpoint
      time_step_interval = 100
      num_files = 2
  []
  [csv]
    type = CSV
    execute_on = 'initial timestep_end'
    time_step_interval = 40
    show = 'full_energy solid_elastic_energy_total solid_kinetic_energy_total dissipated_energy_total full_input_energy damping_work'
  []
[]

###############################Energy Calculation##############################

#fracture energy
###############################################################################
# [Postprocessors]
#   [dissipated_energy_dynamic]
#       type = ElementIntegralVariablePostprocessor
#       variable = dissipated_energy_density
#       execute_on = 'INITIAL TIMESTEP_END'
#   []
# []

# first step dissipated energy
# this postprocessor computes the dissipated energy during the first time step
# and substract in the full energy calculation
# [Postprocessors]
#   [dissipated_energy_first_step]
#     type = FirstStepElementIntegralVariablePostprocessor
#     variable = dissipated_energy_density
#     execute_on = 'TIMESTEP_END'
#   []
# []

#receive the data from subapp
[Postprocessors]
  [./dissipated_energy_dynamic]
    type = Receiver
  [../]
  [./dissipated_energy_first_step]
    type = Receiver
  [../]
[]

[Postprocessors]
  [dissipated_energy_total]
      type = ParsedPostprocessor
      expression = 'dissipated_energy_dynamic - dissipated_energy_first_step'
      pp_names = 'dissipated_energy_dynamic dissipated_energy_first_step'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# input energy
###############################################################################
[Postprocessors]
  [external_work]
    type = FarmsExternalWork
    boundary = '4'
    forces = 'fx fy fz'
  []
  [damping_work]
    type = FarmsExternalWork
    boundary = '3'
    forces = 'fdampx fdampy fdampz'
  []
  [damping_work_top]
    type = FarmsExternalWork
    boundary = '5'
    forces = 'fdampx_top fdampy_top fdampz_top'
  []
  [damping_work_bottom]
    type = FarmsExternalWork
    boundary = '2'
    forces = 'fdampx_bottom fdampy_bottom fdampz_bottom'
  []
[]

[Postprocessors]
  [full_input_energy]
      type = ParsedPostprocessor
      expression = '-1 * external_work - damping_work - damping_work_top - damping_work_bottom'
      pp_names = 'external_work damping_work damping_work_top damping_work_bottom'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# solid kinetic energy
###############################################################################
[AuxVariables]
  [solid_kinetic_energy]
      order = CONSTANT
      family = MONOMIAL
  []
[]

[AuxKernels]
  [solid_kinetic_energy]
      type = ADKineticEnergyAux
      variable = solid_kinetic_energy
      newmark_velocity_x = vel_x
      newmark_velocity_y = vel_y
      newmark_velocity_z = vel_z
      density = density
  []
[]

[Postprocessors]
  [solid_kinetic_energy_total]
      type = ElementIntegralVariablePostprocessor
      variable = solid_kinetic_energy
  []
[]

# solid elastic energy
###############################################################################
[Postprocessors]
  [solid_elastic_energy_dynamic]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie
  []
[]

[Postprocessors]
  [solid_elastic_energy_total]
      type = ParsedPostprocessor
      expression = 'solid_elastic_energy_dynamic'
      pp_names = 'solid_elastic_energy_dynamic'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# Full Energy
###############################################################################
[Postprocessors]
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + dissipated_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
