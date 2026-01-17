# Rate-Dependent Phase Field Fracture - Elastodynamics Input
# Based on Hofacker & Miehe (2012) - IJNME 93:276-301
# This file handles the elastodynamics (displacement field)

E = 50e9
nu = 0.3
Gc_const = 100
density = 2600

# Rate-dependent viscosity parameter [Ns/mm^2]
# Set eta = 0 for rate-independent behavior
# Typical values: 1e-12 (near rate-independent) to 1e-6 (strong rate-dependent)
eta_viscosity = 0

##parametric study on confinement pressure##
confinement_pressure = 10e6

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l = 2e-4
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#finite element properties
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
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
    reduction_type = 'sum'
  []
  [pp_transfer_dissipated_energy_first_step]
    type = MultiAppPostprocessorTransfer
    from_multi_app = 'fracture'
    from_postprocessor = 'dissipated_energy_first_step'
    to_postprocessor = 'dissipated_energy_first_step'
    reduction_type = 'sum'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

#initial damage box
bottom_left1 = '0 -4e-4 0'
top_right1 = '0.0925 4e-4 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../../pf_code2d_puresolid_pulseloadmaxmin/2dmeshfile/fieldscale_test1_2d_amr_fieldscale_debug_enlargedrefine.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '1.0 1.0 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node = true
  []
  [./subdomain_id]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left1}
    top_right = ${top_right1}
    location = INSIDE
    block_id = 1
    input = extranodeset1
  []
  displacements = 'disp_x disp_y'
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
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
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
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  [dissipated_energy_density]
    family = MONOMIAL
    order = CONSTANT
  []
  [fx]
  []
  [fy]
  []
  [fz]
  []
  [fconfinementx]
  []
  [fconfinementy]
  []
  [fconfinementz]
  []
  [fdampx]
  []
  [fdampy]
  []
  [fdampz]
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
  [get_pulse_load_aux]
    type = FunctionAux
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  [mesh_size_aux]
    type = MeshSize
    variable = mesh_size
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
  []
[]

[Kernels]
  [solid_x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y'
  []
  [solid_y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y'
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
[]

[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 4.658e5
    shape_param_beta = 4.661e5
    rise_time = 3e-6
    single_pulse_duration = 2e-5
    EM = 0.25
    gap = 0.008
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    fitting_param_exponent = 0.25
    discharge_center = '0 0 0'
    number_of_pulses = 100
    base_factor = 8000
    peak_pressure = 105e6
  []
[]

[BCs]
  [./Pressure]
    [pressure_inner]
      boundary = 3
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fx
      save_in_disp_y = fy
    []
    [static_pressure_outer]
      boundary = 1
      factor = ${confinement_pressure}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx
      save_in_disp_y = fconfinementy
    []
  []
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
  [damp_outer_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = 1
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
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = 1
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampy
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
    output_properties = 'elastic_strain psie_active psie'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '

  line_search = 'basic'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 30

  end_time = 100e-5

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 5e-8
    iteration_window = 0
    cutback_factor_at_failure = 0.5
    optimal_iterations = 20
    growth_factor = 1.25
    max_time_step_bound = 5e-8
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
    time_step_interval = 20
    show = 'd vel_x vel_y vel_z stress_00 stress_11 stress_01'
  [../]
  [checkpoint]
    type = Checkpoint
    time_step_interval = 100
    num_files = 2
  []
  [csv]
    type = CSV
    execute_on = 'initial timestep_end'
    time_step_interval = 1
    show = 'full_energy solid_elastic_energy_total solid_kinetic_energy_total dissipated_energy_total full_input_energy damping_work'
  []
[]

###############################Energy Calculation##############################

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

[Postprocessors]
  [external_work]
    type = FarmsExternalWork
    boundary = '3'
    forces = 'fx fy fz'
  []
  [confinement_work]
    type = FarmsExternalWork
    boundary = '1'
    forces = 'fconfinementx fconfinementy fconfinementz'
  []
  [damping_work]
    type = FarmsExternalWork
    boundary = '1'
    forces = 'fdampx fdampy fdampz'
  []
[]

[Postprocessors]
  [full_input_energy]
    type = ParsedPostprocessor
    expression = '-1 * external_work - confinement_work - damping_work'
    pp_names = 'external_work confinement_work damping_work'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

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

[Postprocessors]
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + dissipated_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
