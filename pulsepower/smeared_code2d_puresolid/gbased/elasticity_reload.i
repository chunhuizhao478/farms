E = 50e9
nu = 0.373
ft = 150e6 ##computed from pf
# Gc_const = 100
density = 2600
# dx_min = 5e-5

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  2e-4 
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#gradient activity parameters
kappa_i = ${fparse ft / E}
c0 = 1e-12 #minimum value of the gradient activity parameter for the equivalent strain

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0.11
#----------------------------------------------------#

#fieldscale small: dx = 1e-3 < l = 1.64e-3, 3x adaptivity levels

# [Adaptivity]
#   max_h_level = 5
#   marker = 'combo'
#   cycles_per_step = 1
#   [Markers]
#       [./combo]
#         type = FarmsComboMarker
#         markers = 'damage_marker strain_energy_marker'
#         meshsize_marker = 'meshsize_marker'
#       [../]
#       [damage_marker]
#         type = ValueThresholdMarker
#         variable = d
#         refine = 0.01
#       []
#       [strain_energy_marker]
#         type = ValueThresholdMarker
#         variable = psie_active
#         refine = '${fparse 1.0*3/8*Gc_const/l}'
#       []   
#       # if mesh_size > dxmin, refine
#       # if mesh_size < dxmin/100, coarsen (which never happens)
#       # otherwise, do nothing
#       [meshsize_marker]
#         type = ValueThresholdMarker
#         variable = mesh_size
#         refine = '${dx_min}'
#         coarsen = '${fparse dx_min/100}'
#         third_state = DO_NOTHING
#       [] 
#   []
# []

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = nonlocal_subapp.i
    cli_args = 'l=${l};kappa_i=${kappa_i};c0=${c0}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = nonlocal_eqstrain
    source_variable = nonlocal_eqstrain
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = eqstrain_local
    source_variable = eqstrain_local
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

#initial damage box 1
bottom_left1 = '-0.0025 -2e-4 0'
top_right1 = '0.0025 2e-4 0'

#initial damage box 2
bottom_left2 = '-2e-4 -0.0025 0'
top_right2 = '2e-4 0.0025 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../2dmeshfile/fieldscale_test1_2d.msh'
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
    block_id = 1
    input = extranodeset1
  []
  [./subdomain_id2]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left2}
    top_right = ${top_right2}
    location = INSIDE
    block_id = 1
    input = subdomain_id
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
  [./strength]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${fparse ft}
  [../]
  [crack_damage_aux]
    order = FIRST
    family = MONOMIAL
  []
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    order = CONSTANT
    family = MONOMIAL
  []
  [crack_damage_initial]
    family = MONOMIAL
    order = FIRST
  []
  [nonlocal_eqstrain]
      order = FIRST
      family = LAGRANGE
  [] 
  [eqstrain_local]
    family = MONOMIAL
    order = CONSTANT
  []
  [accel_x]
  []
  [accel_y]
  []
  [vel_x]
  []
  [vel_y]
  []
[]

[AuxKernels]
  #
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  #
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
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
  #damage
  # [define_initial_damage_block1]
  #   type = ConstantAux
  #   variable = crack_damage_initial
  #   value = 0.9
  #   block = 1
  #   execute_on = INITIAL
  # []
  # [define_initial_damage_block0]
  #   type = ConstantAux
  #   variable = crack_damage_initial
  #   value = 0
  #   block = '4 5'
  #   execute_on = INITIAL
  # []
  #get eqstrain_local
  [eqstrain_local_aux]
    type = MaterialRealAux
    variable = eqstrain_local
    property = eqstrain_local
    execute_on = 'TIMESTEP_END'
  []
  #get crack damage aux
  [crack_damage_aux]
    type = MaterialRealVectorValueAux
    variable = crack_damage_aux
    property = crack_damage
    component = 0
    execute_on = 'TIMESTEP_END'
  []
[]

[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 4.658e5
    shape_param_beta = 4.661e5
    rise_time = 3e-6
    single_pulse_duration = 1e-5
    EM = 0.03
    gap = 0.001
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    discharge_center = '0 0 0.0005'
    number_of_pulses = 1
    peak_pressure = 130e6 #if peak pressure is specified, the depth variation is ignored
  []
[]

[Physics/SolidMechanics/Dynamic]
  [all]
    add_variables = true
    hht_alpha = ${hht_alpha}
    newmark_beta = ${newmark_beta}
    newmark_gamma = ${newmark_gamma}
    use_automatic_differentiation = false
    # mass_damping_coefficient = 0.1
    # stiffness_damping_coefficient = 0.1
    density = ${density}
    strain = FINITE
  []
[]

# [Kernels]
#   [dispkernel_x]
#     type = DynamicStressDivergenceTensors
#     displacements = 'disp_x disp_y'
#     variable = disp_x
#     component = 0
#     zeta = 1e-8
#     use_displaced_mesh = true
#   []
#   [dispkernel_y]
#     type = DynamicStressDivergenceTensors
#     displacements = 'disp_x disp_y'
#     variable = disp_y
#     component = 1
#     zeta = 1e-8
#     use_displaced_mesh = true
#   []
#   [inertia_x]
#     type = InertialForce
#     variable = disp_x
#     velocity = vel_x
#     acceleration = accel_x
#     beta = 0.25
#     gamma = 0.5
#     use_displaced_mesh = true
#   []
#   [inertia_y]
#     type = InertialForce
#     variable = disp_y
#     velocity = vel_y
#     acceleration = accel_y
#     beta = 0.25
#     gamma = 0.5
#     use_displaced_mesh = true
#   []
# []

[BCs]
  #confinement
  [./Pressure]
    #assign pressure on inner surface
    [pressure_inner]
      boundary = 3
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
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
  #add dampers
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
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [./elastic_stress]
    type = ComputeSmearedCrackingStressDebug
    nonlocal_eqstrain = nonlocal_eqstrain
    damage_evolution_law_span = 0.5
    model = NONLOCAL
    cracking_stress = strength
    initial_crack_damage = crack_damage_initial
    output_properties = 'stress'
    softening_models = abrupt_softening
    cracked_elasticity_type = FULL
    outputs = exodus
  [../]
  # [strain]
  #   type = ComputeFiniteStrain
  #   displacements = 'disp_x disp_y'
  # []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${density}
  [] 
  [./abrupt_softening]
  type = AbruptSoftening
  [../]
  [./exponential_softening]
  type = ExponentialSoftening
  [../]  
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  # petsc_options_value = '101                asm      lu'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  # petsc_options_value = ' lu       mumps       100'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 20

  # dt = 0.5e-7
  end_time = 1e-5

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
  exodus = true
  time_step_interval = 20
  print_linear_residuals = false
  csv = true
  [checkpoint]
      type = Checkpoint
      time_step_interval = 100
      num_files = 2
  []
[]

# [Distributions]
#   #typically for granite
#   #Shape Parameter (k): 5 to 15, commonly around 8 to 12.
#   #Scale Parameter (λ): 5 to 30 MPa, commonly around 10 to 20 MPa.
#   [weibull]
#     type = Weibull
#     shape = 8.0 #k
#     scale = ${ft} #lambda
#     location = 0 
#   []
# [] 

# [ICs]
#   [./strength_var]
#     type =  RandomIC
#     variable = strength
#     distribution = weibull
#   []
# []

[UserObjects]
  [./init_sol_components]
    type = SolutionUserObject
    mesh = ./elasticity_out.e
    system_variables = 'disp_x disp_y crack_damage_aux'
    timestep = LATEST
    force_preaux = true
  [../]
[]

[ICs]
  [disp_x_ic]
    type = SolutionIC
    variable = disp_x
    solution_uo = init_sol_components
    from_variable = disp_x
  []
  [disp_y_ic]
    type = SolutionIC
    variable = disp_y
    solution_uo = init_sol_components
    from_variable = disp_y
  []
  [crack_damage_ic]
    type = SolutionIC
    variable = crack_damage_initial
    solution_uo = init_sol_components
    from_variable = crack_damage_aux
  []
[]