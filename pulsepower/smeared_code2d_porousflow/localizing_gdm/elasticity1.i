#solid properties
#----------------------------------------------------#
E = 50e9
nu = 0.373
ft = 137e6 ##computed from pf
Gc_const = 100
density = 2600
# dx_min = 5e-5
K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  1e-4 
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'
confinement_pressure  = 1e6
#----------------------------------------------------#
#gradient activity parameters
kappa_i = ${fparse ft / E}
c0 = 1e-12 #minimum value of the gradient activity parameter for the equivalent strain
#hydraulic properties
#----------------------------------------------------#
initial_pore_pressure = 0.0965e6
fluid_density = 1000
biot_coefficient = 0.4
fluid_bulk_modulus = 1e+9
viscosity = 1e-3
porosity = 0.008
solid_bulk_modulus_compliance = ${fparse 1.0/K}
# permeability = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
intrinsic_permeability = 5e-19 # m^2

##exponential permeability model
# coeff_b = 10 # coefficient for the exponential function in the effective permeability

##darcy-poiseuille permeability model: ultimate crack opening width
wc = ${fparse Gc_const / ft } # m
perm_exponent = 10 # exponent for the Darcy-Poiseuille model for the effective permeability
#----------------------------------------------------#
#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0
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

# [MultiApps]
#   [fracture]
#     type = TransientMultiApp
#     input_files = nonlocal_subapp2.i
#     cli_args = 'l=${l};kappa_i=${kappa_i};c0=${c0}'
#     execute_on = 'TIMESTEP_END'
#     clone_parent_mesh = true
#   []
# []

# [Transfers]
#   [from_d]
#     type = MultiAppCopyTransfer
#     from_multi_app = 'fracture'
#     variable = nonlocal_eqstrain
#     source_variable = nonlocal_eqstrain
#   []
#   [to_psie_active]
#     type = MultiAppCopyTransfer
#     to_multi_app = 'fracture'
#     variable = eqstrain_local
#     source_variable = eqstrain_local
#   []
# []

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = nonlocal_subapp1.i
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
    variable = 'eqstrain_local crack_damage_aux'
    source_variable = 'eqstrain_local crack_damage_aux'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator #All porous modules must contain
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
    file =  '../../2dmeshfile/fieldscale_test1_2d_small.msh'
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
    scaling = 1e-6
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-6
  [] 
  [pp]
    order = FIRST
    family = LAGRANGE  
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
    family = LAGRANGE
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
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
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
  #
  [effective_perm00_aux]
    family = MONOMIAL
    order = FIRST
  []
  [effective_perm11_aux]
    family = MONOMIAL
    order = FIRST
  []
  [effective_perm01_aux]
    family = MONOMIAL
    order = FIRST
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
  [define_initial_damage_block1]
    type = ConstantAux
    variable = crack_damage_initial
    value = 0.9
    block = 1
    execute_on = INITIAL
  []
  [define_initial_damage_block0]
    type = ConstantAux
    variable = crack_damage_initial
    value = 0
    block = '4 5'
    execute_on = INITIAL
  []
  #get eqstrain_local
  [eqstrain_local_aux]
    type = MaterialRealAux
    variable = eqstrain_local
    property = eqstrain_local
    execute_on = 'INITIAL NONLINEAR TIMESTEP_END'
  []
  #get crack damage aux
  [crack_damage_aux]
    type = MaterialRealAux
    variable = crack_damage_aux
    property = crack_damage
    execute_on = 'TIMESTEP_END'
  []
  ### PorousFlow Aux ###
  #effective permeability
  [effective_permeability_00]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 0
    variable = effective_perm00_aux
  []
  [effective_permeability_11]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 1
    column = 1
    variable = effective_perm11_aux
  []
  [effective_permeability_01]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 1
    variable = effective_perm01_aux
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
    number_of_pulses = 100
    peak_pressure = 150e6 #if peak pressure is specified, the depth variation is ignored
  []
[]

[Kernels]
  #solid
  [inertia_x]
      type = InertialForce
      variable = disp_x
      acceleration = accel_x
      velocity = vel_x
      beta = 0.25
      gamma = 0.5
      eta = 0
  []
  [inertia_y]
      type = InertialForce
      variable = disp_y
      acceleration = accel_y
      velocity = vel_y
      beta = 0.25
      gamma = 0.5
      eta = 0
  []
  [dispkernel_x]
      type = StressDivergenceTensors
      variable = disp_x
      component = 0
  []
  [dispkernel_y]
      type = StressDivergenceTensors
      variable = disp_y
      component = 1
  []
  #pressure coupling on stress tensor
  [poro_x]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_x
      component = 0
  []
  [poro_y]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_y
      component = 1
  []
  #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
  [mass0]
      type = PorousFlowFullySaturatedMassTimeDerivative
      biot_coefficient = ${biot_coefficient}
      coupling_type = HydroMechanical
      multiply_by_density = false
      variable = pp
  []
  #flux * grad(test)
  [flux]
      type = PorousFlowFullySaturatedDarcyBase
      variable = pp
      multiply_by_density = false
      gravity = '0 0 0'
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
      use_displaced_mesh = false
    []
    #assign pressure on outer surface
    [static_pressure_outer]
      boundary = 1
      factor = ${confinement_pressure}
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
  [fix_pp_gauge]
    type = DirichletBC
    variable = pp
    boundary = corner_ptr
    value = ${initial_pore_pressure}
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
  [strain]
    type = ComputeSmallStrain
  []
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [./elastic_stress]
    type = FarmsComputeSmearedCrackingStressGradsSpectral
    nonlocal_eqstrain = nonlocal_eqstrain
    paramA = 0.99
    paramB = 500
    cracking_stress = strength
    initial_crack_damage = crack_damage_initial
    output_properties = 'stress'
    outputs = exodus
    ##---------------------------------------------##
    # porous flow coupling
    ##---------------------------------------------##
    porous_flow_coupling = true
    ##-----darcy_poiseuille_permeability_model-----##
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = ${intrinsic_permeability}
    wc = ${wc}
    perm_exponent = ${perm_exponent}
    ##---------------------------------------------##
  [../]
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${density}
  [] 
  #define initial bulk modulus material property
  #check with youngs_modulus = 50e9, poissons_ratio = 0.373
  [solid_bulk_modulus_compliance]
    type = GenericConstantMaterial
    prop_names = solid_bulk_modulus_compliance
    prop_values = ${solid_bulk_modulus_compliance}
  []
  ##-------------------------------------------------------------------------##
  #porous flow related properties
  ##-------------------------------------------------------------------------##
  [temperature]
    type = PorousFlowTemperature
  []
  [eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  []
  #compute volumetric strain and its rate
  [vol_strain]
    type = PorousFlowVolumetricStrain
    outputs = exodus
  []
  #This Material is used for the fully saturated single-phase situation "
  #"where porepressure is the primary variable", saturation = 1.0
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  #List of variables that represent the mass fractions.
  #If no "variables are provided then num_phases=1=num_components."
  [massfrac]
    type = PorousFlowMassFraction
  []
  #compute porosity
  [porosity]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = ${porosity}
  []
  #comopute permeability
  [permeability] #take effective_perm
    type = ElkPorousFlowPermeabilityDamaged
  []
  #compute biot modulus #include damaged solid compliance
  [biot_modulus]
    type = ElkPorousFlowDamagedBiotModulus
    biot_coefficient = ${biot_coefficient}
    solid_bulk_compliance = ${solid_bulk_modulus_compliance}
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    output_properties = 'PorousFlow_constant_biot_modulus_qp'
    outputs = exodus
  []
  ##----------------------------------------------------------##
  #compute permeability
  # [permeability_constant]
  #     type = PorousFlowPermeabilityConst
  #     permeability = ${permeability}
  # []
  #compute biot modulus
  # [biot_modulus_constant]
  #     type = PorousFlowConstantBiotModulus
  #     biot_coefficient = ${biot_coefficient}
  #     solid_bulk_compliance = ${solid_bulk_modulus_compliance}
  #     fluid_bulk_modulus = ${fluid_bulk_modulus}
  # []  
  ##----------------------------------------------------------##
  #Compute density and viscosity
  [simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
  []
  #define relative permeability as 1 (used in PorousFlowDarcyVelocityComponent)
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
    kr = 1
  []
[]

#provide fluid properties for porous flow 
[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = ${fluid_bulk_modulus}
    density0 = ${fluid_density}
    thermal_expansion = 0
    viscosity = ${viscosity}
  []
[]

#this user object must contain for porous flow
[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [./init_sol_components]
    type = SolutionUserObject
    mesh = ./static_solve_out.e
    system_variables = 'disp_x disp_y pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
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
  [pp_ic]
    type = SolutionIC
    variable = pp
    solution_uo = init_sol_components
    from_variable = pp
  []
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0 */inertia_x */inertia_y */vel_x */vel_y */accel_x */accel_y */damp_outer_x */damp_outer_y */pressure_inner'
    start_time = 0
    end_time = 1e-8 # dt used in the simulation
  []
[../]

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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = ' lu       mumps       100'

  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  # petsc_options_value = 'gmres     hypre  boomeramg True'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 40

  # dt = 0.5e-7
  end_time = 100e-5

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
  time_step_interval = 40
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
