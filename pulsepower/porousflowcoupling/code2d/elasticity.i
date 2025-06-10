#solid properties
#----------------------------------------------------#
E = 50e9 # Young's modulus
nu = 0.373 # Poisson's ratio
Gc_const = 57  # critical energy release rate, N * m
solid_density = 2600 # kg/m^3 
dx_min = 2.5e-5 # minimum mesh size, m
K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  1e-4 # length scale, m
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
confinement_pressure  = 1e6
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
fluid_density = 1000
biot_coefficient = 0.7
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.008
solid_bulk_modulus_compliance = 1.524e-11
# permeability = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
intrinsic_permeability = 5e-19 # m^2
coeff_b = 5.0 # coefficient for the exponential function in the effective permeability
#----------------------------------------------------#

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0.11
#----------------------------------------------------#

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
      [../]
      [damage_marker]
        type = ValueThresholdMarker
        variable = d
        refine = 0.5
      []
      [strain_energy_marker]
        type = ValueThresholdMarker
        variable = psie_active
        refine = '${fparse 1.0*3/8*Gc_const/l}'
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
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator #All porous modules must contain
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../2dmeshfile/fieldscale_test1_2d.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE  
    scaling = 1e-6
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE  
    scaling = 1e-6
  []
  [pp]
    order = FIRST
    family = LAGRANGE  
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
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
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

[Physics/SolidMechanics/Dynamic]
  [all]
    add_variables = true
    hht_alpha = ${hht_alpha}
    newmark_beta = ${newmark_beta}
    newmark_gamma = ${newmark_gamma}
    density = ${solid_density}
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
    gap = 0.001
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    discharge_center = '0 0 0.0005'
    number_of_pulses = 1
    peak_pressure = 100e6 #if peak pressure is specified, the depth variation is ignored
  []
  #strain
  [func_strain_xx]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_00
  [../]
  [func_strain_xy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_01
  [../]
  [func_strain_xz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_02
  [../]
  [func_strain_yy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_11
  [../]
  [func_strain_yz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_12
  [../]
  [func_strain_zz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = elastic_strain_22
  [../]
[]

[Kernels]
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
    variable = pp
  []
  #flux * grad(test)
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = pp
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
[]

[BCs]
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
    density = ${solid_density}
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
    density = ${solid_density}
  []
[]

[Materials]
  [bulk]
    type = GenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    # material property names
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    # decomposition type
    decomposition = SPECTRAL
    # model type
    model_type = AT1
    # constants
    eta = 1e-6
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    # porous flow coupling
    porous_flow_coupling = true
    intrinsic_permeability = ${intrinsic_permeability}
    coeff_b = ${coeff_b}
  []
  [stress]
    type = NDComputeSmallDeformationStress ###
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  #solid properties
  ##-------------------------------------------------------------------------##
  [./initial_strain]
    type = GenericFunctionRankTwoTensor
    tensor_name = static_initial_strain_tensor
    tensor_functions = 'func_strain_xx     func_strain_xy      func_strain_xz 
                        func_strain_xy     func_strain_yy      func_strain_yz
                        func_strain_xz     func_strain_yz      func_strain_zz'
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${solid_density}
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
  # #compute biot modulus #include damaged solid compliance
  # [biot_modulus]
  #   type = ElkPorousFlowDamagedBiotModulus
  #   biot_coefficient = ${biot_coefficient}
  #   solid_bulk_compliance = ${solid_bulk_modulus_compliance}
  #   fluid_bulk_modulus = ${fluid_bulk_modulus}
  # []
  ##----------------------------------------------------------##
  # #compute permeability
  # [permeability_constant]
  #     type = PorousFlowPermeabilityConst
  #     permeability = ${permeability}
  # []
  #compute biot modulus
  [biot_modulus_constant]
      type = PorousFlowConstantBiotModulus
      biot_coefficient = ${biot_coefficient}
      solid_bulk_compliance = ${solid_bulk_modulus_compliance}
      fluid_bulk_modulus = ${fluid_bulk_modulus}
  []  
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

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0'
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

  solve_type = NEWTON

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  # dt = 0.5e-7
  end_time = 6e-5

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
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
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