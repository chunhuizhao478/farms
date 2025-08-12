E = 50e9
nu = 0.373
ft = 25.5e6
Gc_const = 100
solid_density = 2600
dx_min = 5e-5

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  2e-4 
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
confinement_pressure  = 1e6
#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0.11
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
initial_pore_pressure = 0.0965e6
fluid_density = 1000
biot_coefficient = 0.7
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.008
permeability = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
intrinsic_permeability = 5e-19 # m^2
tortosity_value = 1.2

##exponential permeability model
# coeff_b = 10 # coefficient for the exponential function in the effective permeability

##darcy-poiseuille permeability model: ultimate crack opening width
wc = ${fparse 4 * Gc_const / ft } # m
perm_exponent = 50 # exponent for the Darcy-Poiseuille model for the effective permeability
#----------------------------------------------------#

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
    scaling = 1e-6
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-6
  []
  #pore pressure
  [p]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${initial_pore_pressure}
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
  #solid velocity
  [vel_x]
      order = FIRST
      family = LAGRANGE
  []
  [vel_y]
      order = FIRST
      family = LAGRANGE
  []    
  [vel_z]
      order = FIRST
      family = LAGRANGE        
  []
  #solid acceleration
  [accel_x]
      order = FIRST
      family = LAGRANGE
  []
  [accel_y]
      order = FIRST
      family = LAGRANGE        
  []
  [accel_z]
      order = FIRST
      family = LAGRANGE         
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
      beta = 0.25
      execute_on = timestep_end
  []
  [vel_x]
      type = NewmarkVelAux
      variable = vel_x
      acceleration = accel_x
      gamma = 0.5
      execute_on = timestep_end
  []
  [accel_y]
      type = NewmarkAccelAux
      variable = accel_y
      displacement = disp_y
      velocity = vel_y
      beta = 0.25
      execute_on = timestep_end
  []
  [vel_y]
      type = NewmarkVelAux
      variable = vel_y
      acceleration = accel_y
      gamma = 0.5
      execute_on = timestep_end
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
    execute_on = 'TIMESTEP_END'
  []
  ### PorousFlow Aux ###
  #effective permeability
  [effective_permeability_00]
    type = ADMaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 0
    variable = effective_perm00_aux
  []
  [effective_permeability_11]
    type = ADMaterialRealTensorValueAux
    property = effective_perm
    row = 1
    column = 1
    variable = effective_perm11_aux
  []
  [effective_permeability_01]
    type = ADMaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 1
    variable = effective_perm01_aux
  []
[]

[Kernels]
    #inertia force terms
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
    [] 
    #solid stress divergence (sigma * nabla u)
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = false
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
    #poro mechanical coupling (-alpha p grad du)
    #-alpha p grad du = -alpha p du/dx [i] , -alpha p du/dy [j]
    #biot coefficient is called inside
    [poromechanic_ux]
        type = ElkADPoroMechanicsCoupling
        variable = disp_x
        porepressure = p
        component = 0
        multiply_biot_coefficient = true
    []
    [poromechanic_uy]
        type = ElkADPoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
        multiply_biot_coefficient = true
    []    
    #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
    [mass0]
        type = ElkADPorousFlowFullySaturatedMassTimeDerivative
        variable = p
    []
    #flux * grad(test)
    [flux]
        type = ElkADPorousFlowFullySaturatedDarcyBase
        variable = p
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
    number_of_pulses = 10
    peak_pressure = 150e6 #if peak pressure is specified, the depth variation is ignored
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
    type = ADDirichletBC
    variable = disp_x
    boundary = corner_ptr
    value = 0
  []
  [./fix_cptr2_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = corner_ptr
    value = 0
  []
  #fix pressure
  [./fix_pressure]
    type = ADDirichletBC
    variable = p
    boundary = 3
    value = ${initial_pore_pressure}
    use_displaced_mesh = false
  []
  #add dampers
  [damp_outer_x]
    type = ADFarmsNonReflectDashpotBC
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
    type = ADFarmsNonReflectDashpotBC
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
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [strain]
    type = ADComputeSmallStrain
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
    type = SmallDeformationIsotropicElasticityHM
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active'
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
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [porodynamics] # take material property: effective_perm
    type = ElkADPoroDynamicSmearedCrackingMaterials
    rhos_value = ${solid_density}
    rhof_value = ${fluid_density}
    porosity_value = ${porosity}
    tortosity_value = ${tortosity_value} #T in paper
    viscosity_value = ${viscosity}
    bulk_modulus_solid_value = ${K}
    biot_coefficient_value = ${biot_coefficient}
    bulk_modulus_fluid_value = ${fluid_bulk_modulus}
    permeability_value = ${permeability}
  []
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0 */inertia_x */inertia_y */vel_x */vel_y */accel_x */accel_y */damp_outer_x */damp_outer_y'
    start_time = 0
    end_time = 2e-8 # dt used in the simulation
  []
[../]

[Preconditioning]
    [smp]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
      petsc_options_value = ' lu       mumps'
    []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  # petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true
  line_search = 'basic'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  # Add more iterations before failure
  nl_max_its = 30

  # dt = 0.5e-7
  end_time = 10e-5

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = true
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

#this user object must contain for porous flow
# [UserObjects]
#   [./init_sol_components]
#     type = SolutionUserObject
#     mesh = ./static_solve_out.e
#     system_variables = 'disp_x disp_y pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
#     timestep = LATEST
#     force_preaux = true
#   [../]
# []

# [ICs]
#   [disp_x_ic]
#     type = SolutionIC
#     variable = disp_x
#     solution_uo = init_sol_components
#     from_variable = disp_x
#   []
#   [disp_y_ic]
#     type = SolutionIC
#     variable = disp_y
#     solution_uo = init_sol_components
#     from_variable = disp_y
#   []
#   [pp_ic]
#     type = SolutionIC
#     variable = p
#     solution_uo = init_sol_components
#     from_variable = pp
#   []
# []