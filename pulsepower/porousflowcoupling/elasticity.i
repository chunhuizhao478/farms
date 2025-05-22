# dynamic time integration parameters
#----------------------------------------------------#
beta = 0.25
gamma = 0.5
hht_alpha = 0.11
#----------------------------------------------------#

#boundary condition + phase field parameters
#----------------------------------------------------#
confinement_pressure  = 1e6
l = 2e-4 # N * h
#----------------------------------------------------#

#solid properties
#----------------------------------------------------#
E = 50e9
nu = 0.373
solid_density = 2600
Gc_const = 57
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
fluid_density = 1000
biot_coefficient = 0.7
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
#----------------------------------------------------#

#convert properties
#----------------------------------------------------#
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
#----------------------------------------------------#

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
  PorousFlowDictator = dictator #All porous modules must contain
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../meshfile/debug_sample.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '-0.00570169 -0.00612895 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE  
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE  
  []
  [disp_z]
    order = FIRST
    family = LAGRANGE  
  []
  #pore pressure
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
    beta = ${beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${gamma}
    execute_on = 'TIMESTEP_END'
  []
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${gamma}
    execute_on = 'TIMESTEP_END'
  []
  [accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = ${beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = ${gamma}
    execute_on = 'TIMESTEP_END'
  []
  #get pulse load aux
  [get_pulse_load_aux]
    type = FunctionAux 
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  ### PorousFlow Aux ###
  #biot modulus
  [biot_modulus]
    type = MaterialRealAux
    property = PorousFlow_constant_biot_modulus_qp
    variable = biot_modulus_aux
    execute_on = 'TIMESTEP_END'
  []
  #effective permeability
  [effective_permeability_00]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 0
    column = 0
    variable = effective_perm00_aux
  []
  [effective_permeability_11]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 1
    column = 1
    variable = effective_perm11_aux
  []
  [effective_permeability_22]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 2
    column = 2
    variable = effective_perm22_aux
  []
  [effective_permeability_01]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 0
    column = 1
    variable = effective_perm01_aux
  []
  [effective_permeability_02]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 0
    column = 2
    variable = effective_perm02_aux
  []
  [effective_permeability_12]
    type = MaterialRealTensorValueAux
    property = effective_perm_smeared_crack
    row = 1
    column = 2
    variable = effective_perm12_aux
  []  
  #darcy velocity
  [bulk_vel_x]
    type = PorousFlowDarcyVelocityComponent
    variable = darcy_vel_x
    component = x
    fluid_phase = 0
    gravity = '0 0 0'
  []
  [bulk_vel_y]
    type = PorousFlowDarcyVelocityComponent
    variable = darcy_vel_y
    component = y
    fluid_phase = 0
    gravity = '0 0 0'
  []
  [bulk_vel_z]
    type = PorousFlowDarcyVelocityComponent
    variable = darcy_vel_z
    component = z
    fluid_phase = 0
    gravity = '0 0 0'
  []
[]

[Physics/SolidMechanics/Dynamic]
  [all]
    add_variables = true
    hht_alpha = ${hht_alpha}
    newmark_beta = ${beta}
    newmark_gamma = ${gamma}
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
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = ${biot_coefficient}
    variable = disp_z
    component = 2
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
      boundary = 4
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = true
    [] 
    #assign pressure on top surface
    [static_pressure_top]
      boundary = 2
      factor = ${confinement_pressure}
      displacements = 'disp_x disp_y disp_z'
      use_displaced_mesh = false
    []
    #assign pressure on outer surface
    [static_pressure_outer]
      boundary = 3
      factor = ${confinement_pressure}
      displacements = 'disp_x disp_y disp_z'
      use_displaced_mesh = false
    []  
    [static_pressure_bottom]
      boundary = 5
      factor = ${confinement_pressure}
      displacements = 'disp_x disp_y disp_z'
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
  [./fix_cptr3_z]
    type = DirichletBC
    variable = disp_z
    boundary = corner_ptr
    value = 0
  []
[]

[BCs]
  #add dampers
  [damp_outer_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = 3
    beta = ${beta}
    gamma = ${gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${fparse Cs}
    p_wave_speed = ${fparse Cp}
    density = ${solid_density}
  []
  [damp_outer_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = 3
    beta = ${beta}
    gamma = ${gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${fparse Cs}
    p_wave_speed = ${fparse Cp}
    density = ${solid_density}
  []
  [damp_outer_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = 3
    beta = ${beta}
    gamma = ${gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${fparse Cs}
    p_wave_speed = ${fparse Cp}
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
    model_type = AT2
    # constants
    eta = 1e-6 
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [stress]
    type = NDComputeSmallDeformationStress ###
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
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
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [./init_sol_components]
    type = SolutionUserObject
    mesh = ./static_solve_out.e
    system_variables = 'disp_x disp_y disp_z pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
    timestep = LATEST
    force_preaux = true
  [../]
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
  end_time = 1000

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
  time_step_interval = 40
  print_linear_residuals = false
  csv = true
[]

