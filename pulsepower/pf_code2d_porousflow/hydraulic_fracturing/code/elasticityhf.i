fluid_elastic_energy_total_static = 0
solid_elastic_energy_total_static = 0
full_input_energy_static = 0

#external load & pore pressure
#----------------------------------------------------#
maximum_principal_stress = 10e6
minimum_principal_stress = 5e6
initial_pore_pressure = 1e6
#----------------------------------------------------#

#solid properties
#----------------------------------------------------#
E = 30e9 # Young's modulus
nu = 0.3 # Poisson's ratio
Gc_const = 100  # critical energy release rate, N * m
solid_density = 2600 # kg/m^3
K = '${fparse E/3.0/(1.0-2.0*nu)}' #bulk modulus of porous material
K_s = 41.7e9 #bulk modulus of solid grains, material property
G = '${fparse E/2.0/(1.0+nu)}'
l =  2e-2 # length scale, m
ft = '${fparse sqrt(3.0/8.0 * E*Gc_const/l)}'#137 MPa # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
fluid_density = 1000
#biot_coefficient = ${fparse 1 - K/K_s}
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.05
solid_bulk_modulus_compliance = ${fparse 1.0/K} #bulk modulus of porous medium
grain_bulk_modulus = ${fparse K_s} #solid grain bulk modulus derived from alpha: alpha = 1 - K / K_s
# permeability = '5e-19 0 0 0 5e-19 0 0 0 5e-19'

#permeability by regions (domain, main fractures, branch fractures)
intrinsic_permeability_domain = 1e-16 # m^2 #only domain accepts permeability enhancement
intrinsic_permeability_main_fractures = 1e-10 # m^2 #pe
intrinsic_permeability_branch_fractures = 1e-9 #m^2

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

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracturehf.i
    cli_args = 'Gc_const=${Gc_const};l=${l}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = 'd'
    source_variable = 'd'
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active_enhanced mesh_size'
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
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator #All porous modules must contain
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/square_with_hole.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
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
  #
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
  [fx]
  []
  [fy]
  []
  [fz]
  []
  [fx_hole_fractures]
  []
  [fy_hole_fractures]
  []
  [fz_hole_fractures]
  []
  [fconfinementx_top]
  []
  [fconfinementy_top]
  []
  [fconfinementz_top]
  []
  [fconfinementx_bottom]
  []
  [fconfinementy_bottom]
  []
  [fconfinementz_bottom]
  []
  [fconfinementx_left]
  []
  [fconfinementy_left]
  []
  [fconfinementz_left]
  []
  [fconfinementx_right]
  []
  [fconfinementy_right]
  []
  [fconfinementz_right]
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
  [fdampx_left]
  []
  [fdampy_left]
  []
  [fdampz_left]
  []
  [fdampx_right]
  []
  [fdampy_right]
  []
  [fdampz_right]
  []
  #darcy velocity components
  [darcy_vel_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [darcy_vel_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [darcy_vel_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [psie_active_enhanced]
    order = CONSTANT
    family = MONOMIAL
  []
  [biot_modulus_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [biot_coefficient_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [porosity_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  #strain components for energy calculation
  [strain_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_22]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_inc_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_inc_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_inc_22]
    order = CONSTANT
    family = MONOMIAL
  []
  # fluid drainage work density on boundary
  [fluid_drainage_flux_work]
    family = MONOMIAL
    order = CONSTANT
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
  ### Darcy Velocity
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
  #### get enhanced history energy
  [psie_active_enhanced_aux]
    type = MaterialRealAux
    variable = psie_active_enhanced
    property = psie_active_enhanced
    execute_on = 'TIMESTEP_END'
  []
  #### get biot modulus
  [biot_modulus_aux_kernel]
    type = MaterialRealAux
    variable = biot_modulus_aux
    property = PorousFlow_constant_biot_modulus_qp
    execute_on = 'TIMESTEP_END'
  []
  #### get damaged biot coefficient
  [biot_coefficient_aux_kernel]
    type = MaterialRealAux
    variable = biot_coefficient_aux
    property = biot_coefficient_damaged
    execute_on = 'TIMESTEP_END'
  []
  #### get damaged porosity
  [porosity_aux_kernel]
    type = MaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'TIMESTEP_END'
  []
  #### extract elastic strain components
  [extract_strain_00]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_00
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_11
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_22
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []
  #### extract strain increment components
  [extract_strain_inc_00]
    type = RankTwoAux
    rank_two_tensor = strain_increment
    variable = strain_inc_00
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_inc_11]
    type = RankTwoAux
    rank_two_tensor = strain_increment
    variable = strain_inc_11
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_inc_22]
    type = RankTwoAux
    rank_two_tensor = strain_increment
    variable = strain_inc_22
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []
[]

[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 0.3732
    shape_param_beta = 5.891
    rise_time = 0.5
    single_pulse_duration = 10
    EM = 0.0025
    gap = 0.008
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    fitting_param_exponent = 0.25
    discharge_center = '0 0 0'
    number_of_pulses = 100
    base_factor = 8000
    peak_pressure = 10e6 #if peak pressure is specified, the depth variation is ignored
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
      type = ElkPorousFlowEffectiveStressCoupling
      variable = disp_x
      component = 0
      use_damaged_biot = true
  []
  [poro_y]
      type = ElkPorousFlowEffectiveStressCoupling
      variable = disp_y
      component = 1
      use_damaged_biot = true
  []
  #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
  [mass0]
      type = ElkPorousFlowFullySaturatedMassTimeDerivative
      coupling_type = HydroMechanical
      multiply_by_density = false
      variable = pp
      use_damaged_biot = true
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
    [pressure_inner_hole]
      boundary = hole
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fx
      save_in_disp_y = fy
    []
    [pressure_inner_hole_fractures]
      boundary = hole_fracture
      function = func_tri_pulse
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fx_hole_fractures
      save_in_disp_y = fy_hole_fractures
    []
    #assign pressure on top surface
    [static_pressure_top]
      boundary = top
      factor = ${minimum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_top
      save_in_disp_y = fconfinementy_top
    []
    #assign pressure on bottom surface
    [static_pressure_bottom]
      boundary = bottom
      factor = ${minimum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_bottom
      save_in_disp_y = fconfinementy_bottom
    []
    #assign pressure on left surface
    [static_pressure_left]
      boundary = left
      factor = ${maximum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_left
      save_in_disp_y = fconfinementy_left
    []
    #assign pressure on right surface
    [static_pressure_right]
      boundary = right
      factor = ${maximum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_right
      save_in_disp_y = fconfinementy_right
    []
  []
  # add drained pressure
  [./porepressure_drained]
    type = FunctionDirichletBC
    variable = pp
    function = func_tri_pulse
    boundary = hole_fracture
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
  #add dampers on top surfaces
  [damp_top_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = top
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampx_top
  []
  [damp_top_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = top
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampy_top
  []
  #add dampers on bottom surfaces
  [damp_bottom_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = bottom
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampx_bottom
  []
  [damp_bottom_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = bottom
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampy_bottom
  []
  #add dampers on left surfaces
  [damp_left_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = left
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampx_left
  []
  [damp_left_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = left
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampy_left
  []
  #add dampers on right surfaces
  [damp_right_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 0
    boundary = right
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampx_right
  []
  [damp_right_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y'
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    component = 1
    boundary = right
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${solid_density}
    save_in = fdampy_right
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [strain]
    type = ComputeSmallStrain
  []
  [bulk]
    type = GenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    block = domain
    # material property names
    ##---------------------------------------------##
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_inactive = psie_inactive
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    ##---------------------------------------------##
    # decomposition type
    ##---------------------------------------------##
    decomposition = SPECTRAL
    ##---------------------------------------------##
    # model type
    ##---------------------------------------------##
    model_type = AT1
    ##---------------------------------------------##
    # constants
    ##---------------------------------------------##
    eta = 1e-6
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    ##---------------------------------------------##
    # porous flow coupling
    ##---------------------------------------------##
    porous_flow_coupling = true
    ##-----darcy_poiseuille_permeability_model-----##
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = ${intrinsic_permeability_domain}
    wc = ${wc}
    perm_exponent = ${perm_exponent}
    ##---------------------------------------------##
  []
  [elasticity_main_fractures]
    type = NDSmallDeformationIsotropicElasticity
    block = main_fractures
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_inactive = psie_inactive
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    decomposition = SPECTRAL
    model_type = AT1
    eta = 1e-6
    porous_flow_coupling = false
  []
  [elasticity_branch_fractures]
    type = NDSmallDeformationIsotropicElasticity
    block = branch_fractures
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_inactive = psie_inactive
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    decomposition = SPECTRAL
    model_type = AT1
    eta = 1e-6
    porous_flow_coupling = false
  []
  [stress]
    type = NDComputeSmallDeformationStress
    block = domain
    elasticity_model = elasticity
    output_properties = 'stress strain_increment'
    outputs = exodus
  []
  [stress_main_fractures]
    type = NDComputeSmallDeformationStress
    block = main_fractures
    elasticity_model = elasticity_main_fractures
    output_properties = 'stress strain_increment'
    outputs = exodus
  []
  [stress_branch_fractures]
    type = NDComputeSmallDeformationStress
    block = branch_fractures
    elasticity_model = elasticity_branch_fractures
    output_properties = 'stress strain_increment'
    outputs = exodus
  []
  #enhanced history energy with pressure-dependent term (from CMAME paper Appendix A)
  [history_energy_enhanced]
    type = ElkPorousFlowHistoryEnergyEnhanced
    psie_active = psie_active
    pore_pressure = pp
    initial_porosity = ${porosity}
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    bulk_modulus = K
    psie_active_enhanced = psie_active_enhanced
  []
  #solid properties
  ##-------------------------------------------------------------------------##
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
  #damage-dependent porosity (stored in *_damaged properties)
  [porosity_damaged]
    type = ElkPorousFlowDamagedPorosity
    block = domain
    phase_field = d
    initial_porosity = ${porosity}
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
  []
  #compute permeability
  [permeability] #take effective_perm (damage-dependent for domain only)
    type = ElkPorousFlowPermeabilityDamaged
    block = domain
  []
  [permeability_main_fractures]
    type = PorousFlowPermeabilityConst
    block = main_fractures
    permeability = '${intrinsic_permeability_main_fractures} 0 0  0 ${intrinsic_permeability_main_fractures} 0  0 0 ${intrinsic_permeability_main_fractures}'
  []
  [permeability_branch_fractures]
    type = PorousFlowPermeabilityConst
    block = branch_fractures
    permeability = '${intrinsic_permeability_branch_fractures} 0 0  0 ${intrinsic_permeability_branch_fractures} 0  0 0 ${intrinsic_permeability_branch_fractures}'
  []
  #damage-dependent Biot coefficient
  [damaged_biot_coefficient]
    type = ElkPorousFlowDamagedBiotCoefficient
    block = domain
    phase_field = d
    solid_bulk_compliance = ${solid_bulk_modulus_compliance}
    grain_bulk_modulus = ${grain_bulk_modulus}
    minimum_degradation = 1e-6 #this should be eta (shown in fracture.i)
  []
  #constant porosity and biot coefficient for fracture regions (fully open fractures)
  [porosity_biot_fractures]
    type = GenericConstantMaterial
    block = 'main_fractures branch_fractures'
    prop_names = 'PorousFlow_porosity_qp_damaged biot_coefficient_damaged'
    prop_values = '1.0 1.0'
  []
  #compute biot modulus #include damaged solid compliance
  [biot_modulus]
    type = ElkPorousFlowDamagedBiotModulus
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    use_damaged_biot = true
    use_damaged_porosity = true
    #porosity = ${porosity}
    output_properties = 'PorousFlow_constant_biot_modulus_qp'
    outputs = exodus
  []
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
  #Compute flow fluid driving energy
  [flow_fluid_driving_energy]
    type = ElkPorousFlowFluidDrivingEnergy
    use_damaged_biot = true
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
[]

[ICs]
[]

[Controls] # turns off time-dependent terms for the FIRST time step (steady solve with confinement only)
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0 */inertia_x */inertia_y */vel_x */vel_y */accel_x */accel_y
                       */damp_top_x */damp_top_y */damp_bottom_x */damp_bottom_y
                       */damp_left_x */damp_left_y */damp_right_x */damp_right_y
                       */pressure_inner_hole */pressure_inner_hole_fractures'
    start_time = 0
    end_time = 1e-4 # dt used in the simulation
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

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  # petsc_options_value = ' lu       mumps       100'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true
  line_search = 'bt'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 50

  # dt = 0.5e-7
  end_time = 100e-5

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 1e-4
    iteration_window = 0 #the adaptive time stepping happens at number of iterations <-> 'optimal_iterations plus/minus iteration_window'
    cutback_factor_at_failure = 0.5
    optimal_iterations = 20
    growth_factor = 1.25
    max_time_step_bound = 1e-4
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
    show = 'd vel_x vel_y vel_z pp psie_active_enhanced biot_modulus_aux biot_coefficient_aux porosity_aux'
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
    show = 'full_energy full_input_energy solid_elastic_energy_total solid_kinetic_energy_total solid_dissipated_energy_total fluid_elastic_energy_total fluid_kinetic_energy_total fluid_dissipated_energy_total damping_work confinement_work external_work dissipated_energy_first_step dissipated_energy_dynamic'
  []
[]

###############################Energy Calculation##############################

#fracture energy
###############################################################################
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
  [solid_dissipated_energy_total]
      type = ParsedPostprocessor
      expression = 'dissipated_energy_dynamic - dissipated_energy_first_step'
      pp_names = 'dissipated_energy_dynamic dissipated_energy_first_step'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# input energy
###############################################################################
[Postprocessors]
  [external_work_hole]
    type = FarmsExternalWork
    boundary = 'hole'
    forces = 'fx fy fz'
  []
  [external_work_hole_fractures]
    type = FarmsExternalWork
    boundary = 'hole_fracture'
    forces = 'fx_hole_fractures fy_hole_fractures fz_hole_fractures'
  []
  [external_work]
    type = ParsedPostprocessor
    expression = 'external_work_hole + external_work_hole_fractures'
    pp_names = 'external_work_hole external_work_hole_fractures'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Confinement work on each side
  [confinement_work_top]
    type = FarmsExternalWork
    boundary = 'top'
    forces = 'fconfinementx_top fconfinementy_top fconfinementz_top'
  []
  [confinement_work_bottom]
    type = FarmsExternalWork
    boundary = 'bottom'
    forces = 'fconfinementx_bottom fconfinementy_bottom fconfinementz_bottom'
  []
  [confinement_work_left]
    type = FarmsExternalWork
    boundary = 'left'
    forces = 'fconfinementx_left fconfinementy_left fconfinementz_left'
  []
  [confinement_work_right]
    type = FarmsExternalWork
    boundary = 'right'
    forces = 'fconfinementx_right fconfinementy_right fconfinementz_right'
  []
  [confinement_work]
    type = ParsedPostprocessor
    expression = 'confinement_work_top + confinement_work_bottom + confinement_work_left + confinement_work_right'
    pp_names = 'confinement_work_top confinement_work_bottom confinement_work_left confinement_work_right'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Damping work on each side
  [damping_work_top]
    type = FarmsExternalWork
    boundary = 'top'
    forces = 'fdampx_top fdampy_top fdampz_top'
  []
  [damping_work_bottom]
    type = FarmsExternalWork
    boundary = 'bottom'
    forces = 'fdampx_bottom fdampy_bottom fdampz_bottom'
  []
  [damping_work_left]
    type = FarmsExternalWork
    boundary = 'left'
    forces = 'fdampx_left fdampy_left fdampz_left'
  []
  [damping_work_right]
    type = FarmsExternalWork
    boundary = 'right'
    forces = 'fdampx_right fdampy_right fdampz_right'
  []
  [damping_work]
    type = ParsedPostprocessor
    expression = 'damping_work_top + damping_work_bottom + damping_work_left + damping_work_right'
    pp_names = 'damping_work_top damping_work_bottom damping_work_left damping_work_right'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [full_input_energy]
      type = ParsedPostprocessor
      #expression = '-1 * external_work - confinement_work + ${full_input_energy_static} - damping_work + fluid_drainage_work'
      expression = '-1 * external_work - confinement_work + ${full_input_energy_static} - damping_work'
      pp_names = 'external_work confinement_work damping_work'
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
      type = KineticEnergyAux
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
    type = ElementIntegralMaterialProperty
    mat_prop = psie
  []
[]

[Postprocessors]
  [solid_elastic_energy_total]
      type = ParsedPostprocessor
      expression = 'solid_elastic_energy_dynamic + ${solid_elastic_energy_total_static}'
      pp_names = 'solid_elastic_energy_dynamic'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# fluid kinetic energy
###############################################################################
[AuxVariables]
  [fluid_kinetic_energy]
      order = CONSTANT
      family = MONOMIAL
  []
[]

[AuxKernels]
  [fluid_kinetic_energy]
      type = ParsedAux
      variable = fluid_kinetic_energy
      coupled_variables = 'darcy_vel_x darcy_vel_y darcy_vel_z'
      expression = "0.5 * (darcy_vel_x * darcy_vel_x + darcy_vel_y * darcy_vel_y + darcy_vel_z * darcy_vel_z) * ${fluid_density}"
  []
[]

[Postprocessors]
  [fluid_kinetic_energy_total]
      type = ElementIntegralVariablePostprocessor
      variable = fluid_kinetic_energy
  []
[]
###############################################################################

# fluid elastic energy (using damaged biot coefficient for accounting)
###############################################################################
[AuxVariables]
  [fluid_elastic_energy]
      order = CONSTANT
      family = MONOMIAL
  []
[]

[AuxKernels]
  [get_fluid_elastic_energy]
      type = ParsedAux
      variable = fluid_elastic_energy
      coupled_variables = 'strain_00 strain_11 strain_22 pp biot_coefficient_aux'
      expression = "0.5 * biot_coefficient_aux * -pp * (strain_00+strain_11+strain_22)"
  []
[]

[Postprocessors]
  [fluid_elastic_energy_total_dynamic]
      type = ElementIntegralVariablePostprocessor
      variable = fluid_elastic_energy
  []
[]

[Postprocessors]
  [fluid_elastic_energy_total]
      type = ParsedPostprocessor
      expression = '${fluid_elastic_energy_total_static} + fluid_elastic_energy_total_dynamic'
      pp_names = 'fluid_elastic_energy_total_dynamic'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]
###############################################################################

# fluid energy dissipation
###############################################################################
[AuxVariables]
  [fluid_incremental_elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [fluid_incremental_elastic_energy_per_vol]
      type = ParsedAux
      variable = fluid_incremental_elastic_energy
      coupled_variables = 'strain_inc_00 strain_inc_11 strain_inc_22 pp biot_coefficient_aux'
      expression = "biot_coefficient_aux * -pp * (strain_inc_00+strain_inc_11+strain_inc_22)"
  []
[]

[Postprocessors]
  [fluid_incremental_elastic_energy]
      type = ElementIntegralVariablePostprocessor
      variable = fluid_incremental_elastic_energy
  []
  [fluid_incremental_elastic_energy_total]
    type = CumulativeValuePostprocessor
    postprocessor = fluid_incremental_elastic_energy
  []
  [fluid_dissipated_energy_total]
    type = ParsedPostprocessor
    pp_names = 'fluid_incremental_elastic_energy_total fluid_elastic_energy_total'
    expression = "${fluid_elastic_energy_total_static} + fluid_incremental_elastic_energy_total - fluid_elastic_energy_total"
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
###############################################################################

# Full Energy
# Note: Physics uses damaged properties (correct variational formulation from CMAME paper)
# Energy accounting uses incremental work tracking to capture property evolution effects
###############################################################################
[Postprocessors]
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + solid_dissipated_energy_total + fluid_kinetic_energy_total + fluid_elastic_energy_total + fluid_dissipated_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total solid_dissipated_energy_total fluid_kinetic_energy_total fluid_elastic_energy_total fluid_dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
