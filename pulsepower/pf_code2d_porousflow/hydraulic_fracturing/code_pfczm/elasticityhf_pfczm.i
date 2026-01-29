# =============================================================================
# PF-CZM Elasticity Main Application with Porous Flow Coupling
# Based on Wu's Unified Phase-Field Theory (JMPS 2017)
# =============================================================================
#
# Key PF-CZM Parameters (from Wu 2017):
#   - Crack geometric function: alpha(d) = 2d - d^2
#   - Normalization constant: c0 = pi
#   - Degradation: Rational form with explicit failure strength
#
# Material Parameters Required:
#   E    : Young's modulus [Pa]
#   nu   : Poisson's ratio [-]
#   Gc   : Fracture energy [N/m]
#   ft   : Tensile strength [Pa]
#   l    : Length scale [m]
#
# =============================================================================

# Static energy references (from static solve)
fluid_elastic_energy_total_static = 1.973476e+04
solid_elastic_energy_total_static = 4.610345e+05
solid_dissipated_energy_total_static = 1.416913e+03
full_input_energy_static = 4.807692e+05

# =============================================================================
# External Load & Pore Pressure
# =============================================================================
maximum_principal_stress = 10e6
minimum_principal_stress = 5e6

# =============================================================================
# Solid Properties
# =============================================================================
E = 30e9            # Young's modulus [Pa]
nu = 0.3            # Poisson's ratio [-]
Gc_const = 100      # Fracture energy [N/m]
solid_density = 2600 # Density [kg/m^3]

K = '${fparse E/3.0/(1.0-2.0*nu)}'    # Bulk modulus
K_s = 41.7e9        # Grain bulk modulus [Pa]
G = '${fparse E/2.0/(1.0+nu)}'        # Shear modulus
l = 2e-2            # Length scale [m]

# =============================================================================
# PF-CZM Specific Parameters
# =============================================================================
# Tensile strength - KEY PARAMETER for PF-CZM
ft = 10e6            # Tensile strength [Pa] - from material testing

# Irwin's characteristic length
# lch = '${fparse E * Gc_const / ft^2}'

# Normalization constant for alpha(d) = 2d - d^2
c0_val = 3.14159265359  # pi

# Initial slope xi of crack geometric function: xi = d(alpha)/dd at d=0
# For alpha(d) = 2d - d^2: xi = 2
# xi_val = 2

# PF-CZM constitutive parameters
# Linear softening: p=2, a2=-0.5, a3=0
# Cornelissen softening: p=2, a2=1.3868, a3=0.6567
p_deg = 2
a2 = -0.5
a3 = 0
eta = 1e-6

# Critical fracture energy for degradation function
# psic = 3*Gc/(8*l) for linear softening with optimal parameters
psic = '${fparse 3.0 * Gc_const / (8.0 * l)}'

# a1 parameter (computed internally by RationalDegradationFunction)
# a1 = (4/pi) * lch / l = Gc/(psic*c0*l/xi)
a1_check = '${fparse 4.0 / c0_val * E * Gc_const / (ft^2 * l)}'

# Wave speeds
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'

# =============================================================================
# Hydraulic Properties
# =============================================================================
fluid_density = 1000
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.05
solid_bulk_modulus_compliance = ${fparse 1.0/K}
grain_bulk_modulus = ${fparse K_s}

# Permeability by regions
intrinsic_permeability_domain = 1e-16
intrinsic_permeability_main_fractures = 1e-10
intrinsic_permeability_branch_fractures = 1e-9

# Darcy-Poiseuille model
wc = ${fparse Gc_const / ft}
perm_exponent = 10

# =============================================================================
# Finite Element Properties
# =============================================================================
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0.1

# =============================================================================
# MultiApps & Transfers
# =============================================================================
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracturehf_pfczm.i
    cli_args = 'Gc_const=${Gc_const};l=${l};psic=${psic};p_deg=${p_deg};a2=${a2};a3=${a3};eta=${eta}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
    sub_cycling = true
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
  PorousFlowDictator = dictator
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/square_with_hole_quicktest.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '10 10 0'
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
  [fluid_drainage_flux_work]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
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
  [get_pulse_load_aux]
    type = FunctionAux
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  [./max]
    type = ElementLengthAux
    variable = mesh_size
    method = max
    execute_on = TIMESTEP_BEGIN
  [../]
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
  [psie_active_enhanced_aux]
    type = MaterialRealAux
    variable = psie_active_enhanced
    property = psie_active_enhanced
    execute_on = 'TIMESTEP_END'
  []
  [biot_modulus_aux_kernel]
    type = MaterialRealAux
    variable = biot_modulus_aux
    property = PorousFlow_constant_biot_modulus_qp
    execute_on = 'TIMESTEP_END'
  []
  [biot_coefficient_aux_kernel]
    type = MaterialRealAux
    variable = biot_coefficient_aux
    property = biot_coefficient_damaged
    execute_on = 'TIMESTEP_END'
  []
  [porosity_aux_kernel]
    type = MaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'TIMESTEP_END'
  []
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
    shape_param_alpha = 5.500e+03
    shape_param_beta = 4.941e+04
    rise_time = 5e-5
    single_pulse_duration = 1e-3
    EM = 0.0025
    gap = 0.008
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    fitting_param_exponent = 0.25
    discharge_center = '0 0 0'
    number_of_pulses = 100
    base_factor = 8000
    peak_pressure = 10e6
  []
[]

[Kernels]
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
  [mass0]
      type = ElkPorousFlowFullySaturatedMassTimeDerivative
      coupling_type = HydroMechanical
      multiply_by_density = false
      variable = pp
      use_damaged_biot = true
  []
  [flux]
      type = PorousFlowFullySaturatedDarcyBase
      variable = pp
      multiply_by_density = false
      gravity = '0 0 0'
  []
[]

[BCs]
  [./Pressure]
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
    [static_pressure_top]
      boundary = top
      factor = ${minimum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_top
      save_in_disp_y = fconfinementy_top
    []
    [static_pressure_bottom]
      boundary = bottom
      factor = ${minimum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_bottom
      save_in_disp_y = fconfinementy_bottom
    []
    [static_pressure_left]
      boundary = left
      factor = ${maximum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_left
      save_in_disp_y = fconfinementy_left
    []
    [static_pressure_right]
      boundary = right
      factor = ${maximum_principal_stress}
      displacements = 'disp_x disp_y'
      use_displaced_mesh = false
      save_in_disp_x = fconfinementx_right
      save_in_disp_y = fconfinementy_right
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
  # =============================================================================
  # Elasticity Tensor
  # =============================================================================
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

  # =============================================================================
  # PF-CZM Crack Geometric Function: alpha(d) = 2d - d^2
  # =============================================================================
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d - d^2'
    phase_field = d
  []

  # =============================================================================
  # PF-CZM Parameters as Material Properties (non-AD for NDSmallDeformationIsotropicElasticity)
  # a1 = Gc/(psic*c0*l/xi) = (4/pi) * lch / l
  # Note: NDSmallDeformationIsotropicElasticity computes g internally for PF_CZM
  # =============================================================================
  [pfczm_params]
    type = GenericConstantMaterial
    prop_names = 'a1_mat a2_mat a3_mat p_mat'
    prop_values = '${a1_check} ${a2} ${a3} ${p_deg}'
  []

  # =============================================================================
  # Elasticity Model - Using NDSmallDeformationIsotropicElasticity for PF-CZM
  # with porous flow coupling and Darcy-Poiseuille permeability model
  # model_type = PF_CZM uses internal rational degradation function
  # =============================================================================
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    block = domain
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
    model_type = PF_CZM
    a1 = a1_mat
    a2 = a2_mat
    a3 = a3_mat
    p = p_mat
    eta = ${eta}
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    porous_flow_coupling = true
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = ${intrinsic_permeability_domain}
    wc = ${wc}
    perm_exponent = ${perm_exponent}
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
    model_type = PF_CZM
    a1 = a1_mat
    a2 = a2_mat
    a3 = a3_mat
    p = p_mat
    eta = ${eta}
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
    model_type = PF_CZM
    a1 = a1_mat
    a2 = a2_mat
    a3 = a3_mat
    p = p_mat
    eta = ${eta}
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

  # =============================================================================
  # PF-CZM Material Properties
  # Note: xi and c0 are computed internally by CrackGeometricFunction
  # =============================================================================
  [pfczm_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc psic'
    prop_values = '${l} ${Gc_const} ${psic}'
  []

  # =============================================================================
  # Enhanced History Energy (for porous flow coupling)
  # =============================================================================
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

  # =============================================================================
  # Solid Properties
  # =============================================================================
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${solid_density}
  []
  [solid_bulk_modulus_compliance]
    type = GenericConstantMaterial
    prop_names = solid_bulk_modulus_compliance
    prop_values = ${solid_bulk_modulus_compliance}
  []

  # =============================================================================
  # Porous Flow Properties
  # =============================================================================
  [temperature]
    type = PorousFlowTemperature
  []
  [eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
    outputs = exodus
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [porosity_damaged]
    type = ElkPorousFlowDamagedPorosity
    block = domain
    phase_field = d
    initial_porosity = ${porosity}
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
  []
  [permeability]
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
  [damaged_biot_coefficient]
    type = ElkPorousFlowDamagedBiotCoefficient
    block = domain
    phase_field = d
    solid_bulk_compliance = ${solid_bulk_modulus_compliance}
    grain_bulk_modulus = ${grain_bulk_modulus}
    minimum_degradation = ${eta}
  []
  [porosity_biot_fractures]
    type = GenericConstantMaterial
    block = 'main_fractures branch_fractures'
    prop_names = 'PorousFlow_porosity_qp_damaged biot_coefficient_damaged'
    prop_values = '1.0 1.0'
  []
  [biot_modulus]
    type = ElkPorousFlowDamagedBiotModulus
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    use_damaged_biot = true
    use_damaged_porosity = true
    output_properties = 'PorousFlow_constant_biot_modulus_qp'
    outputs = exodus
  []
  [simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
    kr = 1
  []
  [flow_fluid_driving_energy]
    type = ElkPorousFlowFluidDrivingEnergy
    use_damaged_biot = true
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = ${fluid_bulk_modulus}
    density0 = ${fluid_density}
    thermal_expansion = 0
    viscosity = ${viscosity}
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [./init_sol_components]
    type = SolutionUserObject
    mesh = ../static_solve_quicktest/elasticityhf_static_exodus.e
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
  [d_ic_domain]
    type = ConstantIC
    variable = d
    value = 0
    block = domain
  []
  [d_ic_main_fractures]
    type = ConstantIC
    variable = d
    value = 0.9
    block = main_fractures
  []
  [d_ic_branch_fractures]
    type = ConstantIC
    variable = d
    value = 0.9
    block = branch_fractures
  []
[]

[Controls]
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0 */inertia_x */inertia_y */vel_x */vel_y */accel_x */accel_y
                       */damp_top_x */damp_top_y */damp_bottom_x */damp_bottom_y
                       */damp_left_x */damp_left_y */damp_right_x */damp_right_y
                       */pressure_inner_hole */pressure_inner_hole_fractures
                       MultiApps/fracture
                       Transfers/from_d Transfers/to_psie_active
                       Transfers/pp_transfer_dissipated_energy_total Transfers/pp_transfer_dissipated_energy_first_step'
    start_time = 0
    end_time = 1.5e-6
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

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  line_search = 'none'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 50

  end_time = 100

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 1e-6
    iteration_window = 0
    cutback_factor_at_failure = 0.5
    optimal_iterations = 20
    growth_factor = 1.1
    max_time_step_bound = 2e-6
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

# =============================================================================
# Energy Calculation Postprocessors
# =============================================================================
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
      expression = 'dissipated_energy_dynamic + ${solid_dissipated_energy_total_static}'
      pp_names = 'dissipated_energy_dynamic'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

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
      expression = '-1 * external_work - confinement_work + ${full_input_energy_static} - damping_work'
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

[Postprocessors]
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + solid_dissipated_energy_total + fluid_kinetic_energy_total + fluid_elastic_energy_total + fluid_dissipated_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total solid_dissipated_energy_total fluid_kinetic_energy_total fluid_elastic_energy_total fluid_dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
