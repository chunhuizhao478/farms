# =============================================================================
# Stress-Based Driving Force Static Solve with PF-CZM Framework
# Based on Miehe & Mauthe (CMAME 2016) + Wu's PF-CZM (JMPS 2017)
# =============================================================================
#
# This file performs the static solve to establish initial conditions before
# dynamic analysis. Key differences from dynamic solve:
#   - Steady state solver (no time integration)
#   - No inertia terms
#   - No dynamic BCs (dampers, pulse loading)
#
# =============================================================================

# Static energy references (zero for initial static solve)
fluid_elastic_energy_total_static = 0
solid_elastic_energy_total_static = 0
full_input_energy_static = 0

# =============================================================================
# External Load & Pore Pressure
# =============================================================================
maximum_principal_stress = 10e6
minimum_principal_stress = 5e6
initial_pore_pressure = 1e6

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
# Stress-Based Driving Force Parameters (Miehe 2016)
# =============================================================================
# Critical stress - tensile strength
sigma_c = 10e6      # Critical fracture stress [Pa] = tensile strength

# Slope parameter controlling driving force growth rate
zeta = 1.0          # Dimensionless slope parameter

# =============================================================================
# PF-CZM Specific Parameters (Wu 2017)
# =============================================================================
ft = ${sigma_c}
c0_val = 3.14159265359  # pi

# PF-CZM constitutive parameters
p_deg = 2
a2 = -0.5
a3 = 0
eta = 1e-6

# Critical fracture energy for degradation function
psic = '${fparse 3.0 * Gc_const / (8.0 * l)}'

# a1 parameter
a1_check = '${fparse 4.0 / c0_val * E * Gc_const / (ft^2 * l)}'

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
# MultiApps & Transfers
# =============================================================================
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracturehf_static_stsbased.i
    cli_args = 'Gc_const=${Gc_const};l=${l};psic=${psic};p_deg=${p_deg};a2=${a2};a3=${a3};eta=${eta}'
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
  [to_stress_driving_force]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'H_driving mesh_size'
    source_variable = 'stress_driving_force_H mesh_size'
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
  # Stress-based driving force variables
  [max_principal_stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [mid_principal_stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [min_principal_stress]
    order = CONSTANT
    family = MONOMIAL
  []
  # Degradation function g(d) for recovering undegraded stress
  # σ̃_eff = (σ + b*p*I) / g(d)  (Miehe Eq. 95)
  [degradation_g]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_driving_force_D]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_driving_force_H]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
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
  # =============================================================================
  # Degradation Function Extraction for Undegraded Stress Recovery
  # The stress tensor from NDSmallDeformationIsotropicElasticity is degraded:
  #   σ = g(d) * σ̃_eff - b*p*I  (Miehe Eq. 95)
  # To recover undegraded effective stress for driving force (Miehe Eq. 56):
  #   σ̃_eff = (σ + b*p*I) / g(d)
  # =============================================================================
  [get_degradation_g]
    type = MaterialRealAux
    variable = degradation_g
    property = g
    execute_on = 'TIMESTEP_END'
  []
  # =============================================================================
  # Principal Stress Extraction (for stress-based driving force)
  # These are principal stresses of the DEGRADED stress tensor.
  # We recover undegraded effective stress by dividing by g(d) in the
  # stress_driving_force_D calculation.
  # =============================================================================
  [get_max_principal_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = max_principal_stress
    scalar_type = MaxPrincipal
    execute_on = 'TIMESTEP_END'
  []
  [get_mid_principal_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = mid_principal_stress
    scalar_type = MidPrincipal
    execute_on = 'TIMESTEP_END'
  []
  [get_min_principal_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = min_principal_stress
    scalar_type = MinPrincipal
    execute_on = 'TIMESTEP_END'
  []
  # =============================================================================
  # Stress-Based Driving Force D (Miehe Eq. 56)
  # D = zeta * < SUM_a (<σ̃_eff^a>_+ / sigma_c)^2 - 1 >_+
  # =============================================================================
  # IMPORTANT: Must use UNDEGRADED effective stress σ̃_eff (Miehe Eq. 56, 95)
  #
  # The stress tensor from elasticity model is degraded:
  #   σ = g(d) * σ̃_eff - b*p*I  (Miehe Eq. 95)
  #
  # To recover undegraded effective stress:
  #   σ̃_eff = (σ + b*p*I) / g(d)
  # =============================================================================
  [compute_stress_driving_force_D]
    type = ParsedAux
    variable = stress_driving_force_D
    coupled_variables = 'max_principal_stress mid_principal_stress min_principal_stress pp biot_coefficient_aux degradation_g'
    constant_names = 'sigma_c zeta g_min'
    constant_expressions = '${sigma_c} ${zeta} 1e-10'
    expression = 'g_safe := max(degradation_g, g_min);
                  sig_eff1_undeg := (max_principal_stress + biot_coefficient_aux * pp) / g_safe;
                  sig_eff2_undeg := (mid_principal_stress + biot_coefficient_aux * pp) / g_safe;
                  sig_eff3_undeg := (min_principal_stress + biot_coefficient_aux * pp) / g_safe;
                  sig1_pos := max(sig_eff1_undeg, 0);
                  sig2_pos := max(sig_eff2_undeg, 0);
                  sig3_pos := max(sig_eff3_undeg, 0);
                  sum_sq := (sig1_pos/sigma_c)^2 + (sig2_pos/sigma_c)^2 + (sig3_pos/sigma_c)^2;
                  zeta * max(sum_sq - 1, 0)'
    execute_on = 'TIMESTEP_END'
  []
  # =============================================================================
  # History Variable H = max(H_old, D) for irreversibility
  # =============================================================================
  [compute_stress_driving_force_H]
    type = ParsedAux
    variable = stress_driving_force_H
    coupled_variables = 'stress_driving_force_D'
    use_xyzt = false
    expression = 'stress_driving_force_D'
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
    peak_pressure = 5e6
  []
[]

[Kernels]
  # Static solve - no inertia terms
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
  # Static solve - no mass time derivative
  [flux]
      type = PorousFlowFullySaturatedDarcyBase
      variable = pp
      multiply_by_density = false
      gravity = '0 0 0'
  []
[]

[BCs]
  [./Pressure]
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
  # PF-CZM Crack Geometric Function
  # =============================================================================
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d - d^2'
    phase_field = d
  []

  # =============================================================================
  # PF-CZM Parameters
  # =============================================================================
  [pfczm_params]
    type = GenericConstantMaterial
    prop_names = 'a1_mat a2_mat a3_mat p_mat'
    prop_values = '${a1_check} ${a2} ${a3} ${p_deg}'
  []

  # =============================================================================
  # PF-CZM Material Properties
  # =============================================================================
  [pfczm_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc psic'
    prop_values = '${l} ${Gc_const} ${psic}'
  []

  # =============================================================================
  # Elasticity Model
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
  # Enhanced History Energy
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
[]

[ICs]
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
  type = Steady

  solve_type = NEWTON

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  line_search = 'bt'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 50

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  [./exodus]
    type = Exodus
    time_step_interval = 1
    show = 'disp_x disp_y pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22 d max_principal_stress stress_driving_force_H'
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
      expression = 'dissipated_energy_dynamic - dissipated_energy_first_step'
      pp_names = 'dissipated_energy_dynamic dissipated_energy_first_step'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [external_work_hole]
    type = FarmsExternalWork
    boundary = 'hole'
    forces = 'fx fy fz'
    use_displacement_work = true
  []
  [external_work_hole_fractures]
    type = FarmsExternalWork
    boundary = 'hole_fracture'
    forces = 'fx_hole_fractures fy_hole_fractures fz_hole_fractures'
    use_displacement_work = true
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
    use_displacement_work = true
  []
  [confinement_work_bottom]
    type = FarmsExternalWork
    boundary = 'bottom'
    forces = 'fconfinementx_bottom fconfinementy_bottom fconfinementz_bottom'
    use_displacement_work = true
  []
  [confinement_work_left]
    type = FarmsExternalWork
    boundary = 'left'
    forces = 'fconfinementx_left fconfinementy_left fconfinementz_left'
    use_displacement_work = true
  []
  [confinement_work_right]
    type = FarmsExternalWork
    boundary = 'right'
    forces = 'fconfinementx_right fconfinementy_right fconfinementz_right'
    use_displacement_work = true
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
    use_displacement_work = true
  []
  [damping_work_bottom]
    type = FarmsExternalWork
    boundary = 'bottom'
    forces = 'fdampx_bottom fdampy_bottom fdampz_bottom'
    use_displacement_work = true
  []
  [damping_work_left]
    type = FarmsExternalWork
    boundary = 'left'
    forces = 'fdampx_left fdampy_left fdampz_left'
    use_displacement_work = true
  []
  [damping_work_right]
    type = FarmsExternalWork
    boundary = 'right'
    forces = 'fdampx_right fdampy_right fdampz_right'
    use_displacement_work = true
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
