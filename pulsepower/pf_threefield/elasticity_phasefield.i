# =============================================================================
# Three-Field Poroelastodynamics with Phase-Field Damage (2D)
# =============================================================================
# Formulation: Three-field (u, w, p) with AD damage-dependent hydraulic properties
#
# Variables:
#   - disp_x, disp_y: solid skeleton displacement
#   - wf_x, wf_y: fluid relative displacement
#   - p: pore pressure
#
# Coupling: MultiApp phase-field damage with staggered scheme
# =============================================================================

# Static energy values from static solve (to be updated)
fluid_elastic_energy_total_static = 8.082010e-05
solid_elastic_energy_total_static = 5.411602e-03
full_input_energy_static = 5.492422e-03

#solid properties
#----------------------------------------------------#
E = 50e9 # Young's modulus
nu = 0.3 # Poisson's ratio
Gc_const = 40  # critical energy release rate, N/m
solid_density = 2600 # kg/m^3
K = '${fparse E/3.0/(1.0-2.0*nu)}' #bulk modulus of porous material
K_s = 50e9 #bulk modulus of solid grains
G = '${fparse E/2.0/(1.0+nu)}'
l = 2e-4 # length scale, m
ft = '${fparse sqrt(3.0/8.0 * E*Gc_const/l)}'
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
confinement_pressure = 1000000.0
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
fluid_density = 1000
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.008
solid_bulk_modulus_compliance = ${fparse 1.0/K}
grain_bulk_modulus = ${fparse K_s}
intrinsic_permeability = 5e-19 # m^2
tortosity = 1.2

# Darcy-Poiseuille permeability model
wc = ${fparse Gc_const / ft} # m
perm_exponent = 10
#----------------------------------------------------#

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0
#----------------------------------------------------#

# =============================================================================
# MultiApp for Phase-Field
# =============================================================================
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
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

# =============================================================================
# Global Parameters
# =============================================================================
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# =============================================================================
# Mesh
# =============================================================================
#initial damage box
bottom_left1 = '-0.0025 -3e-4 0'
top_right1 = '0.0025 3e-4 0'
bottom_left2 = '-3e-4 -0.0025 0'
top_right2 = '3e-4 0.0025 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = './mesh/fieldscale_test1_2d.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
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

# =============================================================================
# Variables
# =============================================================================
[Variables]
  # Solid displacement
  [disp_x]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-6
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-6
  []
  # Fluid relative displacement
  [wf_x]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-6
  []
  [wf_y]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-6
  []
  # Pore pressure
  [p]
    order = FIRST
    family = LAGRANGE
  []
[]

# =============================================================================
# Auxiliary Variables
# =============================================================================
[AuxVariables]
  # Phase-field damage
  [d]
    family = LAGRANGE
    order = FIRST
  []
  # Solid velocity
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [vel_z]
    # Dummy z-velocity for 2D (required by KineticEnergyAux)
    family = LAGRANGE
    order = FIRST
  []
  # Solid acceleration
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  # Fluid velocity
  [vf_x]
    family = LAGRANGE
    order = FIRST
  []
  [vf_y]
    family = LAGRANGE
    order = FIRST
  []
  # Fluid acceleration
  [af_x]
    family = LAGRANGE
    order = FIRST
  []
  [af_y]
    family = LAGRANGE
    order = FIRST
  []
  # Mesh size for phase-field
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  # Enhanced history energy (with pressure)
  [psie_active_enhanced]
    order = CONSTANT
    family = MONOMIAL
  []
  # Material property outputs
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
  [effective_perm00_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  # Darcy velocity components
  [darcy_vel_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [darcy_vel_y]
    order = CONSTANT
    family = MONOMIAL
  []
  # Strain components for energy calculation
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
  # Pulse load
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  # Force components
  [fx]
    order = SECOND
    family = LAGRANGE
  []
  [fy]
    order = SECOND
    family = LAGRANGE
  []
  [fconfinementx]
    order = SECOND
    family = LAGRANGE
  []
  [fconfinementy]
    order = SECOND
    family = LAGRANGE
  []
  [fdampx]
    order = SECOND
    family = LAGRANGE
  []
  [fdampy]
    order = SECOND
    family = LAGRANGE
  []
[]

# =============================================================================
# Kernels
# =============================================================================
[Kernels]
  # --- Solid Momentum Equation ---
  # Inertia: rho * a^s
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
  # Pore fluid inertia coupling: rho^f * a^f
  [porefluidinertia_x]
    type = ElkADPoreFluidInertialForceCoupling
    variable = disp_x
    fluiddisp = wf_x
    fluidvel = vf_x
    fluidaccel = af_x
    beta = ${newmark_beta}
  []
  [porefluidinertia_y]
    type = ElkADPoreFluidInertialForceCoupling
    variable = disp_y
    fluiddisp = wf_y
    fluidvel = vf_y
    fluidaccel = af_y
    beta = ${newmark_beta}
  []
  # Stress divergence: div(sigma)
  [dispkernel_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
  []
  [dispkernel_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
  # Poromechanics coupling: -alpha * p * grad(test)
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

  # --- Fluid Momentum Equation (Dynamic Darcy Flow) ---
  # rho^f * a^s + rho^f * tau_t / phi * a^f + mu_f / kappa * v^f
  [dynamicdarcyflow_x]
    type = ElkADDynamicDarcyFlowPhaseField
    variable = wf_x
    skeletondisp = disp_x
    skeletonvel = vel_x
    skeletonaccel = accel_x
    fluidvel = vf_x
    fluidaccel = af_x
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    component = 0
  []
  [dynamicdarcyflow_y]
    type = ElkADDynamicDarcyFlowPhaseField
    variable = wf_y
    skeletondisp = disp_y
    skeletonvel = vel_y
    skeletonaccel = accel_y
    fluidvel = vf_y
    fluidaccel = af_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    component = 1
  []
  # Pressure coupling on fluid: -p * grad(test)
  [poromechanic_wx]
    type = ElkADPoroMechanicsCoupling
    variable = wf_x
    porepressure = p
    component = 0
    multiply_biot_coefficient = false
  []
  [poromechanic_wy]
    type = ElkADPoroMechanicsCoupling
    variable = wf_y
    porepressure = p
    component = 1
    multiply_biot_coefficient = false
  []

  # --- Mass Conservation Equation ---
  # Pressure rate: p_dot / M
  [pressure_rate]
    type = ElkADPressureRate
    variable = p
  []
  # Solid skeleton mass conservation: alpha * M * div(v^s)
  [massconservationskeleton]
    type = ElkADMassConservationNewmark
    variable = p
    displacement_x = disp_x
    displacement_y = disp_y
    velocity_x = vel_x
    velocity_y = vel_y
    acceleration_x = accel_x
    acceleration_y = accel_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    multiply_biot_coefficient = true
    plane_strain_correction = true
    poissons_ratio = ${nu}
  []
  # Fluid mass conservation: M * div(v^f)
  [insmass]
    type = ElkADMassConservationNewmark
    variable = p
    displacement_x = wf_x
    displacement_y = wf_y
    velocity_x = vf_x
    velocity_y = vf_y
    acceleration_x = af_x
    acceleration_y = af_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    multiply_biot_coefficient = false
    plane_strain_correction = true
    poissons_ratio = ${nu}
  []
[]

# =============================================================================
# Auxiliary Kernels
# =============================================================================
[AuxKernels]
  # Solid velocity and acceleration (Newmark)
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
  # Fluid velocity and acceleration (Newmark)
  [af_x]
    type = NewmarkAccelAux
    variable = af_x
    displacement = wf_x
    velocity = vf_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vf_x]
    type = NewmarkVelAux
    variable = vf_x
    acceleration = af_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  [af_y]
    type = NewmarkAccelAux
    variable = af_y
    displacement = wf_y
    velocity = vf_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vf_y]
    type = NewmarkVelAux
    variable = vf_y
    acceleration = af_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  # Mesh size for phase-field
  [mesh_size_aux]
    type = ElementLengthAux
    variable = mesh_size
    method = max
    execute_on = TIMESTEP_BEGIN
  []
  # Pulse load aux
  [get_pulse_load_aux]
    type = FunctionAux
    variable = pulse_load_aux
    function = func_tri_pulse
    execute_on = timestep_end
  []
  # Material property outputs
  [psie_active_enhanced_aux]
    type = ADMaterialRealAux
    variable = psie_active_enhanced
    property = psie_active_enhanced
    execute_on = 'TIMESTEP_END'
  []
  [biot_modulus_aux_kernel]
    type = ADMaterialRealAux
    variable = biot_modulus_aux
    property = biot_modulus
    execute_on = 'TIMESTEP_END'
  []
  [biot_coefficient_aux_kernel]
    type = ADMaterialRealAux
    variable = biot_coefficient_aux
    property = biot_coefficient
    execute_on = 'TIMESTEP_END'
  []
  [porosity_aux_kernel]
    type = ADMaterialRealAux
    variable = porosity_aux
    property = porosity
    execute_on = 'TIMESTEP_END'
  []
  [effective_permeability_00]
    type = ADRankTwoAux
    rank_two_tensor = permeability
    variable = effective_perm00_aux
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  # Strain components (AD version)
  [extract_strain_00]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_00
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_11]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_11
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [extract_strain_22]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_22
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []
[]

# =============================================================================
# Functions
# =============================================================================
[Functions]
  [func_tri_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = 4.658e5
    shape_param_beta = 4.661e5
    rise_time = 3e-6
    single_pulse_duration = 1e-5
    EM = 0.00125
    gap = 0.008
    convert_efficiency = 1.0
    fitting_param_alpha = 0.35
    fitting_param_exponent = 0.25
    discharge_center = '0 0 0'
    number_of_pulses = 100
    base_factor = 8000
  []
[]

# =============================================================================
# Boundary Conditions
# =============================================================================
[BCs]
  # Pulse power loading on inner surface
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
  # Drained pressure on inner surface
  #[./porepressure_drained]
  #  type = FunctionDirichletBC
  #  variable = p
  #  function = func_tri_pulse
  #  boundary = 3
  #[]
  # Fix corner point
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
  # Dampers on outer boundary
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
    density = ${solid_density}
    save_in = fdampy
  []
[]

# =============================================================================
# Materials
# =============================================================================
[Materials]
  # Elasticity tensor (AD version)
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [strain]
    type = ADComputeSmallStrain
  []
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []

  # AD Phase-field elasticity with damage-dependent permeability
  [elasticity]
    type = ADSmallDeformationIsotropicElasticityPF
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
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    # Porous flow coupling
    porous_flow_coupling = true
    # Darcy-Poiseuille permeability model
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = ${intrinsic_permeability}
    wc = ${wc}
    perm_exponent = ${perm_exponent}
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []

  # AD Three-field enhanced history energy with pressure, fluid divergence, and kinetic terms
  # Two-field H+:   H+ = max⟨ 2ψ_o^{e+} + Γp² ⟩^+
  # Three-field H+: H+ = max⟨ 2ψ_o^{e+} + Γp² - ΓM²(∇·w)² + (1-φ_o)ρ^f τ_t/φ² |ẇ|² - (1-φ_o)(ρ^f-ρ^s)|u̇|² ⟩^+
  [history_energy_enhanced]
    type = ElkADThreeFieldHistoryEnergyEnhanced
    psie_active = psie_active
    pore_pressure = p
    # Solid velocity components
    solid_velocity_x = vel_x
    solid_velocity_y = vel_y
    # Fluid velocity components
    fluid_velocity_x = vf_x
    fluid_velocity_y = vf_y
    # Fluid displacement (for divergence calculation)
    fluid_disp_x = wf_x
    fluid_disp_y = wf_y
    # Material properties
    initial_porosity = ${porosity}
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    bulk_modulus = K
    # Density and tortosity
    fluid_density = ${fluid_density}
    solid_density = ${solid_density}
    tortosity = ${tortosity}
    # Output property name
    psie_active_enhanced = psie_active_enhanced
    # Formulation options (for testing pressure oscillation)
    # Set to true to use two-field H+ (only elastic + pressure terms)
    use_two_field_formulation = true #for testing
    # Or control individual terms:
    # include_fluid_divergence_term = true
    # include_fluid_kinetic_term = true
    # include_density_diff_kinetic_term = true
  []

  # Note: density is provided by porodynamics material (damage-dependent)

  # Initial bulk modulus compliance - AD version
  [solid_bulk_modulus_compliance]
    type = ADGenericConstantMaterial
    prop_names = solid_bulk_modulus_compliance
    prop_values = ${solid_bulk_modulus_compliance}
  []

  # --- AD Damage-Dependent Hydraulic Properties ---
  # Note: effective_perm is provided by the elasticity material through
  # the Darcy-Poiseuille permeability model (damage-dependent)

  # Damaged Biot coefficient: alpha(d) = 1 - g(d) * K0 / Ks
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = ${solid_bulk_modulus_compliance}
    grain_bulk_modulus = ${grain_bulk_modulus}
    minimum_degradation = 1e-6
  []

  # Damaged porosity: phi(d) = phi_0 + (1 - phi_0) * (1 - g(d))
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = ${porosity}
    porosity_lower_bound = ${porosity}
    porosity_upper_bound = 0.999
  []

  # Damaged Biot modulus: 1/M = phi/Kf + (alpha - phi)/Ks
  [damaged_biot_modulus]
    type = ElkADPorousFlowDamagedBiotModulus
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    use_damaged_biot = true
    use_damaged_porosity = true
  []

  # Poro-dynamic material assembly (three-field with phase-field)
  [porodynamics]
    type = ElkADPoroDynamicPhaseFieldMaterials
    rhos_value = ${solid_density}
    rhof_value = ${fluid_density}
    tortosity_value = ${tortosity}
    viscosity_value = ${viscosity}
    grain_bulk_modulus = ${grain_bulk_modulus}
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    use_damaged_properties = true
    output_properties = 'biot_modulus biot_coefficient porosity density'
    outputs = exodus
  []
[]

# =============================================================================
# UserObjects for Solution Loading
# =============================================================================
[UserObjects]
  [./init_sol_components]
    type = SolutionUserObject
    mesh = ./static_solve_out.e
    system_variables = 'disp_x disp_y pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
    timestep = LATEST
    force_preaux = true
  [../]
[]

# =============================================================================
# Initial Conditions (from static solve)
# =============================================================================
[ICs]
  # Load displacement from static solve
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
  # Load pressure from static solve
  [p_ic]
    type = SolutionIC
    variable = p
    solution_uo = init_sol_components
    from_variable = pp
  []
  # Fluid relative displacement is zero at static equilibrium
  [wf_x_ic]
    type = ConstantIC
    variable = wf_x
    value = 0
  []
  [wf_y_ic]
    type = ConstantIC
    variable = wf_y
    value = 0
  []
[]

# =============================================================================
# Preconditioning
# =============================================================================
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

# =============================================================================
# Executioner
# =============================================================================
[Executioner]
  type = Transient
  solve_type = NEWTON  # Jacobian-free Newton-Krylov (much faster, no Jacobian assembly)

  # Iterative solver with AMG preconditioner
  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'gmres     hypre    boomeramg      300                0.7'

  # Line search
  line_search = 'bt'
  l_max_its = 300
  l_tol = 1e-4

  # Nonlinear solver tolerances
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 25
  nl_div_tol = 1e10

  # Reuse preconditioner to reduce Jacobian computations
  reuse_preconditioner = true

  end_time = 100e-5

  # Fixed point iteration for staggered MultiApp coupling
  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-4
  fixed_point_abs_tol = 1e-6

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = 1e-8
    iteration_window = 2
    cutback_factor_at_failure = 0.5
    optimal_iterations = 12
    growth_factor = 1.2
    max_time_step_bound = 1e-8
  []
  [TimeIntegrator]
    type = NewmarkBeta
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
[]

# =============================================================================
# Controls
# =============================================================================
[Controls]
  [./period0]
    type = TimePeriod
    disable_objects = '*/inertia_x */inertia_y */vel_x */vel_y */accel_x */accel_y */damp_outer_x */damp_outer_y */pressure_inner */porefluidinertia_x */porefluidinertia_y */dynamicdarcyflow_x */dynamicdarcyflow_y'
    start_time = 0
    end_time = 1e-8
  []
[../]

# =============================================================================
# Outputs
# =============================================================================
[Outputs]
  [./exodus]
    type = Exodus
    time_step_interval = 5
    show = 'd vel_x vel_y vf_x vf_y p psie_active_enhanced biot_modulus_aux biot_coefficient_aux porosity_aux'
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
    show = 'full_energy full_input_energy solid_elastic_energy_total solid_kinetic_energy_total solid_dissipated_energy_total fluid_elastic_energy_total fluid_kinetic_energy_total'
  []
[]

# =============================================================================
# Postprocessors (Energy Calculation)
# =============================================================================

# Fracture energy (received from subapp)
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

# Input energy
[Postprocessors]
  [external_work]
    type = FarmsExternalWork
    boundary = '3'
    forces = 'fx fy'
  []
  [confinement_work]
    type = FarmsExternalWork
    boundary = '1'
    forces = 'fconfinementx fconfinementy'
  []
  [damping_work]
    type = FarmsExternalWork
    boundary = '1'
    forces = 'fdampx fdampy'
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

# Solid kinetic energy
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

# Solid elastic energy (use AuxVariable approach for AD material property)
[AuxVariables]
  [solid_elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [solid_elastic_energy_aux]
    type = ADMaterialRealAux
    variable = solid_elastic_energy
    property = psie
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [solid_elastic_energy_dynamic]
    type = ElementIntegralVariablePostprocessor
    variable = solid_elastic_energy
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

# Fluid kinetic energy (three-field: uses fluid velocity)
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
    coupled_variables = 'vf_x vf_y'
    expression = "0.5 * (vf_x * vf_x + vf_y * vf_y) * ${fluid_density}"
  []
[]

[Postprocessors]
  [fluid_kinetic_energy_total]
    type = ElementIntegralVariablePostprocessor
    variable = fluid_kinetic_energy
  []
[]

# Fluid elastic energy
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
    coupled_variables = 'strain_00 strain_11 strain_22 p biot_coefficient_aux'
    expression = "0.5 * biot_coefficient_aux * -p * (strain_00+strain_11+strain_22)"
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

# Full Energy
[Postprocessors]
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + solid_dissipated_energy_total + fluid_kinetic_energy_total + fluid_elastic_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total solid_dissipated_energy_total fluid_kinetic_energy_total fluid_elastic_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
