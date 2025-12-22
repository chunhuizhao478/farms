

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 500
  ny = 250
  xmin = -40
  xmax = 40
  ymin = 0
  ymax = 40
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  gravity = '0 0 0'
  PorousFlowDictator = dictator
  block = 0
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [porepressure]
  []
[]

# ========================================================================
# BOUNDARY CONDITIONS - Apply these functions to your boundaries
# ========================================================================
[BCs]
  # Apply horizontal displacement on bottom boundary (z=0)
  [ux_bottom]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'  # change to match your mesh boundary name
    function = ux_bc
  []


  # Apply vertical displacement on bottom boundary (z=0)
  [uz_bottom]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = uz_bc
  []




  # Apply pore pressure on bottom boundary (z=0)
  [p_bottom]
    type = FunctionDirichletBC
    variable = porepressure  # or 'pressure', depending on your variable name
    boundary = 'bottom'
    function = p_bc
  []
  #  [p_right]
  #  type = DirichletBC
  #  variable = porepressure  # or 'pressure', depending on your variable name
 #   boundary = 'right'
 #   value = 0
 # []

[]



[Functions]
  # ========================================================================
  # HORIZONTAL DISPLACEMENT ux(x,t) = U0 * phi(x) * g(t)
  # ux(x,t) = U0 * exp(-((x-x0)^2)/Lx^2) * (t/t0)*exp(-t/t0)
  # ========================================================================
  [ux_bc]
    type = ParsedFunction
    expression = 'if(t <= 0, 0, U0 * exp(-((x - x0)^2) / (Lx^2)) * (t/t0) * exp(-t/t0))'
    symbol_names = 'U0   x0   Lx   t0'
    symbol_values = '0.001  0.0  3.0  0.5'
    # U0: displacement amplitude (5*U0 from MATLAB)
    # x0: center of Gaussian patch
    # Lx: half-width of localization
    # t0: time scale of pulse (set to ~t_end/5)
  []

  # ========================================================================
  # VERTICAL DISPLACEMENT uz(x,t) = beta*U0 * phi(x) * g(t)
  # uz(x,t) = beta*U0 * exp(-((x-x0)^2)/Lx^2) * (t/t0)*exp(-t/t0)
  # ========================================================================
  [uz_bc]
    type = ParsedFunction
    expression = 'if(t <= 0, 0, U0 * exp(-((x - x0)^2) / (Lx^2)) * (t/t0) * exp(-t/t0))'
    symbol_names = 'U0   x0   Lx   t0'
    symbol_values = '0.00 0.0  3.0   0.5'
    # U0: displacement amplitude (beta*U0 = 0.5*1.0 from MATLAB)
    # x0: center of Gaussian patch
    # Lx: half-width of localization
    # t0: time scale of pulse
  []

  # ========================================================================
  # PORE PRESSURE p(x,t) = P0 * phi(x) * g(t)
  # p(x,t) = P0 * exp(-((x-x0)^2)/Lx^2) * (t/t0)*exp(-t/t0)
  # ========================================================================
  [p_bc]
    type = ParsedFunction
    expression = 'if(t <= 0, 0, P0 * exp(-((x - x0)^2) / (Lx^2)) * (t/t0) * exp(-t/t0))'
    symbol_names = 'P0   x0   Lx   t0'
    symbol_values = '0  0.0  3.0  0.5'
    # P0: pressure amplitude
    # x0: center of Gaussian patch
    # Lx: half-width of localization
    # t0: time scale of pulse
  []
[]



[AuxVariables]
  [stress_xx]
    order = FIRST
    family = MONOMIAL
  []
  [stress_yy]
    order = FIRST
    family = MONOMIAL
  []
  [flux_y]
    order = FIRST
    family = MONOMIAL
  []
  [oneroverbiotmodulus]
    order = FIRST
    family = MONOMIAL
  []
    [grad_porepressure_y]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 1
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [bulk_vel_y]
    type = PorousFlowDarcyVelocityComponent
    variable = flux_y
    component = y
    fluid_phase = 0
  []

  [grad_porepressure_y]
    type = MaterialRealAux
    variable = oneroverbiotmodulus
    property = PorousFlow_constant_biot_modulus_qp
  []

[]

[Kernels]
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_y
    component = 1
  []
  [mass0]
    type = PorousFlowFullySaturatedMassTimeDerivative
    biot_coefficient = 0.6
    coupling_type = HydroMechanical
    variable = porepressure
  []
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = porepressure
    gravity = '0 0 0'
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2.2e9
    density0 = 1
    thermal_expansion = 0
    viscosity = 1e-3
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '4.286e8 1e9'
    # bulk modulus is lambda + 2*mu/3 = 0.5 + 2*0.75/3 = 1
    fill_method = symmetric_isotropic
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.1282
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.6
    solid_bulk_compliance = 9.13043e-10
    fluid_bulk_modulus = 2.2e9
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '3.6691e-17 0 0   0 3.6691e-17 0   0 0 3.6691e-17'
  []
    # Add this material
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
[]

[Preconditioning]
  [andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = 0
  end_time = 1
  dt = 0.005
  automatic_scaling = true

[]

[Outputs]
  exodus = true
[]