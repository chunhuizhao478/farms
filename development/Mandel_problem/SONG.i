[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 400
  ny = 200
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
# BOUNDARY CONDITIONS
# ========================================================================
[BCs]
  # BOTTOM BOUNDARY (y=0)
  # ux = 0.5 * heaviside(-x) * heaviside(t)
  [ux_bottom]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'
    function = ux_bottom_func
  []

  # flux = -3.669e-4 * porepressure / 1e-3 = -0.3669 * porepressure
  [flux_bottom]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = 'bottom'
    pt_vals = '0 1'
    multipliers = '0 -0.3669'
    fluid_phase = 0
    flux_function = 1
    use_mobility = false
    use_relperm = false
  []


[]

[Functions]
  # ux at bottom: 0.5 * heaviside(-x) * heaviside(t)
  [ux_bottom_func]
    type = ParsedFunction
    expression = '0.5 * if(x < 0, 1, 0) * if(t > 0, 1, 0)'
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
    bulk_modulus = 0.8
    density0 = 1
    thermal_expansion = 0
    viscosity = 1
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '0.428571 1'
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
    solid_bulk_compliance = 0.91304
    fluid_bulk_modulus = 0.8
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '3.6691e-4 0 0   0 3.6691e-4 0   0 0 3.6691e-4'
  []
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
  end_time = 10
  dt = 0.001
  automatic_scaling = true
[]

[Outputs]
  exodus = true
[]