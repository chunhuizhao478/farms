

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 50
  xmin = -1
  xmax = 1
  ymin = 0
  ymax = 10
  elem_type = QUAD9
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
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
  [pc]
    type = PorousFlowCapillaryPressureConst
  []
[]

[Variables]
  [disp_x]
  order = SECOND
  family = LAGRANGE
  []
  [disp_y]
  order = SECOND
  family = LAGRANGE
  []
  [porepressure]
  []
[]

[BCs]
  [confinex]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left right'
  []
  [basefixed]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = bottom
  []
  [topdrained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = top
  []
  [topload]
    type = NeumannBC
    variable = disp_y
    value = -1
    boundary = top
  []
[]

[Kernels]
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
   # implicit = false
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
   # implicit = false
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_x
    component = 0
    #implicit = false
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_y
    component = 1
   # implicit = false
  []
  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = porepressure
    fluid_component = 0
  []
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = porepressure
  []
  [flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
    #implicit = false
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 8
    density0 = 1
    thermal_expansion = 0
    viscosity = 0.96
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '2 3'
    # bulk modulus is lambda + 2*mu/3 = 2 + 2*3/3 = 4
    fill_method = symmetric_isotropic
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [ppss]
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityHMBiotModulus
    porosity_zero = 0
    biot_coefficient = 1
    solid_bulk = 4
    constant_fluid_bulk_modulus = 8
    constant_biot_modulus = 1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.5 0 0   0 1.5 0   0 0 1.5'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0 # unimportant in this fully-saturated situation
    phase = 0
  []
[]


[Preconditioning]
  [andy]
    type = SMP
    full = true
 []
[]

[Executioner]
  type = Transient
  solve_type = Newton
   dt = 0.005
  start_time = 0
  end_time = 10
[]
#[Executioner]
#    type = Transient
 #   dt = 0.005
 #   start_time= 0
 #   end_time = 10
 #   [TimeIntegrator]
  #  type = CrankNicolson
  #  []
 #[]

[Outputs]
  execute_on = 'timestep_end'
  [csv]
    type = CSV
  []
  exodus = true
[]