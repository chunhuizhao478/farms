[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 50
    xmin = -4
    xmax = 4
    ymin = 0
    ymax = 4
    #elem_type = QUAD9
  []
  [subdivide_block]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 1
    bottom_left = '-4 0 0'
    top_right = '4 0.4 0'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator
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
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [porepressure]
    order = FIRST
    family = LAGRANGE
  []
[]

[BCs]
  [roller_xmin]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'left right'
  []
  [roller_ymin]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'bottom top'
  []
  [bottom_flux]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = left
    function = flux_fcn
  []
  [top_pressure]
    type = DirichletBC
    variable = porepressure
    boundary = 'top'
    value = 0
  []
[]

[Functions]
  [flux_fcn]
    type = ParsedFunction
    value = 'if(x <= 0.4, 1e-6, 1e-8)'
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
    biot_coefficient = 0.8145
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.8145
    variable = disp_y
    component = 1
  []
  [mass0]
    type = PorousFlowFullySaturatedMassTimeDerivative
    biot_coefficient = 0.8145
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
    bulk_modulus = 0.8145
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
    C_ijkl = '30e9 20e9'
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
    porosity = 0.0274
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.8145
    solid_bulk_compliance = 2.5e-11
    fluid_bulk_modulus = 2e9
  []
  
  # Higher permeability for block 1 (y < 0.4)
  [permeability_lower]
    type = PorousFlowPermeabilityConst
    permeability = '0.025e-9 0 0 0 0.025e-9 0 0 0 0.025e-9'
    block = 1
  []
  
  # Original permeability for block 0 (y >= 0.4)
  [permeability_upper]
    type = PorousFlowPermeabilityConst
    permeability = '0.025e-11 0 0 0 0.025e-11 0 0 0 0.025e-11'
    block = 0
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
  solve_type = PJFNK
  start_time = 0
  end_time = 0.01
  dt = 1e-5
   automatic_scaling = true


 line_search = 'bt'  # Use backtracking line search

[]

[Outputs]
  exodus = true
[]