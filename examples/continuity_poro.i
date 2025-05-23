# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 101
    ny = 101
    xmin = -2000
    xmax = 2000
    ymin = -2000
    ymax = 2000
    elem_type = QUAD4
  []
[]

[GlobalParams]
    displacements = 'disp_x disp_y' 
     PorousFlowDictator = dictator
[]

[Variables]
    [./disp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./porepressure]
        order = FIRST
        family = LAGRANGE
    [../]
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
    pc = 0
  []
[]

[AuxVariables]
  [./vel_x]
    order = FIRST
    family = LAGRANGE
  []
  [./accel_x]
  []
  [./vel_y]
    order = FIRST
    family = LAGRANGE
  []
  [./accel_y]
  []
  # New MMS solution aux variables
  [./mms_disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_p]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_vel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_vel_y]
    order = FIRST
    family = LAGRANGE
  [../]
  
  # Error fields for all variables
  [./error_disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./error_disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./error_p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [velocity_x]
    type = CompVarRate
    variable = vel_x
    coupled = disp_x
  []
  [velocity_y]
    type = CompVarRate
    variable = vel_y
    coupled = disp_y
  [] 
  # MMS solution display - continuous solution
  [./mms_disp_x_aux]
    type = FunctionAux
    function = displacement_x
    variable = mms_disp_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_disp_y_aux]
    type = FunctionAux
    function = displacement_y
    variable = mms_disp_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_p_aux]
    type = FunctionAux
    function = pressure
    variable = mms_p
    execute_on = 'INITIAL TIMESTEP_END'
  [../]

  
  # Error calculation aux kernels
  [./error_disp_x_aux]
    type = ParsedAux
    variable = error_disp_x
    coupled_variables = 'disp_x mms_disp_x'
    expression = 'abs(disp_x - mms_disp_x)'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./error_disp_y_aux]
    type = ParsedAux
    variable = error_disp_y
    coupled_variables = 'disp_y mms_disp_y'
    expression = 'abs(disp_y - mms_disp_y)'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./error_p_aux]
    type = ParsedAux
    variable = error_p
    coupled_variables = 'porepressure mms_p'
    expression = 'abs(porepressure - mms_p)'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]


[Kernels]
  [./stressdiv_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y'
    use_displaced_mesh = false   
  [../]
  [./stressdiv_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y'
    use_displaced_mesh = false
  [../]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
  [../]
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_x
    component = 0
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
  []

  # MMS source terms - using continuous functions
  [./source_term_x]
    type = BodyForce
    variable = disp_x
    function = source_x
  [../]
  [./source_term_y]
    type = BodyForce
    variable = disp_y
    function = source_y
  [../]
  [./source_term_p]
    type = BodyForce
    variable = porepressure
    function = source_p
  [../]
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '20e9 30e9'  # Using lambda=0.5, mu=0.75 from Mandel problem
    fill_method = symmetric_isotropic
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2850  # Using rho=1.0 from Mandel parameters
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
    porosity_zero = 0.1
    biot_coefficient = 0.6
    solid_bulk = 40e9
    constant_fluid_bulk_modulus = 2.25e9
    constant_biot_modulus = 8e9
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.0e-14 0 0   0 1.0e-14 0   0 0 1.0e-14'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0 # unimportant in this fully-saturated situation
    phase = 0
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2.25e9  # Kf from Mandel problem
    density0 = 1000  # Simplified density
    thermal_expansion = 0
    viscosity = 1.0e-3 # Simplified viscosity
  []
[]


[Functions]
  # ================ DISPLACEMENT FUNCTIONS ================
  [./displacement_x]
    type = ParsedFunction
    expression = '0.5*delta*(tanh((t - tbar)/tw) + 1)*exp(-(x + y)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  [./displacement_y]
    type = ParsedFunction
    expression = '0.5*delta*(tanh((t - tbar)/tw) + 1)*exp(-(x + y)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  # ================ PRESSURE FUNCTION ================
  [./pressure]
    type = ParsedFunction
    expression = '0.5*po*(tanh((t - tbar)/tw) + 1)*exp(-(x + y)/width^2)'
    symbol_names = 'po tbar tw width'
    symbol_values = '1 0 0.1 100.0'
  [../]

  # ================ VELOCITY FUNCTIONS ================
  [./initial_velocity_x]
    type = ParsedFunction
    expression = '0.5*delta*exp(-(x + y)/width^2)/(tw*cosh((t - tbar)/tw)^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  [./initial_velocity_y]
    type = ParsedFunction
    expression = '0.5*delta*exp(-(x + y)/width^2)/(tw*cosh((t - tbar)/tw)^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  # ================ SOURCE TERMS ================
  [./source_x]
    type = ParsedFunction
    expression = '-(0.5*alpha*po*tw^2*width^2*(tanh((t - tbar)/tw) + 1) + 1.0*delta*rho*width^4*tanh((t - tbar)/tw)/cosh((t - tbar)/tw)^2 + delta*tw^2*(1.0*lambda + 2.0*mu)*(tanh((t - tbar)/tw) + 1))*exp(-(x + y)/width^2)/(tw^2*width^4)'
    symbol_names = 'delta po tbar tw width lambda mu alpha rho'
    symbol_values = '0.1 1 0 0.1 100.0 0.5 0.75 0.6 1'
  [../]

  [./source_y]
    type = ParsedFunction
    expression = '-(0.5*alpha*po*tw^2*width^2*(tanh((t - tbar)/tw) + 1) + 1.0*delta*rho*width^4*tanh((t - tbar)/tw)/cosh((t - tbar)/tw)^2 + delta*tw^2*(1.0*lambda + 2.0*mu)*(tanh((t - tbar)/tw) + 1))*exp(-(x + y)/width^2)/(tw^2*width^4)'
    symbol_names = 'delta po tbar tw width lambda mu alpha rho'
    symbol_values = '0.1 1 0 0.1 100.0 0.5 0.75 0.6 1'
  [../]

  [./source_p]
    type = ParsedFunction
    expression = '-(1.0*M*alpha*delta*mu_f*width^2/cosh((t - tbar)/tw)^2 + 1.0*M*kappa*po*tw*(tanh((t - tbar)/tw) + 1) - 0.5*mu_f*po*width^4/cosh((t - tbar)/tw)^2)*exp(-(x + y)/width^2)/(M*mu_f*tw*width^4)'
    symbol_names = 'delta po tbar tw width alpha M kappa mu_f'
    symbol_values = '0.1 1 0 0.1 100.0 0.6 4.7059 1.5 1'
  [../]
[]
[BCs]
  # Matched value boundary conditions for all boundaries

  [./disp_x_left]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left'
    function = displacement_x
  [../]
  [./disp_y_left]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'left'
    function = displacement_y
  [../]
  [./p_left]
    type = FunctionDirichletBC
    variable = porepressure
    boundary = 'left'
    function = pressure
  [../]

  # Right boundary (x = 5)
  [./disp_x_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = displacement_x
  [../]
  [./disp_y_right]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'right'
    function = displacement_y
  [../]
  [./p_right]
    type = FunctionDirichletBC
    variable = porepressure
    boundary = 'right'
    function = pressure
  [../]

  # Bottom boundary (y = -5)
  [./disp_x_bottom]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'
    function = displacement_x
  [../]
  [./disp_y_bottom]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = displacement_y
  [../]
  [./p_bottom]
    type = FunctionDirichletBC
    variable = porepressure
    boundary = 'bottom'
    function = pressure
  [../]

  # Top boundary (y = 5)
  [./disp_x_top]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    function = displacement_x
  [../]
  [./disp_y_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = displacement_y
  [../]
  [./p_top]
    type = FunctionDirichletBC
    variable = porepressure
    boundary = 'top'
    function = pressure
  [../]
[]

[ICs]
  # Continuous domain initial conditions
  [./disp_x_ic]
    type = FunctionIC
    variable = disp_x
    function = displacement_x
  [../]
  
  [./disp_y_ic]
    type = FunctionIC
    variable = disp_y
    function = displacement_y
  [../]
  
  [./p_ic]
    type = FunctionIC
    variable = porepressure
    function = pressure
  [../]
[]


[Preconditioning]
    [./smp]
        type = SMP
        full = true

    [../]
[]

[Executioner]
  type = Transient
  dt = 0.001
  end_time = 0.5
  automatic_scaling = true
  [./TimeIntegrator]
    type = CentralDifference
  [../]
[]

[Outputs]
    exodus = true
    time_step_interval = 5
[]

