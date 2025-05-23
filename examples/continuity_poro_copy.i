# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = -500
    xmax = 500
    ymin = -500
    ymax = 500
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
    variable = disp_y
    component = 1
  []
  [mass0]
    type = PorousFlowFullySaturatedMassTimeDerivative
    biot_coefficient = 0.6
    coupling_type = HydroMechanical
    variable = porepressure
    multiply_by_density = false
  []
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = porepressure
    gravity = '0 0 0'
    multiply_by_density = false
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
  [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 0.5
        shear_modulus = 0.75
        use_displaced_mesh = false
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
    prop_values = 1.0  # Using rho=1.0 from Mandel parameters
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
    type = PorousFlowPorosityConst
    porosity = 0.1  # Using porosity=0.1 from Mandel problem
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.6  # Alpha value from both examples
    solid_bulk_compliance = 1.0  # 1/K from Mandel problem (K=1.0)
    fluid_bulk_modulus = 8.0  # Kf=8.0 from Mandel problem
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.5 0 0   0 1.5 0   0 0 1.5'  # k=1.5 from Mandel problem
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 8.0  # Kf from Mandel problem
    density0 = 1.0  # Simplified density
    thermal_expansion = 0
    viscosity = 1.0  # Simplified viscosity
  []
[]


[Functions]
  # ================ DISPLACEMENT FUNCTIONS ================
  [./displacement_x]
    type = ParsedFunction
    expression = 'delta*t^3*exp(-(x^2 + y^2)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  [./displacement_y]
    type = ParsedFunction
    expression = 'delta*t^3*exp(-(x^2 + y^2)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  # ================ PRESSURE FUNCTION ================
  [./pressure]
    type = ParsedFunction
    expression = 'po*t^3*exp(-(x^2 + y^2)/width^2)'
    symbol_names = 'po tbar tw width'
    symbol_values = '1.0e3 0 0.1 100.0'
  [../]

  # ================ VELOCITY FUNCTIONS ================
  [./initial_velocity_x]
    type = ParsedFunction
    expression = '3*delta*t^2*exp(-(x^2 + y^2)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  [./initial_velocity_y]
    type = ParsedFunction
    expression = '3*delta*t^2*exp(-(x^2 + y^2)/width^2)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '0.1 0 0.1 100.0'
  [../]

  # ================ SOURCE TERMS ================
  [./source_x]
    type = ParsedFunction
    expression = '2*t*(-delta*mu*t^2*(-1.0*width^2 + 2.0*y*(x + y)) + 3*delta*rho*width^4 - 2*delta*t^2*x*(lambda*y + x*(lambda + 2*mu)) + t^2*width^2*(-alpha*po*x + delta*(lambda + 2*mu)))*exp(-(x^2 + y^2)/width^2)/width^4'
    symbol_names = 'delta po tbar tw width lambda mu alpha rho'
    symbol_values = '0.1 1.0e3 0 0.1 100.0 0.5 0.75 0.6 1'
  [../]

  [./source_y]
    type = ParsedFunction
    expression = '2*t*(-delta*mu*t^2*(-1.0*width^2 + 2.0*x*(x + y)) + 3*delta*rho*width^4 - 2*delta*t^2*y*(lambda*x + y*(lambda + 2*mu)) + t^2*width^2*(-alpha*po*y + delta*(lambda + 2*mu)))*exp(-(x^2 + y^2)/width^2)/width^4'
    symbol_names = 'delta po tbar tw width lambda mu alpha rho'
    symbol_values = '0.1 1.0e3 0 0.1 100.0 0.5 0.75 0.6 1'
  [../]

  [./source_p]
    type = ParsedFunction
    expression = 't^2*(-6*M*alpha*delta*mu_f*width^2*(x + y) + 4*M*kappa*po*t*(width^2 - x^2 - y^2) + 3*mu_f*po*width^4)*exp(-(x^2 + y^2)/width^2)/(M*mu_f*width^4)'
    symbol_names = 'delta po tbar tw width alpha M kappa mu_f'
    symbol_values = '0.1 1.0e3 0 0.1 100.0 0.6 4.7059 1.5 1'
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
  dt = 0.0025
  end_time = 0.5
  #automatic_scaling = true
  [./TimeIntegrator]
        type = NewmarkBeta
  [../]
[]

[Outputs]
    exodus = true
    time_step_interval = 2
[]

