# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = -1000
    xmax = 1000
    ymin = -1000
    ymax = 1000
    elem_type = QUAD4
  []
[]


[GlobalParams]
    displacements = 'disp_x disp_y' 
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
    [./fluid_vel_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./fluid_vel_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./p]
        order = FIRST
        family = LAGRANGE
    [../]
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
  [./mms_fluid_vel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_fluid_vel_y]
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
  [./error_fluid_vel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./error_fluid_vel_y]
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
  [./mms_vel_x_aux]
    type = FunctionAux
    function = velocity_x
    variable = mms_vel_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_vel_y_aux]
    type = FunctionAux
    function = velocity_y
    variable = mms_vel_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_x_aux]
    type = FunctionAux
    function = darcy_velocity_x
    variable = mms_fluid_vel_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_y_aux]
    type = FunctionAux
    function = darcy_velocity_y
    variable = mms_fluid_vel_y
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
    coupled_variables = 'p mms_p'
    expression = 'abs(p - mms_p)'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./error_fluid_vel_x_aux]
    type = ParsedAux
    variable = error_fluid_vel_x
    coupled_variables = 'fluid_vel_x mms_fluid_vel_x'
    expression = 'abs(fluid_vel_x - mms_fluid_vel_x)'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./error_fluid_vel_y_aux]
    type = ParsedAux
    variable = error_fluid_vel_y
    coupled_variables = 'fluid_vel_y mms_fluid_vel_y'
    expression = 'abs(fluid_vel_y - mms_fluid_vel_y)'
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
  [./skeletoninertia_x]
    type = InertialForce
    variable = disp_x
    use_displaced_mesh = false
  [../]
  [./skeletoninertia_y]
    type = InertialForce
    variable = disp_y
    use_displaced_mesh = false
  [../]
  [./porefluidIFcoupling_x]
    type = CoupledFluidInertialForce
    variable = disp_x
    fluid_vel = fluid_vel_x
    use_displaced_mesh = false
  [../]
  [./porefluidIFcoupling_y]
    type = CoupledFluidInertialForce
    variable = disp_y
    fluid_vel = fluid_vel_y
    use_displaced_mesh = false
  [../]
  [./darcyflow_x]
    type = DynamicDarcyFlow2
    variable = fluid_vel_x
    skeleton_acceleration = disp_x
  [../]
  [./darcyflow_y]
    type = DynamicDarcyFlow2
    variable = fluid_vel_y
    skeleton_acceleration = disp_y
  [../]
  [./poromechskeletoncoupling_x]
    type = PoroMechanicsCoupling
    variable = disp_x
    porepressure = p
    component = 0
  [../]
  [./poromechskeletoncoupling_y]
    type = PoroMechanicsCoupling
    variable = disp_y
    porepressure = p
    component = 1
  [../]
  [./poromechfluidcoupling_x]
    type = PoroMechanicsCoupling2
    variable = fluid_vel_x
    porepressure = p
    component = 0
  [../]
  [./poromechfluidcoupling_y]
    type = PoroMechanicsCoupling2
    variable = fluid_vel_y
    porepressure = p
    component = 1
  [../]
  [./massconservationskeleton]
    type = INSmassSolid
    variable = p
    displacements = 'disp_x disp_y'
  [../]
  [./massconservationpressure]
    type = FluidStorage
    variable = p
  [../]
  [./massconservationfluid]
    type = INSmassFluid
    variable = p
    u = fluid_vel_x
    v = fluid_vel_y
    pressure = p
  [../]
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
    variable = p
    function = source_p
  [../]
  [./source_term_vfx]
    type = BodyForce
    variable = fluid_vel_x
    function = source_vfx
  [../]
  [./source_term_vfy]
    type = BodyForce
    variable = fluid_vel_y
    function = source_vfy
  [../]
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 20e9       # Updated to match symbol_values
    shear_modulus = 30e9  # Updated to match 'mu' in symbol_values
    use_displaced_mesh = false
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [Strain]
    type = ComputeSmallStrain
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2670  # Updated to match 'rho' in symbol_values
  []
  [./rhof]
    type = GenericConstantMaterial
    prop_names = rhof
    prop_values = 1000  # Already matches 'rhof' in symbol_values
  [../]
  [./porosity]
    type = GenericConstantMaterial
    prop_names = porosity
    prop_values = 0.1   # Keep as is (not directly used in MMS functions)
  [../]
  [./hydconductivity]
    type = GenericConstantMaterial
    prop_names = hydconductivity
    prop_values = 1e-12  # Updated to match 'k' in symbol_values
  [../]
  [./hydconductivity_layer]
    type = GenericConstantMaterial
    prop_names = hydconductivity_layer
    prop_values = 1e-12  # Updated to match 'k' in symbol_values
  [../]
  [./biotcoeff]
    type = GenericConstantMaterial
    prop_names = biot_coefficient
    prop_values = 0.6  # Updated to match 'alpha' in symbol_values
  [../]
  [./biotmodulus]
    type = GenericConstantMaterial
    prop_names = biot_modulus
    prop_values = 1e7   # Updated to match 'M' in symbol_values
  [../]
  [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 2.5
  [../]
  [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
  [../]
[]



[Functions]
  # ================ MODEL PARAMETERS ================
  [./model_params]
    type = ParsedFunction
    value = '0' # Dummy function to store parameters
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

  # ================ DISPLACEMENT FUNCTIONS ================
  [./displacement_x]
    type = ParsedFunction
    expression = 'deltax*exp((t - tbar)/tw - (x^2 + y^2)/(2*R^2))'
    symbol_names = 'deltax tbar tw R epsilon'
    symbol_values = '1.0 0.15 0.1 100.0 0.01'
  [../]

  [./displacement_y]
    type = ParsedFunction
    expression = 'deltay*exp((t - tbar)/tw - (x^2 + y^2)/(2*R^2))'
    symbol_names = 'deltay tbar tw R epsilon'
    symbol_values = '0.01 0.15 0.1 100.0 0.01'
  [../]

  # ================ VELOCITY FUNCTIONS ================
  [./velocity_x]
    type = ParsedFunction
    expression = 'deltax*exp((t - tbar)/tw - (x^2 + y^2)/(2*R^2))/tw'
    symbol_names = 'deltax tbar tw R epsilon'
    symbol_values = '1.0 0.15 0.1 100.0 0.01'
  [../]

  [./velocity_y]
    type = ParsedFunction
    expression = 'deltay*exp((t - tbar)/tw - (x^2 + y^2)/(2*R^2))/tw'
    symbol_names = 'deltay tbar tw R epsilon'
    symbol_values = '0.01 0.15 0.1 100.0 0.01'
  [../]

  # ================ PRESSURE FUNCTIONS ================
  [./pressure]
    type = ParsedFunction
    expression = 'p0*x*exp((t - tbar)/tw - (x^2 + y^2)/(2*R^2))'
    symbol_names = 'p0 tbar tw R'
    symbol_values = '1e2 0.25 0.1 100.0'
  [../]

  # ================ DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x]
    type = ParsedFunction
    expression = '-k*(R^2*deltax*muf*rhof*tw - R^2*muf*p0*tw^3 - k*rhop*(R^2*deltax*rhof - R^2*p0*tw^2 + p0*tw^2*x^2) + muf*p0*tw^3*x^2)*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/(R^2*muf^2*tw^3)'
    symbol_names = 'k muf p0 tbar tw R rhof rhop deltax deltay epsilon'
    symbol_values = '1e-15 1e-3 1e2 0.15 0.1 100.0 1000.0 25000.0 1.0 0.01 0.01'
  [../]

  [./darcy_velocity_y]
    type = ParsedFunction
    expression = '-k*(R^2*deltay*muf*rhof*tw - k*rhop*(R^2*deltay*rhof + p0*tw^2*x*y) + muf*p0*tw^3*x*y)*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/(R^2*muf^2*tw^3)'
    symbol_names = 'k muf p0 tbar tw R rhof rhop deltax deltay epsilon'
    symbol_values = '1e-15 1e-3 1e2 0.15 0.1 100.0 1000.0 25000.0 1.0 0.01 0.01'
  [../]

  # ================ STRESS TENSOR COMPONENTS ================
  [./stress_xx]
    type = ParsedFunction
    expression = '(-R^2*alpha*p0*x - 2*deltax*mu*x - lambda*(deltax*x + deltay*y))*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/R^2'
    symbol_names = 'mu lambda alpha p0 deltax deltay tbar tw R epsilon'
    symbol_values = '30e9 20e9 0.6 1e2 1.0 0.01 0.15 0.1 100.0 0.01'
  [../]

  [./stress_yy]
    type = ParsedFunction
    expression = '(-R^2*alpha*p0*x - 2*deltay*mu*y - lambda*(deltax*x + deltay*y))*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/R^2'
    symbol_names = 'mu lambda alpha p0 deltax deltay tbar tw R epsilon'
    symbol_values = '30e9 20e9 0.6 1e2 1.0 0.01 0.15 0.1 100.0 0.01'
  [../]

  [./stress_xy]
    type = ParsedFunction
    expression = '-1.0*mu*(deltax*y + deltay*x)*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/R^2'
    symbol_names = 'mu deltax deltay tbar tw R epsilon'
    symbol_values = '30e9 1.0 0.01 0.15 0.1 100.0 0.01'
  [../]

  # ================ SOURCE TERMS ================
  [./source_x]
    type = ParsedFunction
    expression = '(2*R^4*alpha*muf^2*p0*tw^4 - R^4*deltax*muf^2*rho*tw^2 + R^2*k*rhof*(R^2*deltax*muf*rhof*tw - R^2*muf*p0*tw^3 - k*rhop*(R^2*deltax*rhof - R^2*p0*tw^2 + p0*tw^2*x^2) + muf*p0*tw^3*x^2) + 2*R^2*muf^2*tw^4*(-alpha*p0*x^2 + deltax*mu) - 2*deltax*mu*muf^2*tw^4*x^2 - muf^2*tw^4*(lambda*(-R^2*deltax + x*(deltax*x + deltay*y)) + 1.0*mu*(-R^2*deltax + y*(deltax*y + deltay*x))))*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/(R^4*muf^2*tw^4)'
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

  [./source_y]
    type = ParsedFunction
    expression = '(-R^4*deltay*muf^2*rho*tw^2 + R^2*k*rhof*(R^2*deltay*muf*rhof*tw - k*rhop*(R^2*deltay*rhof + p0*tw^2*x*y) + muf*p0*tw^3*x*y) + 2*R^2*muf^2*tw^4*(-alpha*p0*x*y + deltay*mu) - 2*deltay*mu*muf^2*tw^4*y^2 - muf^2*tw^4*(lambda*(-R^2*deltay + y*(deltax*x + deltay*y)) + 1.0*mu*(-R^2*deltay + x*(deltax*y + deltay*x))))*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/(R^4*muf^2*tw^4)'
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

  [./source_p]
    type = ParsedFunction
    expression = '(M*R^2*alpha*muf^2*tw^2*(deltax*x + deltay*y) + M*k*(-R^2*deltay*muf*rhof*tw*y + R^2*muf*p0*tw^3*x + k*rhop*(R^2*deltay*rhof*y - R^2*p0*tw^2*x + p0*tw^2*x*y^2) - muf*p0*tw^3*x*y^2 + x*(-R^2*deltax*muf*rhof*tw + 3*R^2*muf*p0*tw^3 + k*rhop*(R^2*deltax*rhof - 3*R^2*p0*tw^2 + p0*tw^2*x^2) - muf*p0*tw^3*x^2)) - R^4*muf^2*p0*tw^2*x)*exp((2*R^2*(t - tbar) - tw*(x^2 + y^2))/(2*R^2*tw))/(M*R^4*muf^2*tw^3)'
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

  # ================ DARCY VELOCITY SOURCE TERMS ================
  [./source_vfx]
    type = ParsedFunction
    expression = 'k^2*rhop^2*(R^2*deltax*rhof - R^2*p0*tw^2 + p0*tw^2*x^2)*exp(t/tw - tbar/tw - x^2/(2*R^2) - y^2/(2*R^2))/(R^2*muf^2*tw^4)'
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

  [./source_vfy]
    type = ParsedFunction
    expression = 'k^2*rhop^2*(R^2*deltay*rhof + p0*tw^2*x*y)*exp(t/tw - tbar/tw - x^2/(2*R^2) - y^2/(2*R^2))/(R^2*muf^2*tw^4)'
    symbol_names = 'deltax deltay tbar tw R p0 mu lambda alpha rho rhof rhop k muf M gamma epsilon'
    symbol_values = '1.0 0.01 0.15 0.1 100.0 1e2 30e9 20e9 0.6 2670.0 1000.0 25000.0 1e-15 1e-3 1e7 0.0 0.01'
  [../]

[]


[BCs]
  # Dirichlet boundary conditions for outer boundaries
  # Left boundary 
  [./disp_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./disp_y_left]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  [../]
  
  # Right boundary
  [./disp_x_right]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  [../]
  [./disp_y_right]
    type = DirichletBC
    variable = disp_y
    boundary = right
    value = 0
  [../]

  # Top boundary
  [./disp_x_top]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
  [../]
  [./disp_y_top]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0
  [../]
  
  # Bottom boundary
  [./disp_x_bottom]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  [../]
  [./disp_y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
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
    variable = p
    function = pressure
  [../]
  
  [./fluid_vel_x_ic]
    type = FunctionIC
    variable = fluid_vel_x
    function = darcy_velocity_x
  [../]
  
  [./fluid_vel_y_ic]
    type = FunctionIC
    variable = fluid_vel_y
    function = darcy_velocity_y
  [../]
[]


[Preconditioning]
    [./smp]
        type = SMP
        full = true
        petsc_options = '-snes_converged_reason'
        petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
        petsc_options_value = 'lu superlu_dist 100'
    [../]
[]


[Executioner]
  type = Transient
  dt = 0.0001
  end_time = 0.3
  automatic_scaling = true
  [./TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  [../]
[]

[Outputs]
    exodus = true
    time_step_interval = 10
[]

