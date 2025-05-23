# Verification of Benchmark Problem TPV05-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

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
  [./new_block]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'y>0'
    block_id = 1
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
    split_interface = true
    add_interface_on_two_sides = true
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
    # MMS solution display - Block 0
  [./mms_disp_x_block0]
    type = FunctionAux
    function = displacement_x_block0
    variable = mms_disp_x
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_disp_y_block0]
    type = FunctionAux
    function = displacement_y_block0
    variable = mms_disp_y
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_p_block0]
    type = FunctionAux
    function = pressure_block0
    variable = mms_p
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_vel_x_block0]
    type = FunctionAux
    function = velocity_x_block0
    variable = mms_vel_x
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_vel_y_block0]
    type = FunctionAux
    function = velocity_y_block0
    variable = mms_vel_y
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_x_block0]
    type = FunctionAux
    function = darcy_velocity_x_block0
    variable = mms_fluid_vel_x
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_y_block0]
    type = FunctionAux
    function = darcy_velocity_y_block0
    variable = mms_fluid_vel_y
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  # MMS solution display - Block 1
  [./mms_disp_x_block1]
    type = FunctionAux
    function = displacement_x_block1
    variable = mms_disp_x
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_disp_y_block1]
    type = FunctionAux
    function = displacement_y_block1
    variable = mms_disp_y
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_p_block1]
    type = FunctionAux
    function = pressure_block1
    variable = mms_p
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_vel_x_block1]
    type = FunctionAux
    function = velocity_x_block1
    variable = mms_vel_x
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_vel_y_block1]
    type = FunctionAux
    function = velocity_y_block1
    variable = mms_vel_y
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_x_block1]
    type = FunctionAux
    function = darcy_velocity_x_block1
    variable = mms_fluid_vel_x
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./mms_fluid_vel_y_block1]
    type = FunctionAux
    function = darcy_velocity_y_block1
    variable = mms_fluid_vel_y
    block = 1
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
      # MMS source terms
    [./source_term_x0]
      type = BodyForce
      variable = disp_x
      function = source_x_block0
      block = 0
    [../]
    [./source_term_y0]
      type = BodyForce
      variable = disp_y
      function = source_y_block0
      block = 0
    [../]
    [./source_term_p0]
      type = BodyForce
      variable = p
      function = source_p_block0
      block = 0
    [../]
    [./source_term_vfx0]
      type = BodyForce
      variable = fluid_vel_x
      function = source_vfx_block0
      block = 0
    [../]
    [./source_term_vfy0]
      type = BodyForce
      variable = fluid_vel_y
      function = source_vfy_block0
      block = 0
    [../]

    [./source_term_x1]
      type = BodyForce
      variable = disp_x
      function = source_x_block1
      block = 1
    [../]
    [./source_term_y1]
      type = BodyForce
      variable = disp_y
      function = source_y_block1
      block = 1
    [../]
    [./source_term_p1]
      type = BodyForce
      variable = p
      function = source_p_block1
      block = 1
    [../]
    [./source_term_vfx1]
      type = BodyForce
      variable = fluid_vel_x
      function = source_vfx_block1
      block = 1
    [../]
    [./source_term_vfy1]
      type = BodyForce
      variable = fluid_vel_y
      function = source_vfy_block1
      block = 1
    [../]
    
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 0.5       # Updated to match symbol_values
    shear_modulus = 0.75  # Updated to match 'mu' in symbol_values
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
    prop_values = 2.5  # Updated to match 'rho' in symbol_values
  []
  [./rhof]
    type = GenericConstantMaterial
    prop_names = rhof
    prop_values = 1  # Already matches 'rhof' in symbol_values
  [../]
  [./porosity]
    type = GenericConstantMaterial
    prop_names = porosity
    prop_values = 0.1   # Keep as is (not directly used in MMS functions)
  [../]
  [./hydconductivity]
    type = GenericConstantMaterial
    prop_names = hydconductivity
    prop_values = 1.5  # Updated to match 'k' in symbol_values
  [../]
  [./hydconductivity_layer]
    type = GenericConstantMaterial
    prop_names = hydconductivity_layer
    prop_values = 1.5   # Updated to match 'k' in symbol_values
  [../]
  [./biotcoeff]
    type = GenericConstantMaterial
    prop_names = biot_coefficient
    prop_values = 0.6  # Updated to match 'alpha' in symbol_values
  [../]
  [./biotmodulus]
    type = GenericConstantMaterial
    prop_names = biot_modulus
    prop_values = 4.7059   # Updated to match 'M' in symbol_values
  [../]
  [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 2
  [../]
  [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
  [../]
[]

[Functions]
  # ================ BLOCK 0 (LOWER DOMAIN) DISPLACEMENT FUNCTIONS ================
  [./displacement_x_block0]
    type = ParsedFunction
    expression = '-dx*t^3*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./displacement_y_block0]
    type = ParsedFunction
    expression = '-dy*t^3*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) DISPLACEMENT FUNCTIONS ================
  [./displacement_x_block1]
    type = ParsedFunction
    expression = 'dx*t^3*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./displacement_y_block1]
    type = ParsedFunction
    expression = 'dy*t^3*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) PRESSURE FUNCTION ================
  [./pressure_block0]
    type = ParsedFunction
    expression = '-p*t^3*x*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'p tb tw R'
    symbol_values = '1.0e2 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) PRESSURE FUNCTION ================
  [./pressure_block1]
    type = ParsedFunction
    expression = 'p*t^3*x*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'p tb tw R'
    symbol_values = '1.0e2 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) VELOCITY FUNCTIONS ================
  [./velocity_x_block0]
    type = ParsedFunction
    expression = '-3*dx*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./velocity_y_block0]
    type = ParsedFunction
    expression = '-3*dy*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) VELOCITY FUNCTIONS ================
  [./velocity_x_block1]
    type = ParsedFunction
    expression = '3*dx*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./velocity_y_block1]
    type = ParsedFunction
    expression = '3*dy*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x_block0]
    type = ParsedFunction
    expression = 'k*(R^2*mf*t*(6*dx*rf + p*t^2) - 3*k*rp*(R^2*(2*dx*rf + p*t^2) - 2*p*t^2*x^2) - 2*mf*p*t^3*x^2)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 1.0e2 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block0]
    type = ParsedFunction
    expression = '2*k*(3*R^2*dy*mf*rf*t - 3*k*rp*(R^2*dy*rf - p*t^2*x*y) - mf*p*t^3*x*y)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 1.0e2 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x_block1]
    type = ParsedFunction
    expression = 'k*(-R^2*mf*t*(6*dx*rf + p*t^2) + 3*k*rp*(R^2*(2*dx*rf + p*t^2) - 2*p*t^2*x^2) + 2*mf*p*t^3*x^2)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 1.0e2 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block1]
    type = ParsedFunction
    expression = '2*k*(-3*R^2*dy*mf*rf*t + 3*k*rp*(R^2*dy*rf - p*t^2*x*y) + mf*p*t^3*x*y)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 1.0e2 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block0]
    type = ParsedFunction
    expression = '(-6*R^4*dx*mf^2*r*t + 3*R^2*k*rf*(R^2*mf*(2*dx*rf + p*t^2) - 2*k*p*rp*t*(R^2 - 2*x^2) - 2*mf*p*t^2*x^2) - R^2*mf^2*t^3*(a*p*(R^2 - 2*x^2) + 4*dx*mu) + 8*dx*mf^2*mu*t^3*x^2 + 2*mf^2*t^3*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block0]
    type = ParsedFunction
    expression = '2*(-3*R^4*dy*mf^2*r*t + 3*R^2*k*rf*(R^2*dy*mf*rf + 2*k*p*rp*t*x*y - mf*p*t^2*x*y) + R^2*mf^2*t^3*(a*p*x*y - 2*dy*mu) + 4*dy*mf^2*mu*t^3*y^2 + mf^2*t^3*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp((-x^2 - y^2)/R^2)/(R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block1]
    type = ParsedFunction
    expression = '(6*R^4*dx*mf^2*r*t - 3*R^2*k*rf*(R^2*mf*(2*dx*rf + p*t^2) - 2*k*p*rp*t*(R^2 - 2*x^2) - 2*mf*p*t^2*x^2) + R^2*mf^2*t^3*(a*p*(R^2 - 2*x^2) + 4*dx*mu) - 8*dx*mf^2*mu*t^3*x^2 - 2*mf^2*t^3*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block1]
    type = ParsedFunction
    expression = '2*(3*R^4*dy*mf^2*r*t - 3*R^2*k*rf*(R^2*dy*mf*rf + 2*k*p*rp*t*x*y - mf*p*t^2*x*y) + R^2*mf^2*t^3*(-a*p*x*y + 2*dy*mu) - 4*dy*mf^2*mu*t^3*y^2 - mf^2*t^3*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block0]
    type = ParsedFunction
    expression = '(6*M*R^2*a*mf^2*t^2*(dx*x + dy*y) + 2*M*k*(-R^2*mf*t*(6*dy*rf*y + p*t^2*x) + 3*k*rp*(R^2*(2*dy*rf*y + p*t^2*x) - 2*p*t^2*x*y^2) + 2*mf*p*t^3*x*y^2 + x*(-3*R^2*mf*t*(2*dx*rf + p*t^2) + 3*k*rp*(R^2*(2*dx*rf + 3*p*t^2) - 2*p*t^2*x^2) + 2*mf*p*t^3*x^2)) - 3*R^4*mf^2*p*t^2*x)*exp(-(x^2 + y^2)/R^2)/(M*R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block1]
    type = ParsedFunction
    expression = '(-6*M*R^2*a*mf^2*t^2*(dx*x + dy*y) + 2*M*k*(R^2*mf*t*(6*dy*rf*y + p*t^2*x) - 3*k*rp*(R^2*(2*dy*rf*y + p*t^2*x) - 2*p*t^2*x*y^2) - 2*mf*p*t^3*x*y^2 - x*(-3*R^2*mf*t*(2*dx*rf + p*t^2) + 3*k*rp*(R^2*(2*dx*rf + 3*p*t^2) - 2*p*t^2*x^2) + 2*mf*p*t^3*x^2)) + 3*R^4*mf^2*p*t^2*x)*exp(-(x^2 + y^2)/R^2)/(M*R^4*mf^2)'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block0]
    type = ParsedFunction
    expression = '6*k^2*p*rp^2*t*(-R^2 + 2*x^2)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block0]
    type = ParsedFunction
    expression = '12*k^2*p*rp^2*t*x*y*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.001 1.0e2 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block1]
    type = ParsedFunction
    expression = '6*k^2*p*rp^2*t*(R^2 - 2*x^2)*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block1]
    type = ParsedFunction
    expression = '-12*k^2*p*rp^2*t*x*y*exp(-(x^2 + y^2)/R^2)/(R^2*mf^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 1.0e2 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
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
  [./p_left]
    type = DirichletBC
    variable = p
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
  [./p_right]
    type = DirichletBC
    variable = p
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
  [./p_top]
    type = DirichletBC
    variable = p
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
  [./p_bottom]
    type = DirichletBC
    variable = p
    boundary = bottom
    value = 0
  [../]
  
  # Fluid velocity BCs (zero flux at domain boundaries)
  [./fluid_vel_x_all]
    type = DirichletBC
    variable = fluid_vel_x
    boundary = 'left right top bottom'
    value = 0
  [../]
  [./fluid_vel_y_all]
    type = DirichletBC
    variable = fluid_vel_y
    boundary = 'left right top bottom'
    value = 0
  [../]

  
  # Fluid velocity BCs (zero flux at domain boundaries)
  [ssssss]
    type = MatchedValueBC
    variable = disp_x
    boundary = 'Block1_Block0'
    v = mms_disp_x
  []
  [./ssssacas]
    type = MatchedValueBC
    variable = disp_y
    boundary = 'Block0_Block1'
    v = mms_disp_y
  [../]
  [sssssfevs]
    type = MatchedValueBC
    variable = disp_x
    boundary = 'Block0_Block1'
    v = mms_disp_x
  []
  [./dvf]
    type = MatchedValueBC
    variable = disp_y
    boundary = 'Block1_Block0'
    v = mms_disp_y
  [../]


[]


[ICs]
  # Block 0 (lower domain) initial conditions
  [./disp_x_ic_block0]
    type = FunctionIC
    variable = disp_x
    function = displacement_x_block0
    block = 0
  [../]
  
  [./disp_y_ic_block0]
    type = FunctionIC
    variable = disp_y
    function = displacement_y_block0
    block = 0
  [../]
  
  [./p_ic_block0]
    type = FunctionIC
    variable = p
    function = pressure_block0
    block = 0
  [../]
  
  [./fluid_vel_x_ic_block0]
    type = FunctionIC
    variable = fluid_vel_x
    function = darcy_velocity_x_block0
    block = 0
  [../]
  
  [./fluid_vel_y_ic_block0]
    type = FunctionIC
    variable = fluid_vel_y
    function = darcy_velocity_y_block0
    block = 0
  [../]
  
  # Block 1 (upper domain) initial conditions
  [./disp_x_ic_block1]
    type = FunctionIC
    variable = disp_x
    function = displacement_x_block1
    block = 1
  [../]
  
  [./disp_y_ic_block1]
    type = FunctionIC
    variable = disp_y
    function = displacement_y_block1
    block = 1
  [../]
  
  [./p_ic_block1]
    type = FunctionIC
    variable = p
    function = pressure_block1
    block = 1
  [../]
  
  [./fluid_vel_x_ic_block1]
    type = FunctionIC
    variable = fluid_vel_x
    function = darcy_velocity_x_block1
    block = 1
  [../]
  
  [./fluid_vel_y_ic_block1]
    type = FunctionIC
    variable = fluid_vel_y
    function = darcy_velocity_y_block1
    block = 1
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
   # type = NewmarkBeta
    type = CentralDifference
  [../]
[]

[Outputs]
    exodus = true
    time_step_interval = 5
[]

