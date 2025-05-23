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

    ##damping ratio 
    q = 0

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
    #restoration force (tag after solve)
    [./resid_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    #restoration force for damping (tag after solve)
    [./resid_damp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damp_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    #restoration force for pressure (tag after solve)
    [./resid_press_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_press_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plusminus_scaled_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_scaled_y]
        order = FIRST
        family = LAGRANGE
    []
    #traction
    [traction_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #sliprate
    [sliprate_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [sliprate_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #slip
    [slip_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [slip_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #statevar
    [statevar]
        order = CONSTANT
        family = MONOMIAL
    []
    #velocity
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
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
    #obtain system residuals by tagging
    [restore_x]
        type = TagVectorAux
        vector_tag = 'restore_tag'
        v = 'disp_x'
        variable = 'resid_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_y]
        type = TagVectorAux
        vector_tag = 'restore_tag'
        v = 'disp_y'
        variable = 'resid_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampx]
        type = TagVectorAux
        vector_tag = 'restore_dampx_tag'
        v = 'disp_x'
        variable = 'resid_damp_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampy]
        type = TagVectorAux
        vector_tag = 'restore_dampy_tag'
        v = 'disp_y'
        variable = 'resid_damp_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressx]
        type = TagVectorAux
        vector_tag = 'restore_pressx_tag'
        v = 'disp_x'
        variable = 'resid_press_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressy]
        type = TagVectorAux
        vector_tag = 'restore_pressy_tag'
        v = 'disp_y'
        variable = 'resid_press_y'
        execute_on = 'TIMESTEP_END'
    []
    #calc velocity
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
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

[Problem]
    extra_tag_vectors = 'restore_tag restore_pressx_tag restore_pressy_tag restore_dampx_tag restore_dampy_tag'
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_x disp_y'
            planar_formulation = PLANE_STRAIN
            extra_vector_tags = 'restore_tag'
        [../]
        [../]
    [../]
[]

[Kernels]
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
        extra_vector_tags = 'restore_pressx_tag'
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
        extra_vector_tags = 'restore_pressy_tag'
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
    [./Reactionx]
        type = StiffPropDamping
        variable = disp_x
        component = 0
        extra_vector_tags = 'restore_dampx_tag'
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = disp_y
        component = 1
        extra_vector_tags = 'restore_dampy_tag'
    []
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

[UserObjects]
    #evalute system residual after system solve before auxkernels
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampx]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampy]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampy_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_pressx]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_pressy]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressy_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
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
    symbol_values = '10 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) PRESSURE FUNCTION ================
  [./pressure_block1]
    type = ParsedFunction
    expression = 'p*t^3*x*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'p tb tw R'
    symbol_values = '10 0 0.1 100.0'
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
    expression = '-p*t^3*exp(-(x^2 + y^2)/R^2)/10'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block0]
    type = ParsedFunction
    expression = '-p*t^3*exp(-(x^2 + y^2)/R^2)/50'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x_block1]
    type = ParsedFunction
    expression = 'p*t^3*exp(-(x^2 + y^2)/R^2)/10'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block1]
    type = ParsedFunction
    expression = 'p*t^3*exp(-(x^2 + y^2)/R^2)/50'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block0]
    type = ParsedFunction
    expression = 't*(-3*R^4*(20*dx*r + p*rf*t) - 10*R^2*t^2*(a*p*(R^2 - 2*x^2) + 4*dx*mu) + 80*dx*mu*t^2*x^2 + 20*t^2*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(10*R^4)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block0]
    type = ParsedFunction
    expression = 't*(-3*R^4*(100*dy*r + p*rf*t) + 100*R^2*t^2*(a*p*x*y - 2*dy*mu) + 400*dy*mu*t^2*y^2 + 100*t^2*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(50*R^4)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block1]
    type = ParsedFunction
    expression = 't*(3*R^4*(20*dx*r + p*rf*t) + 10*R^2*t^2*(a*p*(R^2 - 2*x^2) + 4*dx*mu) - 80*dx*mu*t^2*x^2 - 20*t^2*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(10*R^4)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block1]
    type = ParsedFunction
    expression = 't*(3*R^4*(100*dy*r + p*rf*t) + 100*R^2*t^2*(-a*p*x*y + 2*dy*mu) - 400*dy*mu*t^2*y^2 - 100*t^2*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/(50*R^4)'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block0]
    type = ParsedFunction
    expression = 't^2*(M*(150*a*(dx*x + dy*y) + 5*p*t*x + p*t*y) - 75*R^2*p*x)*exp(-(x^2 + y^2)/R^2)/(25*M*R^2)'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block1]
    type = ParsedFunction
    expression = 't^2*(-M*(150*a*(dx*x + dy*y) + 5*p*t*x + p*t*y) + 75*R^2*p*x)*exp(-(x^2 + y^2)/R^2)/(25*M*R^2)'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block0]
    type = ParsedFunction
    expression = 't*(-R^2*k*(60*dx*rf + 3*p*rp*t + 10*p*t^2) - R^2*mf*p*t^2 + 20*k*p*t^2*x^2)*exp(-(x^2 + y^2)/R^2)/(10*R^2*k)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block0]
    type = ParsedFunction
    expression = 't*(-3*R^2*k*(100*dy*rf + p*rp*t) - R^2*mf*p*t^2 + 100*k*p*t^2*x*y)*exp(-(x^2 + y^2)/R^2)/(50*R^2*k)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block1]
    type = ParsedFunction
    expression = 't*(R^2*k*(60*dx*rf + 3*p*rp*t + 10*p*t^2) + R^2*mf*p*t^2 - 20*k*p*t^2*x^2)*exp(-(x^2 + y^2)/R^2)/(10*R^2*k)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block1]
    type = ParsedFunction
    expression = 't*(3*R^2*k*(100*dy*rf + p*rp*t) + R^2*mf*p*t^2 - 100*k*p*t^2*x*y)*exp(-(x^2 + y^2)/R^2)/(50*R^2*k)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]
[]


[BCs]
    #assign displacement boundary condition
    [./matchval_primary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plusminus_scaled_x
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plusminus_scaled_x
        boundary = 'Block1_Block0'
    []
    [./matchval_primary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_plusminus_scaled_y
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_plusminus_scaled_y
        boundary = 'Block1_Block0'
    []
    [./match_vel_primary_y]
        type = MatchedValueBC
        variable = fluid_vel_y
        v = mms_fluid_vel_y
        boundary = 'Block0_Block1'
    []
    [./match_vel_secondary_y]
        type = MatchedValueBC
        variable = fluid_vel_y
        v = mms_fluid_vel_y
        boundary = 'Block1_Block0'
    []
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
  dt =  0.001
  end_time = 0.5
  #automatic_scaling = true
  [./TimeIntegrator]
   # type = NewmarkBeta
    type = CentralDifference
  [../]
[]

[Outputs]
    exodus = true
    time_step_interval = 5
[]


[MultiApps]
    #allocate transfer from mainApp to subApp
    [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'sub.i'
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
  []

[Transfers]
    #get displacement residuals from subApp to mainApp
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'disp_plusminus_sub_scaled_x disp_plusminus_sub_scaled_y traction_sub_strike traction_sub_normal sliprate_sub_strike sliprate_sub_normal slip_sub_strike slip_sub_normal statevar_sub'
        variable = 'disp_plusminus_scaled_x disp_plusminus_scaled_y traction_strike traction_normal sliprate_strike sliprate_normal slip_strike slip_normal statevar'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #push system residual vector from mainApp to subApp
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'resid_x resid_y resid_press_x resid_press_y resid_damp_x resid_damp_y p fluid_vel_x fluid_vel_y'
        variable = 'resid_sub_x resid_sub_y resid_press_sub_x resid_press_sub_y resid_damp_sub_x resid_damp_sub_y p_f fluid_vel_f_x fluid_vel_f_y'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
[]

# [Debug]
#     show_execution_order = ALWAYS
# []