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
    displacements = 'disp_sub_x disp_sub_y'
    
    ##element length (m)
    len = 10
    
    ##rate-and-state coefficients
    f_o = 0.6
    rsf_a = 0.015
    rsf_b = 0.01
    rsf_L = 5e-4
    delta_o = 1

    ##initial normal traction (Pa)
    Tn_o = 120

    ##initial shear traction (Pa)
    Ts_o = 75

    ##initial sliprate (m/s)
    Vini = 1e-6

    ##initial state variable
    statevarini = 0.609124698e7
[]

[Variables]
    [disp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [./fluid_sub_vel_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./fluid_sub_vel_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./p_sub]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    #fault displacement residual from interfacekernel
    [disp_plusminus_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plusminus_sub_scaled_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_scaled_y]
        order = FIRST
        family = LAGRANGE
    []
    #side element volume 
    [element_side_volume]
        order = FIRST
        family = LAGRANGE
    []
    #residual received from mainApp
    [resid_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [resid_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_vel_f_x]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_vel_f_y]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_disp_f_x]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_disp_f_y]
        order = FIRST
        family = LAGRANGE
    []
    [p_f]
        order = FIRST
        family = LAGRANGE
    []
    #residual received from mainApp (damping)
    [resid_damp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [resid_damp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [resid_press_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [resid_press_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []
    #traction
    [traction_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #sliprate
    [sliprate_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [sliprate_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #slip
    [slip_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [slip_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #statevar
    [statevar_sub]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mms_sliprate]
        order = FIRST
        family = LAGRANGE
    [../]
    [./mms_shear]
        order = FIRST
        family = LAGRANGE
    [../]
    [./mms_normal]
        order = FIRST
        family = LAGRANGE
    [../]
    [./mms_pressure]
        order = FIRST
        family = LAGRANGE
    [../]
    [sliprate_sub_strikess]
        order = CONSTANT
        family = MONOMIAL
    []
        [sliprate_sub_strikesss]
        order = CONSTANT
        family = MONOMIAL
    []
    [disp_plus]
        order = FIRST
        family = LAGRANGE
    []
    [disp_minus]
        order = FIRST
        family = LAGRANGE
    []
    
[]

[AuxKernels]
    [./mms_shearss]
        type = FunctionAux
        function = theta_mms
        variable = mms_shear
        execute_on = 'INITIAL TIMESTEP_END'
    [../]
    [./mms_shearspresss]
        type = FunctionAux
        function = source_rsf
        variable = mms_pressure
        execute_on = 'INITIAL TIMESTEP_END'
    [../]
    [./mms_shearnoe]
        type = FunctionAux
        function = normal_stress_interface
        variable = mms_normal
        execute_on = 'INITIAL TIMESTEP_END'
    [../]
    [./mms_shearsssda]
        type = FunctionAux
        function = slip_rate_interface
        variable = mms_sliprate
        execute_on = 'INITIAL TIMESTEP_END'
    [../]
    [./TrapazoidalTimeIntegrator_x]
        type = VariableTimeIntegrationAux
        variable_to_integrate = fluid_vel_f_x
        variable = fluid_disp_f_x
        order = 2
    [../]
    [./TrapazoidalTimeIntegrator_y]
        type = VariableTimeIntegrationAux
        variable_to_integrate = fluid_vel_f_y
        variable = fluid_disp_f_y
        order = 2
    [../]
    #retrieve displacement vector by scaling the residuals
    [scale_disp_residual_x]
        type = ScaleVarAux
        variable = disp_plusminus_sub_scaled_x
        coupled = disp_plusminus_sub_x
        scale = element_side_volume
        execute_on = 'TIMESTEP_END'
    []
    [scale_disp_residual_y]
        type = ScaleVarAux
        variable = disp_plusminus_sub_scaled_y
        coupled = disp_plusminus_sub_y
        scale = element_side_volume
        execute_on = 'TIMESTEP_END'
    []
    #const element_side_volume
    [const_element_side_volume]
        type = ConstantAux
        variable = element_side_volume
        value = 10
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #retrieve fault displacement residual vector using tagging
    [restore_x]
        type = TagVectorAux
        vector_tag = 'restore_tag_x'
        v = 'disp_sub_x'
        variable = 'disp_plusminus_sub_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_y]
        type = TagVectorAux
        vector_tag = 'restore_tag_y'
        v = 'disp_sub_y'
        variable = 'disp_plusminus_sub_y'
        execute_on = 'TIMESTEP_END'
    []
    ##scale the displacement by surface area

    ##output
    [output_traction_strike]
        type = MaterialRealAux
        property = traction_strike_TSN
        variable = traction_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'LINEAR TIMESTEP_END'
    []
    [output_traction_normal]
        type = MaterialRealAux
        property = traction_normal_TSN
        variable = traction_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'LINEAR TIMESTEP_END'
    []
    [output_sliprate_strike]
        type = MaterialRealAux
        property = sliprate_strike
        variable = sliprate_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_sliprate_normal]
        type = MaterialRealAux
        property = sliprate_normal
        variable = sliprate_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_strike]
        type = MaterialRealAux
        property = slip_strike
        variable = slip_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_normal]
        type = MaterialRealAux
        property = slip_normal
        variable = slip_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_statevar]
        type = MaterialRealAux
        property = statevar
        variable = statevar_sub
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
        [output_dispPlus]
        type = MaterialRealAux
        property = alongfaultdisp_strike_plus
        variable = sliprate_sub_strikesss
        boundary = 'Block0_Block1'
        execute_on = 'INITIAL TIMESTEP_END'
    []
        [output_dispmins]
        type = MaterialRealAux
        property = alongfaultdisp_strike_minus
        variable = sliprate_sub_strikess
        boundary = 'Block0_Block1'
        execute_on = 'INITIAL TIMESTEP_END'
    []


[]


[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_sub_x disp_sub_y'
            planar_formulation = PLANE_STRAIN
        [../]
        [../]
    [../]
[]


[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
    [../]
[]

[Kernels]
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_sub_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_sub_y
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_x]
        type = CoupledFluidInertialForce
        variable = disp_sub_x
        fluid_vel = fluid_sub_vel_x
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_y]
        type = CoupledFluidInertialForce
        variable = disp_sub_y
        fluid_vel = fluid_sub_vel_y
        use_displaced_mesh = false
    [../]
    [./darcyflow_x]
        type = DynamicDarcyFlow2
        variable = fluid_sub_vel_x
        skeleton_acceleration = disp_sub_x
    [../]
    [./darcyflow_y]
        type = DynamicDarcyFlow2
        variable = fluid_sub_vel_y
        skeleton_acceleration = disp_sub_y
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_sub_x
        porepressure = p_sub
        component = 0
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_sub_y
        porepressure = p_sub
        component = 1
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling2
        variable = fluid_sub_vel_x
        porepressure = p_sub
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_sub_vel_y
        porepressure = p_sub
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p_sub
        displacements = 'disp_sub_x disp_sub_y'
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p_sub
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p_sub
        u = fluid_sub_vel_x
        v = fluid_sub_vel_y
        pressure = p_sub
    [../]
        [./source_term_x0]
      type = BodyForce
      variable = disp_sub_x
      function = source_x_block0
      block = 0
    [../]
    [./source_term_y0]
      type = BodyForce
      variable = disp_sub_y
      function = source_y_block0
      block = 0
    [../]
    [./source_term_p0]
      type = BodyForce
      variable = p_sub
      function = source_p_block0
      block = 0
    [../]
    [./source_term_vfx0]
      type = BodyForce
      variable = fluid_sub_vel_x
      function = source_vfx_block0
      block = 0
    [../]
    [./source_term_vfy0]
      type = BodyForce
      variable = fluid_sub_vel_y
      function = source_vfy_block0
      block = 0
    [../]

    [./source_term_x1]
      type = BodyForce
      variable = disp_sub_x
      function = source_x_block1
      block = 1
    [../]
    [./source_term_y1]
      type = BodyForce
      variable = disp_sub_y
      function = source_y_block1
      block = 1
    [../]
    [./source_term_p1]
      type = BodyForce
      variable = p_sub
      function = source_p_block1
      block = 1
    [../]
    [./source_term_vfx1]
      type = BodyForce
      variable = fluid_sub_vel_x
      function = source_vfx_block1
      block = 1
    [../]
    [./source_term_vfy1]
      type = BodyForce
      variable = fluid_sub_vel_y
      function = source_vfy_block1
      block = 1
    [../]
[]


[Problem]
    extra_tag_vectors = 'restore_tag_x restore_tag_y'
[]


[InterfaceKernels]
   #apply displacement prediction and retrieve its residuals
    [./ratestate_x]
        type = PoroRateStateInterfaceKernelGlobalx
        variable = disp_sub_x
        neighbor_var = disp_sub_x
        extra_vector_tags = 'restore_tag_x'
        boundary = 'Block0_Block1'
        y_var = disp_sub_y
        vfx_var = fluid_sub_vel_x
       vfy_var = fluid_sub_vel_y
        pressure_var = p_sub
    []
    [./ratestate_y]
        type = PoroRateStateInterfaceKernelGlobaly
        variable = disp_sub_y
        neighbor_var = disp_sub_y
        extra_vector_tags = 'restore_tag_y'
       boundary = 'Block0_Block1'
        x_var = disp_sub_x
        vfx_var = fluid_sub_vel_x
      vfy_var = fluid_sub_vel_y
        pressure_var = p_sub
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
    [./czm_mat]
        type = MMSPoroRateStateFrictionLaw2DAsBC
        reaction_rsf_x  = resid_sub_x
        reaction_rsf_y  = resid_sub_y
        reaction_rsf_pressure_x = resid_press_sub_x
        reaction_rsf_pressure_y = resid_press_sub_y
        reaction_damp_x = resid_damp_sub_x
        reaction_damp_y = resid_damp_sub_y
        reaction_pressdamp_x = resid_damp_sub_x
        reaction_pressdamp_y = resid_damp_sub_y
        interface_pressure = p_f 
        fluid_vel_x = fluid_vel_f_x
        fluid_vel_y = fluid_vel_f_y
        fluid_disp_x = fluid_disp_f_x
        fluid_disp_y = fluid_disp_f_y
        nodal_area = nodal_area
        mms_sliprate = mms_sliprate
        mms_shear = mms_shear
        mms_normal = mms_normal
        mms_pressure = mms_pressure
        boundary = 'Block0_Block1'
        output_properties = 'sliprate_strike slip_strike statevar traction_strike traction_normal alongfaultdisp_strike_plus alongfaultdisp_strike_minus'
        outputs = exodus
    [../]
[]

[UserObjects]
   # compute element side volume (using CONTACT modulus)
     [element_side_volume]
         type = NodalArea
         variable = element_side_volume
         boundary = 'Block0_Block1 Block1_Block0'
         execute_on = 'initial TIMESTEP_BEGIN'
    []
    [recompute_residual_tag_x]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag_x'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_y]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag_y'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = Block0_Block1
      execute_on = 'initial TIMESTEP_BEGIN'
    [../]
[]

[Functions]
  # ================ BLOCK 0 (LOWER DOMAIN) DISPLACEMENT FUNCTIONS ================
  [./displacement_x_block0]
    type = ParsedFunction
    expression = '-dx*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./displacement_y_block0]
    type = ParsedFunction
    expression = '-dy*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) DISPLACEMENT FUNCTIONS ================
  [./displacement_x_block1]
    type = ParsedFunction
    expression = 'dx*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./displacement_y_block1]
    type = ParsedFunction
    expression = 'dy*t^2*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) PRESSURE FUNCTION ================
  [./pressure_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'p tb tw R'
    symbol_values = '10 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) PRESSURE FUNCTION ================
  [./pressure_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'p tb tw R'
    symbol_values = '10 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) VELOCITY FUNCTIONS ================
  [./velocity_x_block0]
    type = ParsedFunction
    expression = '-2*dx*t*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./velocity_y_block0]
    type = ParsedFunction
    expression = '-2*dy*t*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) VELOCITY FUNCTIONS ================
  [./velocity_x_block1]
    type = ParsedFunction
    expression = '2*dx*t*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./velocity_y_block1]
    type = ParsedFunction
    expression = '2*dy*t*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) DARCY VELOCITY FUNCTIONS ================
  [./darcy_velocity_x_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  [./darcy_velocity_y_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'k mf p tb tw R rf rp dx dy'
    symbol_values = '1.5 1 10 0 0.1 100.0 1 20 5.2 0.1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block0]
    type = ParsedFunction
    expression = '2*(-R^4*dx*r - 2*R^2*dx*mu*t^2 + 4*dx*mu*t^2*x^2 + t^2*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp((-x^2 - y^2)/R^2)/R^4'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block0]
    type = ParsedFunction
    expression = '2*(-R^4*dy*r - 2*R^2*dy*mu*t^2 + 4*dy*mu*t^2*y^2 + t^2*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp((-x^2 - y^2)/R^2)/R^4'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR ELASTICITY ================
  [./source_x_block1]
    type = ParsedFunction
    expression = '2*(R^4*dx*r + 2*R^2*dx*mu*t^2 - 4*dx*mu*t^2*x^2 - t^2*(l*(-R^2*dx + 2*x*(dx*x + dy*y)) + mu*(-1.0*R^2*dx + 2.0*y*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/R^4'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  [./source_y_block1]
    type = ParsedFunction
    expression = '2*(R^4*dy*r + 2*R^2*dy*mu*t^2 - 4*dy*mu*t^2*y^2 - t^2*(l*(-R^2*dy + 2*y*(dx*x + dy*y)) + mu*(-1.0*R^2*dy + 2.0*x*(dx*y + dy*x))))*exp(-(x^2 + y^2)/R^2)/R^4'
    symbol_names = 'dx dy p tb tw R l mu a r rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.5 0.75 0.6 2.5 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block0]
    type = ParsedFunction
    expression = '4*a*t*(dx*x + dy*y)*exp(-(x^2 + y^2)/R^2)/R^2'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR PRESSURE ================
  [./source_p_block1]
    type = ParsedFunction
    expression = '-4*a*t*(dx*x + dy*y)*exp(-(x^2 + y^2)/R^2)/R^2'
    symbol_names = 'dx dy p tb tw R a M rf rp k mf'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.6 4.7059 1 20 1.5 1'
  [../]

  # ================ BLOCK 0 (LOWER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block0]
    type = ParsedFunction
    expression = '-2*dx*rf*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block0]
    type = ParsedFunction
    expression = '-2*dy*rf*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  # ================ BLOCK 1 (UPPER DOMAIN) SOURCE TERMS FOR DARCY VELOCITY ================
  [./source_vfx_block1]
    type = ParsedFunction
    expression = '2*dx*rf*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  [./source_vfy_block1]
    type = ParsedFunction
    expression = '2*dy*rf*exp(-(x^2 + y^2)/R^2)'
    symbol_names = 'dx dy p tb tw R mu l a r rf rp k mf M kp'
    symbol_values = '5.2 0.1 10 0 0.1 100.0 0.75 0.5 0.6 2.5 1 20 1.5 1 4.7059 1.5'
  [../]

  # ================ INTERFACE VALUES (y=0) ================
  [./slip_rate_interface]
    type = ParsedFunction
    expression = '4*dx*t*exp(-x^2/R^2)'
    symbol_names = 'dx dy tb tw R'
    symbol_values = '5.2 0.1 0 0.1 100.0'
  [../]

  [./normal_stress_interface]
    type = ParsedFunction
    expression = '(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))*exp(-x^2/R^2)/(R^4*dt^2)'
    symbol_names =  'dx dy mu l tb tw R Lel dt r rf mf k rp p a Tno Tso'
    symbol_values = '5.2 0.1 0.75 0.5 0 0.1 100.0 10 0.001 2.5 1 1 1.5 20 10 0.6 120 75'
  [../]

  [./shear_stress_interface]
    type = ParsedFunction
    expression = '(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))*exp(-x^2/R^2)/(R^4*dt)'
    symbol_names =  'dx dy mu l tb tw R Lel dt r rf mf k rp p a Tno Tso'
    symbol_values = '5.2 0.1 0.75 0.5 0 0.1 100.0 10 0.001 2.5 1 1 1.5 20 10 0.6 120 75'
  [../]

  [./max_pressure_interface]
    type = ParsedFunction
    expression = 'abs(0)'
    symbol_names = 'p tb tw R'
    symbol_values = '10 0 0.1 100.0'
  [../]

  # ================ RATE AND STATE FRICTION SOURCE TERM ================
  [./source_rsf]
    type = ParsedFunction
    expression = 'if(t<1e-6,0,L*(rsfa*(dt*(Lel*R^4*dx*r + 2*Lel*dt*dx*t*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))) + dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))*(0.5*Lel*R^4*dy*r*t + 0.25*Lel*R^4*dy*r*(4*dt + 2.0*t) + 4*Lel*dt^2*dy*t*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)))/(4*rsfa*(-0.125*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2)/2)^2))*cosh(dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))))/sinh(dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2)))) - rsfa/t)*exp((-fo - rsfa*log(4*dx*t*exp(-x^2/R^2)/Vo) + rsfa*log(-2*sinh(dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))))))/rsfb)/(Vo*rsfb) - 1 + 4*dx*t*exp(-x^2/R^2)*exp((-fo - rsfa*log(4*dx*t*exp(-x^2/R^2)/Vo) + rsfa*log(-2*sinh(dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))))))/rsfb)/Vo)'
    symbol_names = 'dx dy mu l tb tw R Lel dt r rf mf k rp p a rsfa rsfb Vo fo L Tno Tso'
    symbol_values = '5.2 0.1 0.75 0.5 0 0.1 100.0 10 0.001 2.5 1 1 1.5 20 10 0.6 0.015 0.01 1 0.6 5e-4 120 75'
  [../]

  [./theta_mms]
    type = ParsedFunction
    expression = 'if(t<1e-6,0,L*exp((-fo - rsfa*log(4*dx*t*exp(-x^2/R^2)/Vo) + rsfa*log(-2*sinh(dt*(Lel*R^4*dx*r*t + Lel*dt*dx*t^2*(6.0*R^2*mu + 2*l*(R^2 - 2*x^2) - 8*mu*x^2) + R^4*Tso*dt*exp(x^2/R^2))/(rsfa*(-0.25*Lel*R^4*dy*r*t*(4*dt + 2.0*t) - 2*Lel*dt^2*dy*t^2*(R^2*(l + 2*mu) + mu*(1.0*R^2 - 2.0*x^2)) - R^4*Tno*dt^2*exp(x^2/R^2))))))/rsfb)/Vo)'
    symbol_names = 'dx dy mu l tb tw R Lel dt r rf mf k rp p a rsfa rsfb Vo fo L Tno Tso'
    symbol_values = '5.2 0.1 0.75 0.5 0 0.1 100.0 10 0.001 2.5 1 1 1.5 20 10 0.6 0.008 0.012 1e-6 0.6 0.02 120 75'
  [../]

[]

[ICs]
  # Block 0 (lower domain) initial conditions
  
  [./disp_x_ic_block0]
    type = FunctionIC
    variable = disp_sub_x
    function = displacement_x_block0
    block = 0
  [../]
  
  [./disp_y_ic_block0]
    type = FunctionIC
    variable = disp_sub_y
    function = displacement_y_block0
    block = 0
  [../]
  
  [./p_ic_block0]
    type = FunctionIC
    variable = p_sub
    function = pressure_block0
    block = 0
  [../]
  
  [./fluid_vel_x_ic_block0]
    type = FunctionIC
    variable = fluid_sub_vel_x
    function = darcy_velocity_x_block0
    block = 0
  [../]
  
  [./fluid_vel_y_ic_block0]
    type = FunctionIC
    variable = fluid_sub_vel_y
    function = darcy_velocity_y_block0
    block = 0
  [../]
  
  # Block 1 (upper domain) initial conditions
  [./disp_x_ic_block1]
    type = FunctionIC
    variable = disp_sub_x
    function = displacement_x_block1
    block = 1
  [../]
  
  [./disp_y_ic_block1]
    type = FunctionIC
    variable = disp_sub_y
    function = displacement_y_block1
    block = 1
  [../]
  
  [./p_ic_block1]
    type = FunctionIC
    variable = p_sub
    function = pressure_block1
    block = 1
  [../]
  
  [./fluid_vel_x_ic_block1]
    type = FunctionIC
    variable = fluid_sub_vel_x
    function = darcy_velocity_x_block1
    block = 1
  [../]
  
  [./fluid_vel_y_ic_block1]
    type = FunctionIC
    variable = fluid_sub_vel_y
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
  #automatic_scaling = true
  [./TimeIntegrator]
   # type = NewmarkBeta
    type = CentralDifference
    solve_type = lumped
  [../]
[]


[Outputs]
    exodus = true
    interval = 1
[]