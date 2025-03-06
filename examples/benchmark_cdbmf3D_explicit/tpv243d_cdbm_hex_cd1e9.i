##########################################################
# TPV205 benchmark with CDBM
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

  [Mesh]
    [msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -15000
      xmax = 15000
      ymin = -20000
      ymax = 0
      zmin = -8000
      zmax = 8000
      nx = 150
      ny = 100
      nz = 80
      subdomain_ids = 1
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -13000 & x <= 13000 & y > -15000 & z < 0'
      block_id = 2
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x > -13000 & x < 13000 & y > -15000 & z > 0'
        block_id = 3
    [] 
    [./split_1]
      type = BreakMeshByBlockGenerator
      input = new_block_2
      split_interface = true
      block_pairs = '2 3'
    []      
    [./sidesets]
      input = split_1
      type = SideSetsFromNormalsGenerator
      normals = '-1 0 0
                  1 0 0
                  0 -1 0
                  0 1 0
                  0 0 -1
                  0 0 1'
      new_boundary = 'left right bottom top back front'
    []       
  []

[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z'
    
    #damping ratio
    q = 0.4
    
    #characteristic length (m) #this gives ~200m for resolve L_f
    Dc = 0.3

    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -1.66
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -1.66
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 1e9

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    # energy ratio
    chi = 0.7

    #
    D = 0
    
[]

[AuxVariables]
    [./resid_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_z]
      order = FIRST
      family = LAGRANGE
    []
    #restoration force for damping (tag after solve)
    [./resid_damp_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_damp_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_z]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_z]
      order = FIRST
      family = LAGRANGE
    []
    [./vel_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_z]
      order = FIRST
      family = LAGRANGE
    []
    [./accel_slipweakening_x]
      order = FIRST
      family = LAGRANGE
    []
    [./accel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_slipweakening_z]
      order = FIRST
      family = LAGRANGE
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_d]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_shear_stress]
        order = FIRST
        family = LAGRANGE
    []
    #
    [./tangent_jump]
      order = CONSTANT
      family = MONOMIAL
    []
    [./tangent_jump_rate]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [./elem_length]
        order = CONSTANT
        family = MONOMIAL
    [../]
    #
    [./jump_x]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_y]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_z]
        order = CONSTANT
        family = MONOMIAL         
    []
    #
    [jump_rate_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [traction_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #obtain parameters from MaterialRealAux, pass parameters to subApp
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [initial_damage_aux]
      order = FIRST
      family = MONOMIAL
    []
    [alpha_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [B_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [./check_function_initial_stress_xx]
      order = FIRST
      family = LAGRANGE      
    []
    [./check_function_initial_stress_xy]
      order = FIRST
      family = LAGRANGE      
    []
    [./check_function_initial_stress_xz]
      order = FIRST
      family = LAGRANGE      
    []
    [./check_function_initial_stress_yy]
      order = FIRST
      family = LAGRANGE      
    []
    [./check_function_initial_stress_yz]
      order = FIRST
      family = LAGRANGE      
    []
    [./check_function_initial_stress_zz]
      order = FIRST
      family = LAGRANGE      
    []
    #
    [cohesion_aux]
      order = FIRST
      family = LAGRANGE         
    []
    [forced_rupture_time_aux]
      order = FIRST
      family = LAGRANGE      
    []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all]
        strain = SMALL
        add_variables = true
        extra_vector_tags = 'restore_tag'
        eigenstrain_names = static_initial_strain_tensor
      [../]
    [../]
  [../]
[]

[Problem]
    extra_tag_vectors = 'restore_tag restore_dampx_tag restore_dampy_tag restore_dampz_tag'
[]

[AuxKernels]
    [Displacment_x]
      type = ProjectionAux
      variable = disp_slipweakening_x
      v = disp_x
      execute_on = 'TIMESTEP_END'
    []
    [Displacement_y]
      type = ProjectionAux
      variable = disp_slipweakening_y
      v = disp_y
      execute_on = 'TIMESTEP_END'
    []
    [Displacement_z]
      type = ProjectionAux
      variable = disp_slipweakening_z
      v = disp_z
      execute_on = 'TIMESTEP_END'
    []
    [Vel_x]
        type = CompVarRate
        variable = vel_slipweakening_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_slipweakening_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_slipweakening_z
      coupled = disp_z
      execute_on = 'TIMESTEP_END'
    []
    #
    [XJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    #
    [XJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []    
    #
    [TractionX]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'        
    []
    [TractionY]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [TractionZ]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []        
    #
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
    [restore_z]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_z'
      variable = 'resid_z'
      execute_on = 'TIMESTEP_END'
    []
    #damping
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
    [restore_dampz]
        type = TagVectorAux
        vector_tag = 'restore_dampz_tag'
        v = 'disp_z'
        variable = 'resid_damp_z'
        execute_on = 'TIMESTEP_END'
    []
    [StaticFricCoeff]
      type = FunctionAux
      variable = mu_s
      function = func_static_friction_coeff_mus
      execute_on = 'INITIAL TIMESTEP_END'
    []
    [DynamicFricCoeff]
      type = FunctionAux
      variable = mu_d
      function = func_dynamic_friction_coeff_mud
      execute_on = 'INITIAL TIMESTEP_END'
    []
    [elem_length]
      type = ConstantAux
      variable = elem_length
      value = 200
    []
    #obtain parameters from MaterialRealAux
    [get_initial_damage]
      type = MaterialRealAux
      variable = initial_damage_aux
      property = initial_damage
    []
    [get_damage]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
    []
    [get_strain_invariant_ratio]
        type = MaterialRealAux
        variable = xi_aux
        property = xi
    []
    [get_B]
        type = MaterialRealAux
        variable = B_aux
        property = B
    []
    #
    #cohesion
    [cohesion]
      type = FunctionAux
      variable = cohesion_aux
      function = func_initial_cohesion
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #forced_rupture_time
    [forced_rupture_time]
      type = FunctionAux
      variable = forced_rupture_time_aux
      function = func_forced_rupture_time
      execute_on = 'INITIAL TIMESTEP_BEGIN'      
    []
    #
    [checkxx]
      type = FunctionAux
      variable = check_function_initial_stress_xx
      function = func_initial_stress_xx
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [checkxy]
      type = FunctionAux
      variable = check_function_initial_stress_xy
      function = func_initial_stress_xy
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [checkxz]
      type = FunctionAux
      variable = check_function_initial_stress_xz
      function = func_initial_stress_xz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [checkyy]
      type = FunctionAux
      variable = check_function_initial_stress_yy
      function = func_initial_stress_yy
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [checkyz]
      type = FunctionAux
      variable = check_function_initial_stress_yz
      function = func_initial_stress_yz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [checkzz]
      type = FunctionAux
      variable = check_function_initial_stress_zz
      function = func_initial_stress_zz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [] 
[]

[Kernels]
    [./inertia_x]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_x
    []
    [./inertia_y]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_y
    []
    [./inertia_z]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_z
    []
    [./Reactionx]
      type = StiffPropDamping
      variable = 'disp_x'
      component = '0'
      extra_vector_tags = restore_dampx_tag
    []
    [./Reactiony]
      type = StiffPropDamping
      variable = 'disp_y'
      component = '1'
      extra_vector_tags = restore_dampy_tag
    []
    [./Reactionz]
      type = StiffPropDamping
      variable = 'disp_z'
      component = '2'
      extra_vector_tags = restore_dampz_tag
    []
[]

[InterfaceKernels]
    [czm_interface_kernel_x]
        type = FarmsCZM
        variable = disp_x
        neighbor_var = disp_x
        boundary = 'Block2_Block3'
    []
    [czm_interface_kernel_y]
        type = FarmsCZM
        variable = disp_y
        neighbor_var = disp_y
        boundary = 'Block2_Block3'
    []
    [czm_interface_kernel_z]
        type = FarmsCZM
        variable = disp_z
        neighbor_var = disp_z
        boundary = 'Block2_Block3'
    []
[]

[Materials]
    #damage breakage model
    [stress_medium]
      type = ComputeDamageBreakageStress3DDynamicCDBMF
      alpha_grad_x = alpha_grad_x
      alpha_grad_y = alpha_grad_y
      alpha_grad_z = alpha_grad_z
      output_properties = 'eps_p eps_e I1 I2 stress xi alpha_damagedvar B'
      outputs = exodus
  [] 
    [density]
        type = GenericConstantMaterial
        prop_names = 'density nonADdensity'
        prop_values = '2700 2700'
    [] 
    [./czm_mat]
        type = FarmsSlipWeakeningCZMcdbm
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        disp_slipweakening_z     = disp_slipweakening_z
        vel_slipweakening_x      = vel_slipweakening_x
        vel_slipweakening_y      = vel_slipweakening_y
        vel_slipweakening_z      = vel_slipweakening_z
        reaction_slipweakening_x = resid_x
        reaction_slipweakening_y = resid_y
        reaction_slipweakening_z = resid_z
        reaction_damp_x = resid_damp_x
        reaction_damp_y = resid_damp_y
        reaction_damp_z = resid_damp_z
        elem_length = elem_length
        mu_d = mu_d
        mu_s = mu_s
        cohesion = cohesion_aux
        forced_rupture_time = forced_rupture_time_aux
        boundary = 'Block2_Block3'
    [../]
    [./static_initial_stress_tensor_slipweakening]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor_slipweakening
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    [../]
    [./static_initial_stress_tensor]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    [../]
    #compute eigenstrain
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      lambda = 30e9
      shear_modulus = 30e9
    [../]
    [./strain_from_initial_stress]
      type = ComputeEigenstrainFromInitialStress
      initial_stress = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                        func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                        func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
      eigenstrain_name = static_initial_strain_tensor
    [../]   
    #dummy material
    [./dummy_mat]
      type = GenericConstantMaterial
      prop_names = 'initial_damage initial_breakage initial_shear_stress damage_perturbation shear_stress_perturbation'
      prop_values = '0.0 0.0 0.0 0.0 0.0' 
    [] 
[]

[Functions]
  [func_static_friction_coeff_mus]
     type = ConstantFunction
     value = 0.18
  []
  #mud constant value: 0.4
  [func_dynamic_friction_coeff_mud]
    type = ConstantFunction
    value = 0.12
  []
  #Note:restrict stress variation along the fault only
  #this function is used in czm only
  [./func_initial_stress_xx]
    type = ParsedFunction
    expression = 'if(-y<15600,  1 * (0.926793 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) ) - (1000 * 9.8 * (-y))), 1 * (-2700 * 9.8 * (-y)))'
    # expression = 'if(-y<15600,  1 * (2.0 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) ) - (1000 * 9.8 * (-y))), 1 * (-2700 * 9.8 * (-y)))'
    # expression = 0
  []
  [./func_initial_stress_xy]
    type = ParsedFunction
    expression = 0
  []
  [./func_initial_stress_xz]
    type = ParsedFunction
    # expression = 'if(-y<15600, 1 * (-0.169029 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) )), 0.0)'
    expression = 'if(-y<15600, 1 * (-0.169029 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) )), 0.0)'
    # type = InitialShearStressCDBM
    # benchmark_type = 'tpv24'
  []  
  [./func_initial_stress_yy]
    type = ParsedFunction
    expression = '-2700 * 9.8 * (-y)'
  []
  [./func_initial_stress_yz]
    type = ParsedFunction
    expression = 0
  []  
  [./func_initial_stress_zz]
    type = ParsedFunction
    # expression = 'if(-y<15600,  1 * (1.073206 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) ) - (1000 * 9.8 * (-y))), 1 * (-2700 * 9.8 * (-y)))'  
    expression = 'if(-y<15600,  1 * (1.073206 * ( (-2700 * 9.8 * (-y)) + (1000 * 9.8 * (-y)) ) - (1000 * 9.8 * (-y))), 1 * (-2700 * 9.8 * (-y)))'  
  []
  [./func_initial_cohesion]
    type = InitialCohesionTPV243D
  []
  [./func_forced_rupture_time]
    type = ForcedRuptureTimeTPV243D
    loc_x = 0
    loc_y = -4000
    loc_z = 0
    r_crit = 4000
    Vs = 3333
  []
[]

[UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    #damping
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
    [recompute_residual_tag_dampz]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampz_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
[]

[Executioner]
    type = Transient
    dt = 0.0025
    end_time = 12.0
    # num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]

[Outputs]
    exodus = true
    time_step_interval = 40
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z alpha_damagedvar B xi'
[]