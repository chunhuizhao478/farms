##########################################################
# TPV205 benchmark with CDBM
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

  [Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -18000
      xmax = 18000
      ymin = -20000
      ymax = 0
      zmin = -10000
      zmax = 10000
      nx = 180
      ny = 100
      nz = 100
      subdomain_ids = 1
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z < 0'
      block_id = 2
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z > 0'
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
    q = 1.0
    
    #characteristic length (m)
    Dc = 0.4

    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 3.204e10
      
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 3.204e10
  
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.99
  
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.99
  
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
  
    ##Compute gamma_damaged_r, xi_1
    #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
    #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
    #check struct_param.m 
  
    #coefficient of damage solid modulus
    gamma_damaged_r = 3.7150e10
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.8248
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.75
    a0 = 7.4289e9
    a1 = -2.214e10
    a2 = 2.0929e10
    a3 = -6.0672e9
  
    #diffusion coefficient #for structural stress coupling
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
    [./B_old]
      order = FIRST
      family = LAGRANGE
    []
    [./xi_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./I2_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./lambda_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./gamma_old]
        order = CONSTANT
        family = MONOMIAL
    []
    #updated alpha, B
    [./alpha_in]
      order = CONSTANT
      family = MONOMIAL
    []
    [./B_in]
      order = CONSTANT
      family = MONOMIAL
    []
    #high-order-dummy
    [./alpha_in_dummy]
      order = FIRST
      family = MONOMIAL
    []
    [./B_in_dummy]
      order = FIRST
      family = MONOMIAL
    []
    #grad_alpha
    [./alpha_grad_x]
      order = CONSTANT
      family = MONOMIAL
    []
    [./alpha_grad_y]
      order = CONSTANT
      family = MONOMIAL
    []
    [./alpha_grad_z]
      order = CONSTANT
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
    [get_xi_old]
      type = MaterialRealAux
      property = xi
      variable = xi_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_I2_old]
        type = MaterialRealAux
        property = I2
        variable = I2_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_mu_old]
        type = MaterialRealAux
        property = shear_modulus
        variable = mu_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_lambda_old]
        type = MaterialRealAux
        property = lambda
        variable = lambda_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_gamma_old]
        type = MaterialRealAux
        property = gamma_damaged
        variable = gamma_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
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
      type = ComputeDamageBreakageStress3D
      alpha_in = alpha_in_dummy
      B_in = B_in_dummy
      alpha_grad_x = alpha_grad_x
      alpha_grad_y = alpha_grad_y
      alpha_grad_z = alpha_grad_z
      output_properties = 'eps_p eps_e eps_total I1'
      outputs = exodus
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    [elasticity]
      type = ComputeIsotropicElasticityTensor
      shear_modulus = 32.04e9
      lambda = 32.04e9
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
[]

[Functions]
  [func_static_friction_coeff_mus]
     type = InitialStaticFrictionCoeffTPV243D
  []
  #mud constant value: 0.4
  [func_dynamic_friction_coeff_mud]
      type = ConstantFunction
      value = 0.12
  []
  #Note:restrict stress variation along the fault only
  #this function is used in czm only
  [./func_initial_stress_xx]
    type = InitialStressTPV243D
    i = 1
    j = 1
  []
  [./func_initial_stress_xy]
    type = InitialStressTPV243D
    i = 1
    j = 2
  []
  [./func_initial_stress_xz]
    type = InitialStressTPV243D
    i = 1
    j = 3
  []  
  [./func_initial_stress_yy]
    type = InitialStressTPV243D
    i = 2
    j = 2
  []
  [./func_initial_stress_yz]
    type = InitialStressTPV243D
    i = 2
    j = 3
  []  
  [./func_initial_stress_zz]
    type = InitialStressTPV243D
    i = 3
    j = 3
  []
  [./func_initial_cohesion]
    type = InitialCohesionTPV243D
  []
  [./func_forced_rupture_time]
    type = ForcedRuptureTimeTPV243D
    loc_x = 0
    loc_y = -7500
    loc_z = 0
    r_crit = 4000
    Vs = 3464
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
    dt = 0.0002
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
    time_step_interval = 50
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z traction_x traction_y traction_z jump_x jump_y jump_z jump_rate_x jump_rate_y jump_rate_z mu_s alpha_in B_in xi_old'
[]

[MultiApps]
  [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'tpv2053d_cdbm_hex_sub.i'
      execute_on = 'TIMESTEP_BEGIN'
  [../]
[]

[Transfers]
  [pull_resid]
      type = MultiAppCopyTransfer
      from_multi_app = sub_app
      source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub alpha_checked_dummy B_checked_dummy'
      variable = 'alpha_in B_in alpha_grad_x alpha_grad_y alpha_in_dummy B_in_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
  #we actually don't need to pass alpha and B
  [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old alpha_in_dummy B_in_dummy'
      variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old alpha_old_dummy B_old_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
[]
