##########################################################
# TPV205 benchmark
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

[Mesh]
    [./msh]
      type = FileMeshGenerator
    #   file =  '../../../meshgenerator/tpv205/tpv2053d_xyplane.msh'
      file =  '../../../meshgenerator/tpv205/tpv2053d_local_xyplane.msh'
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
[]

[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z'
    
    #damping ratio
    q = 1.0
    
    #characteristic length (m)
    Dc = 0.4
    
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
    [./resid_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_z]
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
    [./resid_damp_z]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_sw_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_damp_sw_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_sw_z]
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
[]

[Modules]
    [./TensorMechanics]
    [./Master]
        [./all]
            strain = SMALL
            add_variables = true
            generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
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
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_y]
      type = ProjectionAux
      variable = disp_slipweakening_y
      v = disp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_z]
      type = ProjectionAux
      variable = disp_slipweakening_z
      v = disp_z
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_x]
        type = CompVarRate
        variable = vel_slipweakening_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_slipweakening_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_slipweakening_z
      coupled = disp_z
      execute_on = 'TIMESTEP_BEGIN'
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
    [Residual_x]
      type = ProjectionAux
      variable = resid_slipweakening_x
      v = resid_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
      type = ProjectionAux
      variable = resid_slipweakening_y
      v = resid_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_z]
      type = ProjectionAux
      variable = resid_slipweakening_z
      v = resid_z
      execute_on = 'TIMESTEP_BEGIN'
    []
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
    #
    [Residual_damp_x]
      type = ProjectionAux
      variable = resid_damp_sw_x
      v = resid_damp_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_damp_y]
      type = ProjectionAux
      variable = resid_damp_sw_y
      v = resid_damp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_damp_z]
      type = ProjectionAux
      variable = resid_damp_sw_z
      v = resid_damp_z
      execute_on = 'TIMESTEP_BEGIN'
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
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [DynamicFricCoeff]
      type = FunctionAux
      variable = mu_d
      function = func_dynamic_friction_coeff_mud
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [elem_length]
      type = ConstantAux
      variable = elem_length
      value = 200
    []
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

# [InterfaceKernels]
#     [czm_interface_kernel_x]
#         type = FarmsCZM
#         variable = disp_x
#         neighbor_var = disp_x
#         boundary = 'Block2_Block3'
#     []
#     [czm_interface_kernel_y]
#         type = FarmsCZM
#         variable = disp_y
#         neighbor_var = disp_y
#         boundary = 'Block2_Block3'
#     []
#     [czm_interface_kernel_z]
#         type = FarmsCZM
#         variable = disp_z
#         neighbor_var = disp_z
#         boundary = 'Block2_Block3'
#     []
# []

[Materials]
    [stress]
        type = ComputeLinearElasticStress
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
        type = FarmsSlipWeakeningCZM
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        disp_slipweakening_z     = disp_slipweakening_z
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        reaction_slipweakening_z = resid_slipweakening_z
        reaction_damp_x = resid_damp_sw_x
        reaction_damp_y = resid_damp_sw_y
        reaction_damp_z = resid_damp_sw_z
        elem_length = elem_length
        mu_d = mu_d
        mu_s = mu_s
        boundary = 'Block2_Block3'
    [../]
    [./static_initial_stress_tensor_slipweakening]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor_slipweakening
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    [../]
[]

[Functions]
    [func_static_friction_coeff_mus]
       type = InitialStaticFrictionCoeff
    []
    #mud constant value
    [func_dynamic_friction_coeff_mud]
        type = ConstantFunction
        value = 0.525
    []
    #Note:restrict stress variation along the fault only
    #this function is used in czm only
    [./func_initial_stress_xx]
        type = ConstantFunction
        value = 0
    []
    [./func_initial_stress_xy]
        type = ConstantFunction
        value = 0
    []
    [./func_initial_stress_xz]
      type = InitialShearStress
      benchmark_type = tpv205
    []
    [./func_initial_stress_yy]
      type = ConstantFunction
      value = 0
    []
    [./func_initial_stress_yz]
        type = ConstantFunction
        value = 0
    []
    [./func_initial_stress_zz]
      type = ConstantFunction
      value = -120e6
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
    # num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]

[Outputs]
    exodus = true
    interval = 40
    # show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z traction_x traction_y traction_z jump_x jump_z mu_s tangent_jump_rate'
[]
