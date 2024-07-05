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
        add_interface_on_two_sides = true
    []     
[]

[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_sub_x disp_sub_y disp_sub_z'
    
    #characteristic length (m)
    Dc = 0.4
    
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
    [disp_sub_z]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    #fault displacement residual from interfacekernel
    [disp_plusminus_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_z]
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
    [disp_plusminus_sub_scaled_z]
        order = FIRST
        family = LAGRANGE
    []
    #
    [disp_sw_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sw_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sw_sub_z]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_sub_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_sub_z]
      order = FIRST
      family = LAGRANGE
    []
    #restoration force for damping (tag after solve)
    [./resid_damp_sub_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_damp_sub_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_sub_z]
        order = FIRST
        family = LAGRANGE
    [../] 
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
    [./jump_sub_x]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_sub_y]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_sub_z]
        order = CONSTANT
        family = MONOMIAL         
    []
    #
    [jump_rate_sub_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_sub_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_sub_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [traction_sub_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_sub_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_sub_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [displacement_plus_x]
        order = CONSTANT
        family = MONOMIAL    
    [] 
    [displacement_plus_x_secondary]
        order = CONSTANT
        family = MONOMIAL    
    []
    [displacement_minus_x]
        order = CONSTANT
        family = MONOMIAL    
    []
    [displacement_minus_x_nodal]
        order = FIRST
        family = LAGRANGE       
    []
    [displacement_minus_x_nodal_secondary]
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
        [../]
    [../]
    [../]
[]

[Problem]
    extra_tag_vectors = 'restore_tag_x restore_tag_y restore_tag_z'
[]

[AuxKernels]
    #retrieve displacement vector by scaling the residuals
    [scale_disp_residual_x]
        type = FarmsScaleVarAux
        variable = disp_plusminus_sub_scaled_x
        coupled = disp_plusminus_sub_x
        execute_on = 'TIMESTEP_END'
    []
    [scale_disp_residual_y]
        type = FarmsScaleVarAux
        variable = disp_plusminus_sub_scaled_y
        coupled = disp_plusminus_sub_y
        execute_on = 'TIMESTEP_END'
    []
    [scale_disp_residual_z]
        type = FarmsScaleVarAux
        variable = disp_plusminus_sub_scaled_z
        coupled = disp_plusminus_sub_z
        execute_on = 'TIMESTEP_END'
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
    [restore_z]
        type = TagVectorAux
        vector_tag = 'restore_tag_z'
        v = 'disp_sub_z'
        variable = 'disp_plusminus_sub_z'
        execute_on = 'TIMESTEP_END'
    []    
    #
    [XJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_sub_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_sub_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_sub_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    #
    [XJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_sub_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_sub_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_sub_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []    
    #
    [TractionX]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_sub_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'        
    []
    [TractionY]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_sub_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [TractionZ]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_sub_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    [] 
    #
    [XDispPlus]
        type = MaterialRealVectorValueAux
        property = displacements_plus_global
        variable = displacement_plus_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [XDispMinus]
        type = MaterialRealVectorValueAux
        property = displacements_minus_global
        variable = displacement_minus_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    #      
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
[]

[InterfaceKernels]
    #apply displacement prediction and retrieve its residuals
    [./czm_x]
        type = FarmsInterfaceKernelGlobalx
        variable = disp_sub_x
        neighbor_var = disp_sub_x
        extra_vector_tags = 'restore_tag_x'
        boundary = 'Block2_Block3'
    []
    [./czm_y]
        type = FarmsInterfaceKernelGlobaly
        variable = disp_sub_y
        neighbor_var = disp_sub_y
        extra_vector_tags = 'restore_tag_y'
        boundary = 'Block2_Block3'
    []
    [./czm_z]
        type = FarmsInterfaceKernelGlobalz
        variable = disp_sub_z
        neighbor_var = disp_sub_z
        extra_vector_tags = 'restore_tag_z'
        boundary = 'Block2_Block3'
    []
[]

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
        disp_slipweakening_x     = disp_sw_sub_x
        disp_slipweakening_y     = disp_sw_sub_y
        disp_slipweakening_z     = disp_sw_sub_z
        reaction_slipweakening_x = resid_sub_x
        reaction_slipweakening_y = resid_sub_y
        reaction_slipweakening_z = resid_sub_z
        reaction_damp_x = resid_damp_sub_x
        reaction_damp_y = resid_damp_sub_y
        reaction_damp_z = resid_damp_sub_z
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
    [recompute_residual_tag_z]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag_z'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
[]

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
[]
