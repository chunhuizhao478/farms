##########################################################
# TPV14 benchmark
##########################################################

[Mesh]
    [./msh]
      type = FileMeshGenerator
      # file =  '../../meshgenerator/tpv14/tpv143d_400m.msh'
      file =  '../../meshgenerator/tpv14/tpv143d_200m.msh'
      # file =  '../../meshgenerator/tpv14/tpv143d_100m.msh'
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -16000 & x <= 12000 & y < 0 & z >= -15000'
      block_id = 2
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x >= -16000 & x <= 12000 & y > 0 & z >= -15000'
        block_id = 3
    []    
    [./new_block_3]
        type = ParsedSubdomainMeshGenerator
        input = new_block_2
        combinatorial_geometry = '0.5773505384 * x + y < 0 & z > -15000 & x >= 100 & y <= -100 & x <= 10393 & y >= -6000'
        block_id = 4
    []
    [./new_block_4]
        type = ParsedSubdomainMeshGenerator
        input = new_block_3
        combinatorial_geometry = '0.5773505384 * x + y > 0 & z > -15000 & x >= 100 & y <= -100 & x <= 10393 & y >= -6000'
        block_id = 5
    []    
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_4
        split_interface = true
        block_pairs = '2 3; 4 5'
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
      [./tria_area_aux]
          order = CONSTANT
          family = MONOMIAL
      []
      [./nodal_area]
          order = CONSTANT
          family = MONOMIAL
      [../]
      #
      [./check_function_initial_stress_xx]
        order = FIRST
        family = LAGRANGE      
      []
      [./check_function_initial_stress_xy]
        order = FIRST
        family = LAGRANGE      
      []
      [./check_function_initial_stress_yy]
        order = FIRST
        family = LAGRANGE      
      []
      [./check_function_initial_stress_zz]
        order = FIRST
        family = LAGRANGE      
      []
[]
    
[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block2_Block3 Block4_Block5'
        strain = SMALL
        generate_output='traction_x traction_y traction_z jump_x jump_y jump_z'
    [../]
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
      [XJump_rate]
        type = FDCompVarRate
        variable = tangent_jump_rate
        coupled = tangent_jump
        execute_on = 'TIMESTEP_BEGIN'
      []
      [XJump_var]
          type = CompVar
          variable = tangent_jump
          coupled = jump_x
          execute_on = 'TIMESTEP_END'
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
      [checkxx]
          type = FunctionAux
          variable = check_function_initial_stress_xx
          function = func_initial_stress_xx
          execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
      [checkxy]
        type = FunctionAux
        variable = check_function_initial_stress_xy
        function = func_initial_strike_shear_stress
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
      [checkyy]
        type = FunctionAux
        variable = check_function_initial_stress_yy
        function = func_initial_stress_yy
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
      [checkzz]
        type = FunctionAux
        variable = check_function_initial_stress_zz
        function = func_initial_stress_zz
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
      [initial_shear_stress_check]
        type = FunctionAux
        variable = ini_shear_stress
        function = func_initial_strike_shear_stress
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
      #fault length
      [fault_len]
        type = ConstantAux
        variable = nodal_area
        value = 200
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
      #SlipWeakeningMultifaults ONLY supports TRIA currently!
      [./czm_mat]
          type = SlipWeakeningMultifaults3D
          disp_slipweakening_x     = disp_slipweakening_x
          disp_slipweakening_y     = disp_slipweakening_y
          disp_slipweakening_z     = disp_slipweakening_z
          reaction_slipweakening_x = resid_slipweakening_x
          reaction_slipweakening_y = resid_slipweakening_y
          reaction_slipweakening_z = resid_slipweakening_z
          reaction_damp_x = resid_damp_sw_x
          reaction_damp_y = resid_damp_sw_y
          reaction_damp_z = resid_damp_sw_z
          nodal_area = nodal_area
          mu_d = mu_d
          mu_s = mu_s
          boundary = 'Block2_Block3 Block4_Block5'
      [../]
      [./static_initial_stress_tensor_slipweakening]
          type = GenericFunctionRankTwoTensor
          tensor_name = static_initial_stress_tensor_slipweakening
          tensor_functions = 'func_initial_stress_xx                func_initial_strike_shear_stress      func_initial_stress_00 
                              func_initial_strike_shear_stress      func_initial_stress_yy                func_initial_stress_00
                              func_initial_stress_00                func_initial_stress_00                func_initial_stress_zz'
      [../]
    []
    
    [Functions]
      [func_static_friction_coeff_mus]
        type = ConstantFunction
        value = 0.677
      []
      #mud constant value
      [func_dynamic_friction_coeff_mud]
          type = ConstantFunction
          value = 0.525
      []
      #Note:restrict stress variation along the fault only
      #this function is used in czm only
      [func_initial_strike_shear_stress]
        type = InitialShearStress
        benchmark_type = tpv14
      []
      [./func_initial_stress_00]
        type = ConstantFunction
        value = 0.0
      []
      [./func_initial_stress_yy]
        type = ConstantFunction
        value = -120e6
      []
      [./func_initial_stress_xx]
        type = ConstantFunction
        value = 0
      []
      [./func_initial_stress_zz]
        type = ConstantFunction
        value = 0
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
    end_time = 1.0
    # num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]
    
#for cluster run
[Outputs]
    exodus = true
    interval = 40
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z traction_x traction_y traction_z jump_x jump_y jump_z mu_s tangent_jump_rate'
[]