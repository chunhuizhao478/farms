##########################################################
# TPV24 benchmark setup 
##########################################################

[Mesh]
    [./msh]
      type = FileMeshGenerator
      file =  '../../meshgenerator/cdbm/planarfault3D/TPV24.msh'
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x > -16000 & x < 12000 & y < 0 & z > -15000'
      block_id = 1
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x > -16000 & x < 12000 & y > 0 & z > -15000'
        block_id = 2
    []    
    [./new_block_3]
        type = ParsedSubdomainMeshGenerator
        input = new_block_2
        combinatorial_geometry = '0.5773505384 * x + y < 0 & z > -15000 & x >= 150 & y <= -150 & x <= 10417 & y >= -6025'
        block_id = 3
    []
    [./new_block_4]
        type = ParsedSubdomainMeshGenerator
        input = new_block_3
        combinatorial_geometry = '0.5773505384 * x + y > 0 & z > -15000 & x >= 150 & y <= -150 & x <= 10417 & y >= -6025'
        block_id = 4
    []    
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_4
        split_interface = true
        block_pairs = '1 2; 3 4'
    []     
[]
    
[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z'
    
    #damping ratio
    q = 0.2
    
    #characteristic length (m)
    Dc = 0.3
    
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
      [./tangent_jump_rate]
          order = CONSTANT
          family = MONOMIAL
      []
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
    
[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block1_Block2 Block3_Block4'
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
    extra_tag_vectors = 'restore_tag'
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
        value = 150
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
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
      [./Reactionx]
        type = StiffPropDamping
        variable = 'disp_x'
        component = '0'
      []
      [./Reactiony]
        type = StiffPropDamping
        variable = 'disp_y'
        component = '1'
      []
    []
    
    [Materials]
      #damage breakage model
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
          nodal_area = nodal_area
          mu_d = mu_d
          mu_s = mu_s
          tria_area = tria_area_aux
          cohesion = cohesion_aux
          forced_rupture_time = forced_rupture_time_aux
          boundary = 'Block1_Block2 Block3_Block4'
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
          value = 0.18
      []
      #mud constant value: 0.4
      [func_dynamic_friction_coeff_mud]
          type = ConstantFunction
          value = 0.12
      []
      #Note:restrict stress variation along the fault only
      #this function is used in czm only
      [func_initial_strike_shear_stress]
        type = InitialStresscontmfbfs3D
        i = 1
        j = 2
      []
      #this function is used in medimum
      [func_initial_stress_xy_const]
        type = InitialStresscontmfbfs3D
        i = 1
        j = 2
      []
      [./func_initial_stress_00]
        type = ConstantFunction
        value = 0.0
      []
      [./func_initial_stress_yy]
        type = InitialStresscontmfbfs3D
        i = 2
        j = 2
      []
      [./func_initial_stress_xx]
        type = InitialStresscontmfbfs3D
        i = 1
        j = 1
      []
      [./func_initial_stress_zz]
        type = InitialStresscontmfbfs3D
        i = 3
        j = 3
      []
      [./func_initial_cohesion]
        type = InitialCohesion3D
      []
      [./func_forced_rupture_time]
        type = ForcedRuptureTime
        loc_x = -8000
        loc_y = 0
        loc_z = -10000
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
[]
    
[Executioner]
    type = Transient
    dt = 0.0025
    end_time = 12.0
    # num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]
    
#for cluster run
[Outputs]
    exodus = true
    interval = 40
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z traction_x traction_y traction_z jump_x jump_y jump_z'
[]