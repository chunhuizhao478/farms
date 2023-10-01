[Mesh]
    [generated1]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 201
        ny = 200
        xmin = -10050
        xmax = 10050
        ymin = -10000
        ymax = 10000
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = generated1
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    [interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = new_block
        primary_block = 0
        paired_block = 1
        new_boundary = 'Block0_Block1'
    []
    [interface2]
        type = SideSetsBetweenSubdomainsGenerator
        input = interface
        primary_block = 1
        paired_block = 0
        new_boundary = 'Block1_Block0'
    []
[]

[Variables]
    #block0
    [./disp_x_b0]
        order = FIRST
        family = LAGRANGE
        block = 0
    [../]
    [./disp_y_b0]
        order = FIRST
        family = LAGRANGE
        block = 0
    [../]
    #block1
    [./disp_x_b1]
        order = FIRST
        family = LAGRANGE
        block = 1
    [../]
    [./disp_y_b1]
        order = FIRST
        family = LAGRANGE
        block = 1
    [../]
[]

[GlobalParams]
    q = 0.1
[]

[AuxVariables]
    #restoration force (tag after solve)
    [./resid_x_b0]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y_b0]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_x_b1]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y_b1]
        order = FIRST
        family = LAGRANGE
    [../]
    #restoration force for damping (tag after solve)
    [./resid_damp_x_b0]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damp_y_b0]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_x_b1]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damp_y_b1]
        order = FIRST
        family = LAGRANGE
    [../]
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plus_scaled_x]
        order = FIRST
        family = LAGRANGE
        initial_condition = 1
    []
    [disp_minus_scaled_x]
        order = FIRST
        family = LAGRANGE
        initial_condition = -1
    []
    [disp_plus_scaled_y]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
    []
    [disp_minus_scaled_y]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
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
    #velocity
    [vel_x_b0]
        order = FIRST
        family = LAGRANGE
    []
    [vel_x_b1]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y_b0]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y_b1]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    #obtain system residuals by tagging
    #sts
    #block 0
    [restore_stsx_b0]
        type = TagVectorAux
        vector_tag = 'restore_stsx_b0_tag'
        v = 'disp_x_b0'
        variable = 'resid_x_b0'
        execute_on = 'TIMESTEP_END'
        block = 0
    []
    [restore_stsy_b0]
        type = TagVectorAux
        vector_tag = 'restore_stsy_b0_tag'
        v = 'disp_y_b0'
        variable = 'resid_y_b0'
        execute_on = 'TIMESTEP_END'
        block = 0
    []
    #block 1
    [restore_stsx_b1]
        type = TagVectorAux
        vector_tag = 'restore_stsx_b1_tag'
        v = 'disp_x_b1'
        variable = 'resid_x_b1'
        execute_on = 'TIMESTEP_END'
        block = 1
    []
    [restore_stsy_b1]
        type = TagVectorAux
        vector_tag = 'restore_stsy_b1_tag'
        v = 'disp_y_b1'
        variable = 'resid_y_b1'
        execute_on = 'TIMESTEP_END'
        block = 1
    []
    #damping
    #block 0
    [restore_dampx_b0]
        type = TagVectorAux
        vector_tag = 'restore_dampx_b0_tag'
        v = 'disp_x_b0'
        variable = 'resid_damp_x_b0'
        execute_on = 'TIMESTEP_END'
        block = 0
    []
    [restore_dampy_b0]
        type = TagVectorAux
        vector_tag = 'restore_dampy_b0_tag'
        v = 'disp_y_b0'
        variable = 'resid_damp_y_b0'
        execute_on = 'TIMESTEP_END'
        block = 0
    []
    #damping
    #block 1
    [restore_dampx_b1]
        type = TagVectorAux
        vector_tag = 'restore_dampx_b1_tag'
        v = 'disp_x_b1'
        variable = 'resid_damp_x_b1'
        execute_on = 'TIMESTEP_END'
        block = 1
    []
    [restore_dampy_b1]
        type = TagVectorAux
        vector_tag = 'restore_dampy_b1_tag'
        v = 'disp_y_b1'
        variable = 'resid_damp_y_b1'
        execute_on = 'TIMESTEP_END'
        block = 1
    []
    #calc velocity
    #block 0
    [Vel_x_b0]
        type = CompVarRate
        variable = vel_x_b0
        coupled = disp_x_b0
        execute_on = 'TIMESTEP_BEGIN'
        block = 0
    []
    [Vel_y_b0]
        type = CompVarRate
        variable = vel_y_b0
        coupled = disp_y_b0
        execute_on = 'TIMESTEP_BEGIN'
        block = 0
    []
    #block 1
    [Vel_x_b1]
        type = CompVarRate
        variable = vel_x_b1
        coupled = disp_x_b1
        execute_on = 'TIMESTEP_BEGIN'
        block = 1
    []
    [Vel_y_b1]
        type = CompVarRate
        variable = vel_y_b1
        coupled = disp_y_b1
        execute_on = 'TIMESTEP_BEGIN'
        block = 1
    []
[]

[Kernels]
    #block 0
    [./disp_b0_x]
        type = StressDivergenceTensors
        variable = disp_x_b0
        displacements = 'disp_x_b0 disp_y_b0'
        component = 0
        extra_vector_tags = 'restore_stsx_b0_tag'
        block = 0
    []
    [./disp_b0_y]
        type = StressDivergenceTensors
        variable = disp_y_b0
        displacements = 'disp_x_b0 disp_y_b0'
        component = 1
        extra_vector_tags = 'restore_stsy_b0_tag'
        block = 0
    []
    #block 1
    [./disp_b1_x]
        type = StressDivergenceTensors
        variable = disp_x_b1
        displacements = 'disp_x_b1 disp_y_b1'
        component = 0
        extra_vector_tags = 'restore_stsx_b1_tag'
        block = 1
    []
    [./disp_b1_y]
        type = StressDivergenceTensors
        variable = disp_y_b1
        displacements = 'disp_x_b1 disp_y_b1'
        component = 1
        extra_vector_tags = 'restore_stsy_b1_tag'
        block = 1
    []
    #block 0
    [./inertia_x_b0]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x_b0
        block = 0
    []
    [./inertia_y_b0]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y_b0
        block = 0
    []
    #block 1
    [./inertia_x_b1]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x_b1
        block = 1
    []
    [./inertia_y_b1]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y_b1
        block = 1
    []
    #block 0
    [./Reactionx_b0]
        type = StiffPropDamping
        displacements = 'disp_x_b0 disp_y_b0'
        variable = disp_x_b0
        component = 0
        extra_vector_tags = 'restore_dampx_b0_tag'
        block = 0
    []
    [./Reactiony_b0]
        type = StiffPropDamping
        displacements = 'disp_x_b0 disp_y_b0'
        variable = disp_y_b0
        component = 1
        extra_vector_tags = 'restore_dampy_b0_tag'
        block = 0
    []
    #block 1
    [./Reactionx_b1]
        type = StiffPropDamping
        displacements = 'disp_x_b1 disp_y_b1'
        variable = disp_x_b1
        component = 0
        extra_vector_tags = 'restore_dampx_b1_tag'
        block = 1
    []
    [./Reactiony_b1]
        type = StiffPropDamping
        displacements = 'disp_x_b1 disp_y_b1'
        variable = disp_y_b1
        component = 1
        extra_vector_tags = 'restore_dampy_b1_tag'
        block = 1
    []
[]

[Materials]
    [elasticity_b0]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
        block = 0
    []
    [elasticity_b1]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
        block = 1
    []
    [stress_b0]
        type = ComputeLinearElasticStress
        block = 0
    []
    [stress_b1]
        type = ComputeLinearElasticStress
        block = 1
    []
    [strain_b0]
        type = ComputeSmallStrain
        displacements = 'disp_x_b0 disp_y_b0'
        block = 0
    []
    [strain_b1]
        type = ComputeSmallStrain
        displacements = 'disp_x_b1 disp_y_b1'
        block = 1
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
[]

[BCs]
    #assign displacement boundary condition
    #block0_block1
    [./matchval_primary_xt]
        type = MatchedValueBC
        variable = disp_x_b0
        v = disp_plus_scaled_x
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_xb]
        type = MatchedValueBC
        variable = disp_x_b1
        v = disp_minus_scaled_x
        boundary = 'Block1_Block0'
    []
    #block1_block0
    [./matchval_primary_yt]
        type = MatchedValueBC
        variable = disp_y_b0
        v = disp_plus_scaled_y
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_yb]
        type = MatchedValueBC
        variable = disp_y_b1
        v = disp_minus_scaled_y
        boundary = 'Block1_Block0'
    []
[]

[UserObjects]
    #evalute system residual after system solve before auxkernels
    #sts
    #block 0
    [recompute_residual_tag_stsx_b0]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_stsx_b0_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_stsy_b0]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_stsy_b0_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    #block 1
    [recompute_residual_tag_stsx_b1]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_stsx_b1_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_stsy_b1]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_stsy_b1_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    #damping
    #block 0 
    [recompute_residual_tag_dampx_b0]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_b0_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampy_b0]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampy_b0_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    #block 1 
    [recompute_residual_tag_dampx_b1]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_b1_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampy_b1]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampy_b1_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
[]

[Executioner]
    type = Transient
    dt = 1
    end_time = 5.0
    num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]

[Problem]
    extra_tag_vectors = 'restore_stsx_b0_tag restore_stsy_b0_tag restore_stsx_b1_tag restore_stsy_b1_tag restore_dampx_b0_tag restore_dampy_b0_tag restore_dampx_b1_tag restore_dampy_b1_tag'
[]

# [MultiApps]
#     #allocate transfer from mainApp to subApp
#     [./sub_app]
#       type = TransientMultiApp
#       positions = '0 0 0'
#       input_files = 'sub.i'
#       execute_on = 'INITIAL TIMESTEP_BEGIN'
#     [../]
# []

# [Transfers]
#     #get displacement residuals from subApp to mainApp
#     [pull_resid]
#         type = MultiAppCopyTransfer
#         from_multi_app = sub_app
#         source_variable = 'disp_plusminus_sub_scaled_x disp_plusminus_sub_scaled_y traction_sub_strike traction_sub_normal sliprate_sub_strike sliprate_sub_normal slip_sub_strike slip_sub_normal statevar_sub'
#         variable = 'disp_plusminus_scaled_x disp_plusminus_scaled_y traction_strike traction_normal sliprate_strike sliprate_normal slip_strike slip_normal statevar'
#         execute_on = 'INITIAL TIMESTEP_BEGIN'
#     []
#     #push system residual vector from mainApp to subApp
#     [push_disp]
#         type = MultiAppCopyTransfer
#         to_multi_app = sub_app
#         source_variable = 'resid_x resid_y resid_damp_x resid_damp_y'
#         variable = 'resid_sub_x resid_sub_y resid_damp_sub_x resid_damp_sub_y'
#         execute_on = 'INITIAL TIMESTEP_BEGIN'
#     []
# []
