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
    [break_boundary]
        type = BreakBoundaryOnSubdomainGenerator
        input = new_block
    []
    [interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = break_boundary
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

[GlobalParams]
    displacements = 'disp_x disp_y'

    ##damping ratio 
    q = 0.5

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
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plus_scaled_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_minus_scaled_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plus_scaled_y]
        order = FIRST
        family = LAGRANGE
    []
    [disp_minus_scaled_y]
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
    #
    [alongfaultdisp_strike_plus]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_strike_minus]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_normal_plus]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_normal_minus]
        order = CONSTANT
        family = MONOMIAL
    []
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
    #
    # [project_disp_plus_x]
    #     type = ProjectionAux
    #     variable = disp_plus_scaled_x
    #     v = alongfaultdisp_strike_plus
    #     boundary = 'Block0_Block1'
    # []
    # [project_disp_minus_x]
    #     type = ProjectionAux
    #     variable = disp_minus_scaled_x
    #     v = alongfaultdisp_strike_minus
    #     boundary = 'Block1_Block0'
    # []
    # [project_disp_plus_y]
    #     type = ProjectionAux
    #     variable = disp_plus_scaled_y
    #     v = alongfaultdisp_normal_plus
    #     boundary = 'Block0_Block1'
    # []
    # [project_disp_minus_y]
    #     type = ProjectionAux
    #     variable = disp_minus_scaled_y
    #     v = alongfaultdisp_normal_minus
    #     boundary = 'Block1_Block0'
    # []
[]

[Problem]
    extra_tag_vectors = 'restore_tag restore_dampx_tag restore_dampy_tag'
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
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
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
[]

[BCs]
    #assign displacement boundary condition
    [./matchval_primary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plus_scaled_x
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_minus_scaled_x
        boundary = 'Block1_Block0'
    []
    [./matchval_primary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_plus_scaled_y
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_minus_scaled_y
        boundary = 'Block1_Block0'
    []
[]

[Executioner]
    type = Transient
    dt = 0.00125
    end_time = 5.0
  #  num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 20
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
        source_variable = 'traction_sub_strike traction_sub_normal sliprate_sub_strike sliprate_sub_normal slip_sub_strike slip_sub_normal statevar_sub alongfaultdisp_strike_plus_sub alongfaultdisp_strike_minus_sub alongfaultdisp_normal_plus_sub alongfaultdisp_normal_minus_sub'
        variable = 'traction_strike traction_normal sliprate_strike sliprate_normal slip_strike slip_normal statevar alongfaultdisp_strike_plus alongfaultdisp_strike_minus alongfaultdisp_normal_plus alongfaultdisp_normal_minus'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #push system residual vector from mainApp to subApp
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'resid_x resid_y resid_damp_x resid_damp_y'
        variable = 'resid_sub_x resid_sub_y resid_damp_sub_x resid_damp_sub_y'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
[]

# [Debug]
#     show_execution_order = ALWAYS
# []