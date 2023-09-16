[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 800
        ny = 800
        xmin = -10000
        xmax = 10000
        ymin = -10000
        ymax = 10000
    [../]
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    #add "Block0_Block1" and "Block1_Block0" interfaces
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
    q = 0.2

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
     #restoration force damping (tag after solve)
     [./resid_damp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damp_y]
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
    #sliprate
    [sliprate_mag]
        order = CONSTANT
        family = MONOMIAL
    []
    #slip
    [slip_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    #statevar
    [statevar]
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
    #damping
    [restore_damp_x]
        type = TagVectorAux
        vector_tag = 'restore_dampx_tag'
        v = 'disp_x'
        variable = 'resid_damp_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_damp_y]
        type = TagVectorAux
        vector_tag = 'restore_dampy_tag'
        v = 'disp_y'
        variable = 'resid_damp_y'
        execute_on = 'TIMESTEP_END'
    []
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
    [recompute_dampx_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_dampy_residual_tag]
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
[]

[Executioner]
    type = Transient
    dt = 0.00125
    end_time = 3.0
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
      sub_cycling = true
    [../]
  []

[Transfers]
    #get displacement residuals from subApp to mainApp
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'disp_plusminus_sub_scaled_x disp_plusminus_sub_scaled_y traction_sub_strike sliprate_sub_mag slip_sub_strike statevar_sub'
        variable = 'disp_plusminus_scaled_x disp_plusminus_scaled_y traction_strike sliprate_mag slip_strike statevar'
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