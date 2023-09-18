[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 3
        nx = 50
        ny = 50
        nz = 75
        xmin = -5000
        xmax = 5000
        ymin = -5000
        ymax = 5000
        zmin = -15000
        zmax = 0
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
    displacements = 'disp_x disp_y disp_z'

    ##damping ratio 
    q = 0.1

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
    [./disp_z]
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
    [./resid_z]
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
    [disp_plusminus_scaled_z]
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
    [traction_dip]
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
    [sliprate_dip]
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
    #velocity
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    #
    [tangent_jump]
        order = CONSTANT
        family = MONOMIAL
    []
    [tangent_jump_rate]
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
    [restore_z]
        type = TagVectorAux
        vector_tag = 'restore_tag'
        v = 'disp_z'
        variable = 'resid_z'
        execute_on = 'TIMESTEP_END'
    []
    #calc velocity
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_BEGIN'
    []
    #get jump from uo - aux
    [tangent_jump_rate_aux]
        type = InterfaceValueUserObjectAux
        interface_uo_name = tangent_jump_rate_uo
        variable = tangent_jump_rate
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_x disp_y disp_z'
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
    [./inertia_z]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_z
    []
    [./Reactionx]
        type = StiffPropDamping
        variable = disp_x
        component = 0
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = disp_y
        component = 1
    []
    [./Reactionz]
        type = StiffPropDamping
        variable = disp_z
        component = 2
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
    [tangent_jump_rate_uo]
        type = InterfaceQpValueUserObject
        var = disp_x
        var_neighbor = disp_x
        boundary = 'Block0_Block1'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        interface_value_type = jump_primary_minus_secondary
        value_type = rate
    []
[]

[Problem]
    extra_tag_vectors = 'restore_tag'
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
    [./matchval_primary_z]
        type = MatchedValueBC
        variable = disp_z
        v = disp_plusminus_scaled_z
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_z]
        type = MatchedValueBC
        variable = disp_z
        v = disp_plusminus_scaled_z
        boundary = 'Block1_Block0'
    []
[]

[Executioner]
    type = Transient
    dt = 0.0025
    end_time = 1.5
    # num_steps = 10
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
      input_files = 'sub_200m.i'
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
  []

[Transfers]
    #get displacement residuals from subApp to mainApp
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'disp_plusminus_sub_scaled_x disp_plusminus_sub_scaled_y disp_plusminus_sub_scaled_z traction_sub_strike traction_sub_normal traction_sub_dip sliprate_sub_strike sliprate_sub_normal sliprate_sub_dip slip_sub_strike statevar_sub'
        variable = 'disp_plusminus_scaled_x disp_plusminus_scaled_y disp_plusminus_scaled_z traction_strike traction_normal traction_dip sliprate_strike sliprate_normal sliprate_dip slip_strike statevar'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #push system residual vector from mainApp to subApp
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'resid_x resid_y resid_z'
        variable = 'resid_sub_x resid_sub_y resid_sub_z'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
[]
