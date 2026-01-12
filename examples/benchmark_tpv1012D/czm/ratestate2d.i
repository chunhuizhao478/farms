#This is the main file for solving for rate-and-state friction

#The following main quantities are defined as material property and are declared in "CZMComputeLocalTractionBaseRSF2D":

#CZMComputeLocalTractionBaseRSF2D <- CZMComputeLocalTractionTotalBaseRSF2D <- RateStateFrictionLaw2Dv7

#traction_strike (traction along strike direction)
#sliprate_strike (sliprate along strike direction)
#slip_strike     (slip along strike direction)
#statevar        (statevar)

#The information stores in the blocks attached to the primary surface

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 400
        ny = 400
        xmin = -40000
        xmax = 40000
        ymin = -40000
        ymax = 40000
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y>0 & x>-18000 & x<18000'
      block_id = 1
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'y<0 & x>-18000 & x<18000'
      block_id = 2
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block_2
      split_interface = true
      block_pairs = '1 2'
    []
[]

[GlobalParams]
    ##primary variable
    displacements = 'disp_x disp_y'
    
    ##damping ratio 
    q = 0.1
    
    ##element length (m)
    len = 200
    
    ##rate-and-state coefficients
    f_o = 0.6
    rsf_a = 0.008
    rsf_b = 0.012
    rsf_L = 0.02
    delta_o = 1e-6

    ##initial normal traction (Pa)
    Tn_o = 120e6

    ##initial shear traction (Pa)
    Ts_o = 75e6

    ##initial sliprate (m/s)
    Vini = 1e-12

    ##initial state variable
    statevarini = 1.606238999213454e9
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
    ##restoration force (pass into material kernel)
    [./resid_rsf_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_rsf_y]
        order = FIRST
        family = LAGRANGE
    [../]

    [./resid_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    [../]
    
    ##velocity vector (output)
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_y]
        order = FIRST
        family = LAGRANGE
    []

    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []

    #tangent jump
    [./tangent_jump]
        order = CONSTANT
        family = MONOMIAL
    []

    #tangent jump rate
    [./jump_x_rate]
        order = CONSTANT
        family = MONOMIAL
    []

    #slip rate
    [./sliprate_strike]
        order = FIRST
        family = MONOMIAL
    []

[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block1_Block2'
        strain = SMALL
        generate_output = 'traction_x traction_y jump_x jump_y'
    [../]
[]

[Physics]
    [SolidMechanics]
        [QuasiStatic]
            [all]
                strain = SMALL
                add_variables = true
                planar_formulation = PLANE_STRAIN
                generate_output = 'stress_xx stress_yy stress_xy'
                extra_vector_tags = 'restore_tag'
            []
        []
    []
[]

[Problem]
    extra_tag_vectors = 'restore_tag'
[]

[AuxKernels]
    ##Compute velocity (output) 
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

    ##Compute tangent jump rate
    [TJump_rate]
        type = FDCompVarRate
        variable = jump_x_rate
        coupled = jump_x
        execute_on = 'TIMESTEP_END'
    []

    ##pass restoration force
    [Residual_x]
      type = ProjectionAux
      variable = resid_rsf_x
      v = resid_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
      type = ProjectionAux
      variable = resid_rsf_y
      v = resid_y
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

    ##initial shear traction
    [StrikeShearStress]
        type = FunctionAux
        variable = ini_shear_stress_perturb
        function = func_initial_strike_shear_stress
        execute_on = 'TIMESTEP_BEGIN'
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
    [./czm_mat]
        type = RateStateFrictionLaw2Dv7
        reaction_rsf_x = resid_rsf_x
        reaction_rsf_y = resid_rsf_y
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'Block1_Block2'
    [../]
[]

[Functions]
    [./func_initial_strike_shear_stress]
        type = InitialStrikeShearStressPerturbRSF2D
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
    end_time = 3.0
    # num_steps = 5
    [TimeIntegrator]
        type = CentralDifference
        solve_type = consistent
    []
[]

[Outputs]
    exodus = true
    interval = 10
[]

[BCs]
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = right
    []
[]