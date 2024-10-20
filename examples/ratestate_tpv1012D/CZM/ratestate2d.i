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
        nx = 300
        ny = 300
        xmin = -15000
        xmax = 15000
        ymin = -15000
        ymax = 15000
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
    []
[]

[GlobalParams]
    ##primary variable
    displacements = 'disp_x disp_y'
    
    ##damping ratio 
    q = 0.1
    
    ##element length (m)
    len = 100
    
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
    [./tangent_jump_rate]
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
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output = 'tangent_jump'
    [../]
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          planar_formulation = PLANE_STRAIN
          displacements = 'disp_x disp_y'
        [../]
      [../]
    [../]
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
        variable = tangent_jump_rate
        coupled = tangent_jump
        execute_on = 'TIMESTEP_END'
    []

    ##pass restoration force
    ##Pass value at time t ##it may not be different using either "BEGIN" or "END"
    [Residual_x_tsend]
        type = CompVar
        variable = resid_rsf_x
        coupled = resid_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y_tsend]
        type = CompVar
        variable = resid_rsf_y
        coupled = resid_y
        execute_on = 'TIMESTEP_BEGIN'
    []

    ##initial shear traction
    [StrikeShearStress]
        type = FunctionAux
        variable = ini_shear_stress_perturb
        function = func_initial_strike_shear_stress
        execute_on = 'TIMESTEP_BEGIN'
    []

    ##output
    [outputsliprate]
        type = MaterialRealAux
        property = sliprate_strike
        variable = sliprate_strike
        boundary = 'Block0_Block1'
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
        boundary = 'Block0_Block1'
    [../]
[]

[Functions]
    [./func_initial_strike_shear_stress]
        type = InitialStrikeShearStressPerturbRSF2D
    []
[]

[Executioner]
    type = Transient
    dt = 0.0025
    end_time = 3.0
    # num_steps = 5
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 10
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'ratestate2d_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'resid_sub_x resid_sub_y'
        variable = 'resid_x resid_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'disp_x disp_y'
        variable = 'disp_sub_x disp_sub_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]