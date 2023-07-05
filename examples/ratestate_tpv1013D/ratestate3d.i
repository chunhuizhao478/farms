#This is the main file for solving for rate-and-state friction

#The following main quantities are defined as material property and are declared in "CZMComputeLocalTractionBaseRSF3D":

#CZMComputeLocalTractionBaseRSF3D <- CZMComputeLocalTractionTotalBaseRSF3D <- RateStateFrictionLaw3Dv4

#traction_strike (traction along strike direction)
#sliprate_strike (sliprate along strike direction)
#slip_strike     (slip along strike direction)
#statevar        (statevar)

#The information stores in the blocks attached to the primary surface

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 3
        nx = 50
        ny = 50
        nz = 50
        xmin = -5000
        xmax = 5000
        ymin = -5000
        ymax = 5000
        zmin = -5000
        zmax = 5000
    [../]
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
    displacements = 'disp_x disp_y disp_z'
    
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
    [./disp_z]
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
    [./resid_rsf_z]
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
    [./resid_z]
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
    [./vel_z]
        order = FIRST
        family = LAGRANGE
    []

    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []

    [./tangent_jump]
        order = CONSTANT
        family = MONOMIAL
    []

    [./tangent_jump_rate]
        order = CONSTANT
        family = MONOMIAL
    []

[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='normal_jump tangent_jump normal_traction tangent_traction'
    [../]
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          displacements = 'disp_x disp_y disp_z'
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
    [Vel_z]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []

    ##pass restoration force
    ##Pass value at time t ##it may not be different using either "BEGIN" or "END"
    ##This operation is not necessary, just to be clear that "resid_rsf" is get updated and feed into material object
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
    [Residual_z_tsend]
        type = CompVar
        variable = resid_rsf_z
        coupled = resid_z
        execute_on = 'TIMESTEP_BEGIN'
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
    [./inertia_z]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_z
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
    [./Reactionz]
        type = StiffPropDamping
        variable = 'disp_z'
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
        type = RateStateFrictionLaw3Dv4
        reaction_rsf_x = resid_rsf_x
        reaction_rsf_y = resid_rsf_y
        reaction_rsf_z = resid_rsf_z
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'Block0_Block1'
        output_properties = 'sliprate_strike slip_strike statevar traction_strike traction_normal traction_dip'
        outputs = exodus
    [../]
[]

[Functions]
    [./func_initial_strike_shear_stress]
        type = InitialStrikeShearStressPerturbRSF3D
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
    interval = 10
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'ratestate3d_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'resid_sub_x resid_sub_y resid_sub_z'
        variable = 'resid_x resid_y resid_z'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'disp_x disp_y disp_z'
        variable = 'disp_sub_x disp_sub_y disp_sub_z'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]