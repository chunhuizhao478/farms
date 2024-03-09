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
    [./accel_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./accel_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./vel_y]
        order = FIRST
        family = LAGRANGE
    [../]
    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        use_displaced_mesh = false
        use_automatic_differentiation = true
        generate_output='normal_jump tangent_jump normal_traction tangent_traction'
    [../]
[]

[Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
    []
[]

# [InterfaceKernels]
#     [tsperturb]
#         type = InitialShearTractionKernel
#         variable = disp_x
#         neighbor_var = disp_x
#         ini_shear_sts_perturb = ini_shear_stress_perturb
#         boundary = 'Block0_Block1'
#     []
# []

[AuxKernels]
    [./accel_x] 
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = TIMESTEP_END
    [../]
    [./vel_x] 
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = TIMESTEP_END
    [../]
    [./accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = TIMESTEP_END
    [../]
    [./vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = TIMESTEP_END
    [../]
    ##initial shear traction
    [StrikeShearStress]
        type = FunctionAux
        variable = ini_shear_stress_perturb
        function = func_initial_strike_shear_stress
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[Materials]
    [Elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ADComputeLinearElasticStress
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2670'
    []
    [./czm_mat]
        type = ADRateStateFrictionLaw2D
        boundary = 'Block0_Block1'
    [../]
[]  

[Executioner]
    type = Transient
    solve_type = 'PJFNK'
    start_time = 0
    end_time = 1.5
    num_steps = 1
    l_tol = 1e-6
    nl_rel_tol = 1e-4
    nl_max_its = 100
    nl_abs_tol = 1e-6
    dt = 0.005
    timestep_tolerance = 1e-6
    petsc_options_iname = '-pc_type -ksp_gmres_restart'
    petsc_options_value = 'lu       101'
    automatic_scaling = true
    nl_forced_its = 2
    line_search = 'none'
[]

[Functions]
    [./func_initial_strike_shear_stress]
        type = InitialStrikeShearStressPerturbRSF2D
    []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]