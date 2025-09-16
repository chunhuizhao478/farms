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
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = ' -15000 -15000 0'
        new_boundary = corner_ptr
        input = split
    []    
[]

[GlobalParams]
    ##primary variable
    displacements = 'disp_x disp_y'
    
    # ##rate-and-state coefficients
    # fo = 0.6
    # a = 0.008
    # b = 0.012
    # length_scale_ref = 0.02
    # slip_rate_ref = 1e-6

    # ##initial normal traction (Pa)
    # T2_o = -120e6

    # ##initial shear traction (Pa)
    # T1_o = 75e6

    # ##initial sliprate (m/s)
    # slip_rate_ini = 1e-12

    # ##initial state variable
    # state_variable_ini = 1.606238999213454e9
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
    ##velocity vector (output)
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
    []
    [accel_y]
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

[]

# [Modules/TensorMechanics/CohesiveZoneMaster]
#     [./czm_ik]
#         boundary = 'Block0_Block1'
#         strain = SMALL
#         generate_output = 'traction_x traction_y jump_x jump_y'
#     [../]
# []

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

# [AuxKernels]
#     [accel_x]
#         type = NewmarkAccelAux
#         variable = accel_x
#         displacement = disp_x
#         velocity = vel_x
#         beta = 0.25
#         execute_on = 'TIMESTEP_END'
#     []
#     [vel_x]
#         type = NewmarkVelAux
#         variable = vel_x
#         acceleration = accel_x
#         gamma = 0.5
#         execute_on = 'TIMESTEP_END'
#     []
#     [accel_y]
#         type = NewmarkAccelAux
#         variable = accel_y
#         displacement = disp_y
#         velocity = vel_y
#         beta = 0.25
#         execute_on = 'TIMESTEP_END'
#     []
#     [vel_y]
#         type = NewmarkVelAux
#         variable = vel_y
#         acceleration = accel_y
#         gamma = 0.5
#         execute_on = 'TIMESTEP_END'
#     []
#     ##initial shear traction
#     [StrikeShearStress]
#         type = FunctionAux
#         variable = ini_shear_stress_perturb
#         function = func_initial_strike_shear_stress
#         execute_on = 'TIMESTEP_BEGIN'
#     []
# []

[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
    []
    # [./inertia_x]
    #     type = InertialForce
    #     use_displaced_mesh = false
    #     variable = disp_x
    #     acceleration = accel_x
    #     velocity = vel_x
    #     beta = 0.25
    #     gamma = 0.5
    #     eta = 0
    # []
    # [./inertia_y]
    #     type = InertialForce
    #     use_displaced_mesh = false
    #     variable = disp_y
    #     acceleration = accel_y
    #     velocity = vel_y
    #     beta = 0.25
    #     gamma = 0.5
    #     eta = 0
    # []       
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
    # [density]
    #     type = GenericConstantMaterial
    #     prop_names = density
    #     prop_values = 2670
    # []
    # [./czm_mat]
    #     type = RateStateFriction2d
    #     T1_perturb = ini_shear_stress_perturb
    #     boundary = 'Block0_Block1'
    # [../]
[]

# [Functions]
#     [./func_initial_strike_shear_stress]
#         type = InitialStrikeShearStressPerturbRSF2D
#     []
# []

[BCs]
    [compression_top_y]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -120e6
    []
    [compression_bottom_y]
        type = NeumannBC
        variable = disp_y
        boundary = bottom
        value = -120e6
    []
    [shear_top_x]
        type = NeumannBC
        variable = disp_x
        boundary = top
        value = 75e6
    []
    [shear_bottom_x]
        type = NeumannBC
        variable = disp_x
        boundary = bottom
        value = -75e6
    []
    [fix_ptr_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [fix_ptr_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
    type = Steady
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    # start_time = 0
    # end_time = 12
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 8
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    # automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    # dt = 0.005
    # [./TimeIntegrator]
    #     type = NewmarkBeta
    #     beta = 0.25
    #     gamma = 0.5
    # [../]
[]

[Outputs]
    exodus = true
    interval = 1
[]