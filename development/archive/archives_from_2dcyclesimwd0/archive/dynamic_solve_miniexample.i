[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 100
        ny = 100
        xmin = -60000
        xmax = 60000
        ymin = -60000
        ymax = 60000
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '0 -60000 0'
        new_boundary = corner_ptr
        input = msh
    []
    displacements = 'disp_x disp_y'
    use_displaced_mesh = true
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
[]


[Variables]
    [disp_x]
        order = FIRST
        family = LAGRANGE     
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE    
    []
[]

[AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = timestep_end
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = timestep_end
    []
[]

[Kernels]
    [dispkernel_x]
        type = TotalLagrangianStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
    [dispkernel_y]
        type = TotalLagrangianStressDivergence
        variable = disp_y
        component = 1
        large_kinematics = true
    []
    [./inertia_x]
        type = InertialForce
        use_displaced_mesh = true
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = true
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
    []   
[]  

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    start_time = -1e-8
    end_time = 1e100
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 20
    nl_abs_tol = 1e-10
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'bt'
    dt = 1e-8
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
        inactive_tsteps = 1
    [../]
[]

[Controls] # turns off inertial terms for the first time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */accel_x */accel_y */inertia_x */inertia_y'
    start_time = 0.0
    end_time = 2e-8 # dt used in the simulation
  [../]
[../]

[Outputs] 
    exodus = true
    time_step_interval = 1
[]


[Functions]
    [func_top_bc]
        type = ParsedFunction
        expression = 'if (t>dt, 1e-9 * t, 0)'
        symbol_names = 'dt'
        symbol_values = '2e-8'
    []
[]

[BCs]
    [bc_load_top_x]
        type = FunctionDirichletBC
        variable = disp_x
        function = func_top_bc
        boundary = top
    []    
    [bc_load_top_y]
        type = DirichletBC
        variable = disp_y
        value = -337.5
        boundary = top
    []
    [bc_fix_bottom_x]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    [bc_load_left_x]
        type = DirichletBC
        variable = disp_x
        value = 213.75
        boundary = left
    []
    [bc_load_right_x]
        type = DirichletBC
        variable = disp_x
        value = -213.75
        boundary = right
    []   
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr2_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
[]