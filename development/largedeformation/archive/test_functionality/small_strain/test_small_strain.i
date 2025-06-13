#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/tabulardamagezone_small_refine.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0
                    0 0 -1
                    0 0 1'
        new_boundary = 'left right bottom top back front'
    [] 
[]

[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
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
    [disp_z]
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
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        static_initialization = true
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        static_initialization = true
    []
    [dispkernel_z]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        static_initialization = true
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_z]
        type = ADInertialForce
        variable = disp_z
        velocity = vel_z
        acceleration = accel_z
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
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
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = timestep_end
    []
[]

[Materials]
    [damagestress]
        type = ADComputeLinearElasticStress
        outputs = exodus
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [./elasticity_volume_1]
        type = ADComputeIsotropicElasticityTensor
        shear_modulus = 7.7376e9
        lambda = 30e9
    [../]   
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
      petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '      
    []
  []
  
[Executioner]
    type = Transient
    solve_type = Newton
    start_time = 0
    # end_time = 800
    num_steps = 1000
    l_max_its = 10
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 10
    nl_abs_tol = 1e-6
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1
    [TimeIntegrator]
        type = NewmarkBeta
    []
[]  

[Outputs]
    exodus = true
    interval = 10
    # show = 'disp_x disp_y disp_z'
    print_linear_residuals=true
[]

#We assume the simulation is loaded with compressive pressure and shear stress
[BCs]
    [fix_front_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = front
    []
    [fix_front_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = front
    [] 
    [pressure_front]
        type = ADFunctionDirichletBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        function = '2.0 / 1000 * t'
    []   
    [fix_back_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = back
    []
    [fix_back_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = back
    []
    [fix_back_z]
        type = ADDirichletBC
        variable = disp_z
        value = 0
        boundary = back
    []
[]