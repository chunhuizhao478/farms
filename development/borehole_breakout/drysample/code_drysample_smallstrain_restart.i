[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh.msh'
    [] 
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 19.9e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 15.6e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.5
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.9
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 2e3

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.05 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
    
    #
    D = 0
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
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
[]

[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
[]

[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
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
    # [./inertia_z]
    #     type = InertialForce
    #     use_displaced_mesh = false
    #     variable = disp_z
    #     acceleration = accel_z
    #     velocity = vel_z
    #     beta = 0.25
    #     gamma = 0.5
    #     eta = 0
    # [] 
[]

[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        # outputs = exodus
    [] 
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2640'
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBM
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 xi stress'
        block = '3'
        outputs = exodus
    [] 
    [stress_elastic]
        type = ComputeLinearElasticStress
        block = '1 2'
        output_properties = 'elastic_strain stress'
        outputs = exodus
    []
    [elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 48.5e9
        poissons_ratio = 0.22
    []
    [dummy_matprop]
        type = GenericConstantMaterial
        prop_names = 'initial_damage initial_breakage shear_stress_perturbation'
        prop_values = '0.0 0.0 0.0'  
    []    
[]  

#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-2.6477e-5 - 3.3e-7 * t'
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
    # solve_type = 'NEWTON'
    solve_type = 'PJFNK'
    start_time = 1550
    end_time = 4000 #extend the time
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 20
    nl_abs_tol = 1e-8
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1e-1
    [./TimeIntegrator]
        type = ImplicitEuler
        # type = BDF2
        # type = CrankNicolson
    [../]
    # [TimeStepper]
    #     type = FarmsIterationAdaptiveDT
    #     dt = 0.1
    #     cutback_factor_at_failure = 0.5
    #     optimal_iterations = 8
    #     growth_factor = 1.5
    #     max_time_step_bound = 10
    # []
[]

[Outputs] 
    exodus = true
    time_step_interval = 50
    show = 'stress_22 B alpha_damagedvar xi eps_e_22 vel_x vel_y vel_z'
    [./csv]
        type = CSV
        time_step_interval = 1
        show = 'strain_z react_z'
    [../]
    [out]
        type = Checkpoint
        time_step_interval = 50
        num_files = 2
    []
[]

[BCs]
    #fix bottom boundary
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 7
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 7
        value = 0
    []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = 7
        value = 0
    []
    #applied load on top boundary
    [applied_top_z_dispload]
        type = FunctionDirichletBC
        variable = disp_z
        boundary = 6
        function = applied_load_top
    [] 
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = 4
          factor = 17.2e6
          displacements = 'disp_x disp_y'
        [../]
    []
[]

#compute the reaction force on the top boundary
[Postprocessors]
    [./react_z]
      type = SidesetReaction
      direction = '0 0 1'
      stress_tensor = stress
      boundary = 6
    [../]
    [./strain_z]
        type = FunctionValuePostprocessor
        function = applied_load_top
    []
[]

[Problem]
    #Note that the suffix is left off in the parameter below.
    restart_file_base = code_drysample_smallstrain_out_cp/LATEST  # You may also use a specific number here
[]