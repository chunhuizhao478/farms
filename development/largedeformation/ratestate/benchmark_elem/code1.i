[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 3
        nx = 1
        ny = 1
        nz = 1
        xmin = 0
        xmax = 1
        ymin = 0
        ymax = 1
        zmin = 0
        zmax = 1
    [] 
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
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
    Cd_constant = 300 #

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 500

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.01 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 5e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.7
    
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
    [xi_computed]
        order = CONSTANT
        family = MONOMIAL
    []
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
    [shear_strain]
        order = CONSTANT
        family = MONOMIAL
    []
    [shear_strain_rate]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [state_variable_rate]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [compute_xi]
        type = CompXi3D
        variable = xi_computed
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    #
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    [Vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    #
    [shear_strain]
        type = MaterialRankTwoTensorAux
        property = total_lagrange_strain
        i = 0
        j = 1
        variable = shear_strain
        block = 0
    []
    [shear_strain_rate]
        type = TimeDerivativeAux
        variable = shear_strain_rate
        functor = shear_strain
        block = 0
    []
    #
    [state_variable_rate]
        type = MaterialRateRealAux
        variable = state_variable_rate
        property = state_variable
        block = 0
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
    [dispkernel_z]
        type = TotalLagrangianStressDivergence
        variable = disp_z
        component = 2
        large_kinematics = true
    []  
[]

[Materials]
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar velgrad_L shear_modulus a0 a1 a2 a3' 
        use_state_var_evolution = true
        const_A = 2e9
        const_B = 1.2e9
        const_theta_o = 1e4
        initial_theta0 = 1e4
        outputs = exodus
        block = 0
    [] 
    [initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = 0
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain state_variable state_variable_tensor strain_invariant_ratio'
        outputs = exodus
    []
[]  

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '1e-4 * t'
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
    start_time = 0
    end_time = 300000 #extend the time
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 5
    nl_abs_tol = 1e-10
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type'
    # petsc_options_value = 'gmres     hypre'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1
    [./TimeIntegrator]
        type = ImplicitEuler
        # type = BDF2
        # type = CrankNicolson
    [../]
[]

[Outputs] 
    exodus = true
    # show = 'pk2_stress_01 total_lagrange_strain_01 state_variable shear_modulus shear_strain_rate state_variable_rate pk2_stress_11'
    time_step_interval = 1
[]

[BCs]
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = bottom
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = bottom
        value = 0
    []
    [fix_top_y]
        type = DirichletBC
        variable = disp_y
        boundary = top
        value = 0
    []
    [applied_top_x]
        type = FunctionDirichletBC
        variable = disp_x
        boundary = top
        function = applied_load_top
    [] 
    [fix_back_z]
        type = DirichletBC
        variable = disp_z
        boundary = back
        value = 0
    []
    [./Pressure]
        [static_pressure_left]
            boundary = left
            factor = 10e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []
        [static_pressure_right]
            boundary = right
            factor = 10e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []  
        [static_pressure_front]
            boundary = front
            factor = 10e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []          
    []
[]