#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mirrormesh.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '0 -30000 0'
        new_boundary = corner_ptr
        input = sidesets
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
        
    ##----continuum damage breakage model----##
    #initial lambda value (SECOND lame constant) [Pa]
    lambda_o = 10e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 10e9
    
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
    Cd_constant = 300

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 500

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 10

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 3

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.1

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.01 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.7
        
[]

[Variables]
    [disp_x]
        order = SECOND
        family = LAGRANGE
    []
    [disp_y]
        order = SECOND
        family = LAGRANGE
    []
[]

[AuxVariables]
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [correlated_shear_modulus_o]
        order = SECOND
        family = LAGRANGE
    []
[]

[AuxKernels]
    [get_initial_damage]
        type = MaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
    []
    [get_shear_modulus_o]
        type = FunctionAux
        variable = correlated_shear_modulus_o
        function = node_shear_modulus
        execute_on = 'INITIAL'
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
        output_properties = 'alpha_damagedvar B_damagedvar shear_modulus_o_mat shear_modulus'
        outputs = exodus
        block = '1 2'
        use_shear_modulus_o_aux = true
        shear_modulus_o_aux = correlated_shear_modulus_o
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain deviatroic_stress'
        outputs = exodus
        block = '1 2'
    []
    # elastic
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
        block = 3
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
        block = 3
    []
    [initial_damage_surround]
        type = InitialDamageCycleSim2DDebug
        len_of_fault = 1000
        sigma = 5e2
        peak_val = 0.7
        output_properties = 'initial_damage'
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
    type = Steady
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    # start_time = -1e-12
    # end_time = 1e100
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 20
    nl_abs_tol = 1e-10
    petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    # dt = 1e-8
    # steady_state_detection = true
[]  

[Outputs]
    exodus = true       
[]

[BCs]
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    [] 
    [./Pressure]
        [static_pressure_top]
            boundary = top
            factor = 120e6
            displacements = 'disp_x disp_y'
        []    
        [static_pressure_left]
            boundary = left
            factor = 135e6
            displacements = 'disp_x disp_y'
        []  
        [static_pressure_right]
            boundary = right
            factor = 135e6
            displacements = 'disp_x disp_y'
        []     
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
    #add initial shear stress
    [./initial_shear_stress]
        type = NeumannBC
        variable = disp_x
        value = 20e6
        boundary = top
    []
[]

[UserObjects]
    [reader_node_shear_modulus]
        type = PropertyReadFile
        prop_file_name = 'mapped_vonkarman_field.csv'
        read_type = 'node'
        nprop = 3 # number of columns in CSV
    []
[]

[Functions]
    [node_shear_modulus]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'reader_node_shear_modulus'
        read_type = 'node'
        # 0-based indexing
        column_number = '2'
    []
[]