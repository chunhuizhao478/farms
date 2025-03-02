#continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_large.msh'
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
        new_boundary = 'left right front back bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = ' -12000 -10000 -20000;
                   12000 -10000 -20000;
                   12000 10000  -20000;
                  -12000 10000  -20000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
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
    xi_d = -0.8
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    # energy ratio
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
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [correlated_randalpha_o]
        order = FIRST
        family = LAGRANGE
    []
    [initial_cd_aux]
        order = FIRST
        family = MONOMIAL
    []
[]

[AuxKernels]
    [get_initial_damage]
        type = MaterialRealAux
        variable = initial_damage_aux
        property = initial_damage_time_dependent_material
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
    [gravity]
        type = Gravity
        variable = disp_z
        value = -9.81
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
    # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar shear_modulus_o_mat shear_modulus'
        outputs = exodus
        # use initial damage time dependent
        build_param_use_initial_damage_time_dependent_mat = true
        build_param_peak_value = 0.7
        build_param_sigma = 5e2
        build_param_len_of_fault = 14000
        build_param_use_initial_damage_3D = true
        build_param_len_of_fault_dip = 10000
        build_param_center_point = '0 0 -10000'
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Debug
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain deviatroic_stress strain_invariant_ratio'
        outputs = exodus
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
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
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-10
    nl_max_its = 20
    nl_abs_tol = 1e-12
    # this is very robust, use as default
    petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
[]  

[Outputs]
    exodus = true    
[]

[BCs]
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []
    #Note: use neuamnnBC gives minimum waves than pressureBC  
    [static_pressure_left]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = left
        function = func_pos_xx_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = right
        function = func_neg_xx_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    #
    [static_pressure_front]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = front
        function = func_pos_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = back
        function = func_neg_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []
    #
    [static_pressure_front_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = front
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_left_shear]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = left
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right_shear]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = right
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []   
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = DirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]

[Functions]
    [func_pos_yy_stress]
        type = ParsedFunction      
        expression = 'if(-z<15600, -1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_yy_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'  
    []
    [func_pos_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600, -1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'
    []
    [func_pos_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, -1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, -1 * (-0.8 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
    [func_neg_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, 1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, 1 * (-0.8 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
[]