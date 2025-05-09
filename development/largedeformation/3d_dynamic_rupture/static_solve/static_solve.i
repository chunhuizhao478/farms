#implicit continuum damage-breakage model dynamics

#material properties
# lambda_o = 32.04e9
# shear_modulus_o = 32.04e9
# xi_o = -0.8
# xi_d = -0.9
# chi = 0.8
fluid_density = 1000   
solid_density = 2700
gravity_pos = 9.81
gravity_neg = -9.81

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/mesh_large.msh'
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
        coord = ' -120000 -120000 -120000;
                   120000 -120000 -120000;
                   120000 120000  -120000;
                  -120000 120000  -120000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9
        
    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9
    
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
    Cd_constant = 1e5

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 1e-4

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
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
    
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
    [xi_output]
        order = FIRST
        family = MONOMIAL
    []
    [I2_output]
        order = FIRST
        family = MONOMIAL
    []
    [alpha_damagedvar_output]
        order = FIRST
        family = MONOMIAL
    []
    [B_damagedvar_output]
        order = FIRST
        family = MONOMIAL
    []
[]

[AuxKernels]
    [get_xi]
        type = MaterialRealAux
        variable = xi_output
        property = strain_invariant_ratio
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_output
        property = second_elastic_strain_invariant
    []
    [get_alpha_damagedvar]
        type = MaterialRealAux
        variable = alpha_damagedvar_output
        property = alpha_damagedvar
    []
    [get_B_damagedvar]
        type = MaterialRealAux
        variable = B_damagedvar_output
        property = B_damagedvar
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
        value = ${gravity_neg}
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
        # output_properties = 'alpha_damagedvar B_damagedvar shear_modulus_o_mat shear_modulus'
        # outputs = exodus
        # use initial damage time dependent
        build_param_use_initial_damage_time_dependent_mat = true
        build_param_peak_value = 0.7
        build_param_sigma = 5e2
        build_param_len_of_fault = 28000
        build_param_use_initial_damage_3D = true
        build_param_len_of_fault_dip = 15000
        build_param_center_point = '0 0 -7500'
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Debug
        large_kinematics = true
        # output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain deviatroic_stress strain_invariant_ratio second_elastic_strain_invariant'
        # outputs = exodus
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
    []
    #shear stress perturbation
    [damage_perturbation]
        type = PerturbationRadial
        nucl_center = '0 0 0'
        peak_value = 0
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'shear_stress'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
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
    nl_rel_tol = 1e-10
    nl_max_its = 20
    nl_abs_tol = 1e-12
    # petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    # petsc_options_value = 'gmres     hypre  True'
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
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
    show = 'disp_x disp_y disp_z xi_output I2_output alpha_damagedvar_output B_damagedvar_output'
[]

#parameters for the initial stress field
################################################
bxx = 0.926793
byy = 1.073206
bxy = -0.169029
linear_variation_cutoff_distance = 15600
################################################
[Functions]
    [func_pos_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
[]

[BCs]
    [fix_bottom_z]
        type = ADDirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []
    #Note: use neuamnnBC gives minimum waves than pressureBC  
    [static_pressure_left]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = left
        function = func_pos_xx_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = right
        function = func_neg_xx_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    #
    [static_pressure_front]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = front
        function = func_pos_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = back
        function = func_neg_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []
    #
    [static_pressure_front_shear]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = front
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_left_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = left
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = right
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []   
    # fix ptr
    [./fix_cptr1_x]
        type = ADDirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = ADDirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = ADDirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]