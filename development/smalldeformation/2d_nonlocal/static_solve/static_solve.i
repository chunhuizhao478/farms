#implicit continuum damage-breakage model dynamics

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
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '-600000 -600000 0'
        new_boundary = corner_ptr
        input = sidesets
    []
    [./extranodeset2]
        type = ExtraNodesetGenerator
        coord = '600000 -600000 0'
        new_boundary = corner_ptr
        input = extranodeset1
    []
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
        type = ADMaterialRealAux
        variable = xi_output
        property = xi_initial
    []
    [get_I2]
        type = ADMaterialRealAux
        variable = I2_output
        property = I2_initial
    []
    [get_alpha_damagedvar]
        type = ADMaterialRealAux
        variable = alpha_damagedvar_output
        property = initial_damage
    []
    [get_B_damagedvar]
        type = ADMaterialRealAux
        variable = B_damagedvar_output
        property = initial_breakage
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
    []
[]

[Materials]
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y'
        outputs = exodus
    [] 
    [stress]
        type = ADComputeDamageStressStaticDistributionDynamicCDBM
        lambda_o = 32.04e9
        shear_modulus_o = 32.04e9
        chi = 0.8
        xi_o = -0.8
        xi_d = -0.9
        outputs = exodus
        block = '1 3'
    [] 
    [elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
    []
    [stress_elastic]
        type = ADComputeLinearElasticStress
        block = '2'
    []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    [initial_damage_surround]
        type = ADInitialDamageCycleSim2D
        len_of_fault = 18000
        sigma = 5e2
        peak_val = 0.7
        output_properties = 'initial_damage'      
        outputs = exodus
    []
    [dummy_material]
        type = ADGenericConstantMaterial
        prop_names = 'initial_breakage'
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
[]

[BCs]
    #add initial shear stress
    [initial_shear_stress_top]
        type = ADNeumannBC
        variable = disp_x
        value = 14e6
        boundary = top
    [] 
    [initial_shear_stress_bottom]
        type = ADNeumannBC
        variable = disp_x
        value = -14e6
        boundary = bottom
    []
    # [initial_shear_stress_left]
    #     type = NeumannBC
    #     variable = disp_y
    #     value = -12e6
    #     boundary = top
    # [] 
    # [initial_shear_stress_right]
    #     type = NeumannBC
    #     variable = disp_y
    #     value = 12e6
    #     boundary = bottom
    # []
    # 
    [static_pressure_top]
        type = ADNeumannBC
        variable = disp_y
        boundary = top
        value = -50e6
        displacements = 'disp_x disp_y'
    []   
    [static_pressure_bottom]
        type = ADNeumannBC
        variable = disp_y
        boundary = bottom
        value = 50e6
        displacements = 'disp_x disp_y'
    []  
    [static_pressure_left]
        type = ADNeumannBC
        variable = disp_x
        boundary = left
        value = 50e6
        displacements = 'disp_x disp_y'
    []  
    [static_pressure_right]
        type = ADNeumannBC
        variable = disp_x
        boundary = right
        value = -50e6
        displacements = 'disp_x disp_y'
    []       
    # fix left ptr
    [./fix_cptr1_x]
        type = ADDirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr2_y]
        type = ADDirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []   
    # fix right ptr
    [./fix_cptr4_y]
        type = ADDirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
[]