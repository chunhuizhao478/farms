#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../../mesh/mesh_small.msh'
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
        type = ADMaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
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
        lambda_o = 30e9
        shear_modulus_o = 30e9
        xi_o = -0.8
        chi = 0.7
        xi_d = -0.9
        outputs = exodus
    []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    [initial_damage_surround]
        type = ADInitialDamageCycleSim2D
        len_of_fault = 8000
        sigma = 5e2
        peak_val = 0.7
        use_damage_perturb = true
        damage_perturb = 'damage_perturb'
        output_properties = 'initial_damage'      
        outputs = exodus
        block = 1
    []
    [define_damage_perturb]
        type = ADDamagePerturbationSquare2D
        nucl_center = '0 0'
        e_damage = 0.3
        length = 400
        thickness = 200
        duration = 1.0
        sigma = 1318.02
        block = 1
    []
    [const_damage_b2]
        type = ADGenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
        block = 2
    []
    [const_damage_b3]
        type = ADGenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
        block = 3
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
    nl_rel_tol = 1e-8
    nl_max_its = 20
    nl_abs_tol = 1e-10
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
    [bc_fix_bottom_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []  
    [bc_fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = bottom
    []
    # 
    [./Pressure]
        [static_pressure_top]
            boundary = top
            factor = 50e6
            displacements = 'disp_x disp_y'
        []    
        [static_pressure_left]
            boundary = left
            factor = 50e6
            displacements = 'disp_x disp_y'
        []  
        [static_pressure_right]
            boundary = right
            factor = 50e6
            displacements = 'disp_x disp_y'
        []     
    []        
    # fix ptr
    # [./fix_cptr1_x]
    #     type = ADDirichletBC
    #     variable = disp_x
    #     boundary = corner_ptr
    #     value = 0
    # []
    # [./fix_cptr2_y]
    #     type = ADDirichletBC
    #     variable = disp_y
    #     boundary = corner_ptr
    #     value = 0
    # []
    #
    #add initial shear stress
    [./initial_shear_stress]
        type = ADNeumannBC
        variable = disp_x
        value = 12e6
        boundary = top
    []    
[]