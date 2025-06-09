#continuum damage-breakage model dynamics

#material properties
lambda_o = 15.62e9
shear_modulus_o = 19.92e9
xi_o = -0.8073
xi_d = -0.8073
chi = 0.8

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_adaptive.msh'
    []  
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
    [initial_breakage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_cd_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    #stress function checks
    [stress_xx_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_yy_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_xy_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_zz_initial]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [get_initial_damage]
        type = ADMaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
    []
    [get_initial_breakage]
        type = ADMaterialRealAux
        variable = initial_breakage_aux
        property = initial_breakage
    []
    [get_initial_I2]
        type = ADMaterialRealAux
        variable = initial_I2_aux
        property = I2_initial
    []
    [get_initial_xi]
        type = ADMaterialRealAux
        variable = initial_xi_aux
        property = xi_initial
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
[]

[Materials]
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        outputs = exodus
    [] 
    [stress]
        type = ADComputeDamageStressStaticDistributionDynamicCDBM
        lambda_o = ${lambda_o}
        shear_modulus_o = ${shear_modulus_o}
        xi_o = ${xi_o}
        chi = ${chi}
        xi_d = ${xi_d}
        outputs = exodus
        block = '3'
    [] 
    [stress_elastic]
        type = ADComputeLinearElasticStress
        output_properties = 'stress'
        outputs = exodus
        block = '1 2'
    []
    [elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = ${lambda_o}
        shear_modulus = ${shear_modulus_o}
    []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    [initial_damage_surround]
        type = ADGenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
    []
    [initial_breakage_surround]
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
        type = DirichletBC
        variable = disp_z
        boundary = 6
        value = -2.6477e-5
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