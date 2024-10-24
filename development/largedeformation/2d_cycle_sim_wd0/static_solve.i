#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = './meshfile/tpv2052dm.msh'
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
    [xi_computed]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        variable = disp_y
        component = 1
    []
[]

[Materials]
    [stress_nucleation]
        type = ADComputeDamageStressStaticDistribution
        lambda_o = 10e9
        shear_modulus_o = 10e9
        xi_o = -0.8
        gamma_damaged_r = 1.1595e10
        outputs = exodus
    []
    [compute_xi]
        type = ADComputeXi
    []
    [initial_damage_strip]
        type = ADGenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.7'
        block = 1
        output_properties = 'initial_damage'
        outputs = exodus
    []
    [initial_damage_surround]
        type = ADInitialDamageCycleSim2D
        output_properties = 'initial_damage'
        outputs = exodus
        block = 3
    []
    [initial_damage_zero]
        type = ADGenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0'
        block = '2'
        output_properties = 'initial_damage'
        outputs = exodus
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y'
        outputs = exodus
    []     
[] 
  
[Executioner]
    type = Steady
    solve_type = 'NEWTON'
    automatic_scaling = true
[]  

[Outputs]
    exodus = true   
[]

[BCs]
    [bc_load_top_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = top
    []
    [bc_fix_bottom_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = bottom
    []
    [bc_fix_bottom_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    [bc_fix_left_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = left
    []
    [bc_fix_right_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = right
    []
    [./neumann_top_y]
        type = ADNeumannBC
        variable = disp_y
        boundary = top
        value = -120e6
    [../]
    [./neumann_left_x]
        type = ADNeumannBC
        variable = disp_x
        boundary = left
        value = 135e6
    [../]
    [./neumann_right_x]
        type = ADNeumannBC
        variable = disp_x
        boundary = right
        value = -135e6
    [../]
[]