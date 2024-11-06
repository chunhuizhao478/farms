[Mesh]
    [msh]
        type = GeneratedMeshGenerator
        dim = 1
        nx = 1000
        xmin = 0
        xmax = 10
    []
[]

[GlobalParams]
    displacements = 'disp_x'
    use_displaced_mesh = true
[]


[Variables]
    [disp_x]
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
[]

[Kernels]
    # [dispkernel_x]
    #     type = StressDivergenceTensors
    #     variable = disp_x
    #     component = 0
    # []
    [dispkernel_x]
        type = TotalLagrangianStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
[]

[Materials]
    # small strain
    # [strain]
    #     type = ComputeSmallStrain
    # []
    # [compute_stress]
    #     type = ComputeLinearElasticStress
    # []
    # [compute_elasticity]
    #     type = ComputeIsotropicElasticityTensor
    #     youngs_modulus = 25e9
    #     poissons_ratio = 0.2
    # []
    # finite strain
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
    []
    [compute_elasticity]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 25e9
        poissons_ratio = 0.2
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
    nl_abs_tol = 1e-12
    nl_rel_tol = 1e-8
[]

[Outputs]
    [./exodus]
      type = Exodus
    [../]
[]

[BCs]
    [fix_left_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = left
    []
    [applied_force_right_x]
        type = NeumannBC
        variable = disp_x
        value = 1e3
        boundary = right
    []
[]