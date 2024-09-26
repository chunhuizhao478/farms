[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 3
        nx = 10
        ny = 10
        nz = 10
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
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
    []
    [compute_stress]
        type = ComputeLagrangianLinearElasticStress
        large_kinematics = true
        outputs = exodus
    []
    [compute_strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
    []
[]  

[Executioner]
    type = Steady
[]

[Outputs] 
    exodus = true
[]

[BCs]
    [applied_disp_front]
        type = DirichletBC
        variable = disp_z
        value = -1e-4
        boundary = front
    []
    [applied_disp_back]
        type = DirichletBC
        variable = disp_z
        value = 1e-4
        boundary = back
    []
    [applied_disp_left]
        type = DirichletBC
        variable = disp_x
        value = 1e-4
        boundary = left
    []
    [applied_disp_right]
        type = DirichletBC
        variable = disp_x
        value = -1e-4
        boundary = right
    []
    [applied_disp_top]
        type = DirichletBC
        variable = disp_y
        value = -1e-4
        boundary = top
    []
    [applied_disp_bottom]
        type = DirichletBC
        variable = disp_y
        value = 1e-4
        boundary = bottom
    []
[]