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
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '0 0 1'
        new_boundary = corner_ptr
        input = msh
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

[AuxVariables]
    [xi_computed]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [compute_xi]
        type = CompXi3D
        variable = xi_computed
        execute_on = 'TIMESTEP_END'
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
        type = ComputeStVenantKirchhoffStress
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
    [fix_bottom_z]
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
        type = DirichletBC
        variable = disp_x
        boundary = top
        value = 0
    []
    [./Pressure]
        [static_pressure_back]
            boundary = back
            factor = 80e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []  
        [static_pressure_front]
            boundary = front
            factor = 80e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []    
        [static_pressure_left]
            boundary = left
            factor = 80e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []  
        [static_pressure_right]
            boundary = right
            factor = 80e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []         
    []    
[]