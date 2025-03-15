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
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
    []
[]

[AuxKernels]
    [compute_xi]
        type = CompXi3D
        variable = xi_computed
        execute_on = 'TIMESTEP_END'
    []
[]

[Materials]
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 10e9
        shear_modulus = 10e9
    []
    [compute_stress]
        type = ComputeLinearElasticStress
        output_properties = 'stress'
        outputs = exodus
      []
    [compute_strain]
        type = ComputeSmallStrain
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
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = top
    []
    [bc_fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = bottom
    []
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    [bc_fix_left_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = left
    []
    [bc_fix_right_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = right
    []
    [./neumann_top_y]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -120e6
    [../]
    [./neumann_left_x]
        type = NeumannBC
        variable = disp_x
        boundary = left
        value = 135e6
    [../]
    [./neumann_right_x]
        type = NeumannBC
        variable = disp_x
        boundary = right
        value = -135e6
    [../]
[]