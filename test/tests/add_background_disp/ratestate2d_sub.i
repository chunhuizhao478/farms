#This file is only used to retrieve residual vector (elastic force) using "save_in" option at TIMESTEP_BEGIN of each time step

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 300
        ny = 300
        xmin = -15000
        xmax = 15000
        ymin = -15000
        ymax = 15000
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
    []
[]

[Variables]
    [disp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    [./resid_sub_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_sub_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y'
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          displacements = 'disp_sub_x disp_sub_y'
          save_in = 'resid_sub_x resid_sub_y'
        [../]
      [../]
    [../]
[]

[Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
[]

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
[]