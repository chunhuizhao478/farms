[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole_narrownotch.msh'
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
  []

  [AuxVariables]
    [./strength]
      order = CONSTANT
      family = MONOMIAL
      initial_condition = 6.43e6
    [../]
    [crack_damage_aux]
      order = FIRST
      family = MONOMIAL
    []
    [./xi]
      order = CONSTANT
      family = MONOMIAL         
    []
  []

  [AuxKernels]
    [get_xi]
        type = CompXi3D
        variable = xi
    []
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4'
    []
  []

  [GlobalParams]
    displacements = 'disp_x disp_y'
  [] 

  [Physics/SolidMechanics/QuasiStatic]
    [./all]
      strain = small
      add_variables = true
    [../]
  []
  
  [BCs]
    [fix_support_x]
      type = DirichletBC
      variable = disp_x
      boundary = right
      value = 0
    []
    [fix_support_y]
      type = DirichletBC
      variable = disp_y
      boundary = right
      value = 0
    []
    [apply_load]
      type = FunctionDirichletBC
      variable = disp_x
      boundary = left
      function = func_loading
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 48.5e9
      poissons_ratio = 0.22
    [../]
    [./elastic_stress]
      type = ComputeLinearElasticStress
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
  solve_type = Newton
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  petsc_options_value = '101                asm      lu'
  automatic_scaling = true
  line_search = 'none'
  # num_steps = 1
  l_max_its = 100
  nl_max_its = 10
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_tol = 1e-5
[]
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []