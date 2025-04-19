[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole_narrownotch.msh'
  []
[]

[GlobalParams]

  displacements = 'disp_x disp_y'
    
  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 15.62e9
      
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 19.92e9
  
  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.9
  
  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.9
  
  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8
  
  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #if option 2, use Cd_constant
  Cd_constant = 80

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = 100

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  # CBCBH_multiplier = 0.0
  CBH_constant = 0

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = 0

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = 0.05

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = 0.05 #1e-3
  
  #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  C_g = 1e-12 #
  
  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  m1 = 10
  
  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  m2 = 1
  
  #coefficient of energy ratio Fb/Fs = chi < 1
  chi = 0.8
  
  #
  D = 0
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4 * t'
    []
    #strain
    [func_strain_xx]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = elastic_strain_00
    [../]
    [func_strain_xy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_01
    [../]
    [func_strain_xz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_02
    [../]
    [func_strain_yy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_11
    [../]
    [func_strain_yz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_12
    [../]
    [func_strain_zz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_22
    [../] 
  []

  [Physics/SolidMechanics/QuasiStatic]
    [./all]
      strain = small
      add_variables = true
      eigenstrain_names = static_initial_strain_tensor
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
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBM
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 xi stress'
        outputs = exodus
    []
    [./initial_strain]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_strain_tensor
      tensor_functions = 'func_strain_xx     func_strain_xy      func_strain_xz 
                          func_strain_xy     func_strain_yy      func_strain_yz
                          func_strain_xz     func_strain_yz      func_strain_zz'
    []
    [dummy_matprop]
      type = GenericConstantMaterial
      prop_names = 'initial_damage initial_breakage shear_stress_perturbation'
      prop_values = '0.0 0.0 0.0'  
    []  
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []
  
  [Executioner]
    type = Transient
    solve_type = Newton
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type'
    # petsc_options_value = 'lu       mumps NONZERO'
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # automatic_scaling = true
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 30
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    l_tol = 1e-5
    start_time = 0
    end_time = 4000
    dt = 0.01
    automatic_scaling = true
    [./TimeIntegrator]
        type = ImplicitEuler
    [../]
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []

  [ICs]
    [disp_x_ic]
        type = SolutionIC
        variable = disp_x
        solution_uo = init_sol_components
        from_variable = disp_x
    []
    [disp_y_ic]
        type = SolutionIC
        variable = disp_y
        solution_uo = init_sol_components
        from_variable = disp_y
    []
  []

  [UserObjects]
    [./init_sol_components]
        type = SolutionUserObject
        mesh = ./cdb_narrownotch_static_out.e
        system_variables = 'disp_x disp_y elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
        timestep = LATEST
        force_preaux = true
    [../]
  []