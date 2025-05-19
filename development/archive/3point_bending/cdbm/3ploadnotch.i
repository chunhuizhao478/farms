[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  'meshnotch.msh'
  []
  [./elastic_region_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    top_right = '0.15 0.1 0'
    block_id = 2
  []
  [./elastic_region_2]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_1
    bottom_left = '0.30 0 0'
    top_right = '0.45 0.1 0'
    block_id = 2
  []
  [./extranodeset_1]
    type = ExtraNodesetGenerator
    coord = '   0 0 0'
    new_boundary = support1
    input = elastic_region_2
  []
  [./extranodeset_2]
    type = ExtraNodesetGenerator
    coord =  '0.45 0 0'
    new_boundary = support2
    input = extranodeset_1
  []
  [./extranodeset_3]
    type = ExtraNodesetGenerator
    coord = '0.22 0.1 0'
    new_boundary = load
    input = extranodeset_2
  [] 
  displacements = 'disp_x disp_y'
[]

[GlobalParams]

  displacements = 'disp_x disp_y'
    
  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 5.56e9
      
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 8.33e9
  
  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.8
  
  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.8
  
  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8
  
  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #if option 2, use Cd_constant
  Cd_constant = 1e4

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
    [alpha_damagedvar_aux]
      order = FIRST
      family = MONOMIAL
    []
    [B_damagedvar_aux]
      order = FIRST
      family = MONOMIAL
    []
  []

  [AuxKernels]
    [get_damage]
      type = MaterialRealAux
      variable = alpha_damagedvar_aux
      property = alpha_damagedvar
      block = '0'
    []
    [get_Bdamage]
      type = MaterialRealAux
      variable = B_damagedvar_aux
      property = B_damagedvar
      block = '0'
    []
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4 * t'
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
  []
  
  [BCs]
    [fix_support_y]
      type = DirichletBC
      variable = disp_y
      boundary = support1
      value = 0
    []
    [fix_support2_x]
      type = DirichletBC
      variable = disp_x
      boundary = support2
      value = 0
    []
    [fix_support2_y]
      type = DirichletBC
      variable = disp_y
      boundary = support2
      value = 0
    []
    [apply_load]
      type = FunctionDirichletBC
      variable = disp_y
      boundary = load
      function = func_loading
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 20e9
      poissons_ratio = 0.2
    [../]
    [damage_mat]
      type = DamageBreakageMaterial
      output_properties = 'alpha_damagedvar B_damagedvar'
      outputs = exodus
      block = '1'
    []
    [stress_medium]
      type = ComputeLagrangianDamageBreakageStressPK2
      large_kinematics = true
      output_properties = 'pk2_stress strain_invariant_ratio'
      outputs = exodus
      block = '1'
    []
    [elastic_stress]
      type = ComputeStVenantKirchhoffStress
      block = 2
    []
    [dummy_matprop]
      type = GenericConstantMaterial
      prop_names = 'initial_damage initial_breakage shear_stress_perturbation'
      prop_values = '0.0 0.0 0.0'  
    []
    [strain]
      type = ComputeLagrangianStrain
      large_kinematics = true
      output_properties = 'deformation_gradient'
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
    type = Transient
    solve_type = Newton
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type'
    # petsc_options_value = 'lu       mumps NONZERO'
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # automatic_scaling = true
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 30
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0
    end_time = 4000
    dt = 0.1
    automatic_scaling = true
    [./TimeIntegrator]
        type = ImplicitEuler
    [../]
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []