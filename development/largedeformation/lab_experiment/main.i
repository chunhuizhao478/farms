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
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 10e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 10e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
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
    Cd_constant = 300

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 3e5

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 3e6

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.001 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-5
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.5
    
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
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
    []
    # damage
    [damage_mat]
        type = DamageBreakageMaterial
    [] 
    [initial_damage]
        type = GenericConstantMaterial
        prop_names = initial_damage
        prop_values = 0
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2
        large_kinematics = true
        output_properties = 'pk1_stress'
        outputs = exodus
    []
    # elastic
    # [elastic_tensor]
    #     type = ComputeIsotropicElasticityTensor
    #     lambda = 1e10
    #     shear_modulus = 1e10
    # []
    # [compute_stress]
    #     type = ComputeLagrangianLinearElasticStress
    #     large_kinematics = true
    #     outputs = exodus
    # []
[]  

[Functions]
[]
  
[Executioner]
    type = Transient
    solve_type = Newton
    dt = 1
    num_steps = 1
    abort_on_solve_fail = true
    nl_abs_tol = 1e-6
    nl_rel_tol = 1e-8
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

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = 'initial_load_check_out.e'
      system_variables = 'disp_x disp_y disp_z'
      timestep = LATEST
      force_preaux = true
    [../]
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
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]