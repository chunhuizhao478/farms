#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 400
        ny = 80
        xmin = 0
        xmax = 0.05
        ymin = 0
        ymax = 0.01
    []
    [./box]
        type = SubdomainBoundingBoxGenerator
        input = msh
        block_id = 1
        bottom_left = '0 0 0'
        top_right = '0.05 0.004 0'
    []
    [./box2]
        type = SubdomainBoundingBoxGenerator
        input = box
        block_id = 0
        bottom_left = '0 0.004 0'
        top_right = '0.05 0.006 0'
    []
    [./box3]
        type = SubdomainBoundingBoxGenerator
        input = box2
        block_id = 2
        bottom_left = '0 0.006 0'
        top_right = '0.05 0.01 0'
    []
[]

[GlobalParams]

    ##----continuum damage breakage model----##
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9

    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9

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

    #if option 2, use Cd_constant #specify by auxiliary variable
    Cd_constant = 1e6

    #strain rate dependent Cd options
    m_exponent = 0.8
    strain_rate_hat = 1e-8
    cd_hat = 10

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd #specify by auxiliary variable
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    CBH_constant = 1e6

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf" #specify by auxiliary variable
    C_1 = 30

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3

    #diffusion parameter #close the gradient
    D_diffusion = 0

[]

#this sub-app solves damage/breakage evolution equations
[Variables]
    [alpha_damagedvar_sub]
        order = FIRST
        family = LAGRANGE
    []
    [B_damagedvar_sub]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    [xi_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    [I2_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_damage_sub_aux]
        order = FIRST
        family = LAGRANGE
    []
    #deviatroic_strain_rate
    [deviatroic_strain_rate_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [bounds_dummy]
        order = FIRST
        family = LAGRANGE
    []
    #
    [structural_stress_coefficient_sub]
        order = FIRST
        family = MONOMIAL
    []
    #
    [Cd_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xio_aux]
        order = FIRST
        family = LAGRANGE
    []
    [xid_aux]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    #damagevar
    [time_derivative_alpha]
        type = TimeDerivative
        variable = alpha_damagedvar_sub
    []
    [diffusion_alpha]
        type = DamageEvolutionDiffusion
        variable = alpha_damagedvar_sub
        coupled = B_damagedvar_sub
        block = '0'
    []
    [forcing_term_alpha]
        type = DamageEvolutionConditionalForcing
        variable = alpha_damagedvar_sub
        coupled = B_damagedvar_sub
        block = '0'
    []
    # [perturb_source_alpha]
    #     type = PerturbationSource
    #     variable = alpha_damagedvar_sub
    #     damage_source = 'damage_perturbation'
    #     block = '1 3'
    # []
    #breakagevar
    [time_derivative_B]
        type = TimeDerivative
        variable = B_damagedvar_sub
    []
    [forcing_term_B]
        type = BreakageEvolutionConditionalForcing
        variable = B_damagedvar_sub
        coupled = alpha_damagedvar_sub
        block = '0'
    []
    # [perturb_source_b]
    #     type = PerturbationSource
    #     variable = B_damagedvar_sub
    #     damage_source = 'damage_perturbation'
    #     block = '1 3'
    # []
[]

[Bounds]
    [alpha_damagedvar_upper_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = alpha_damagedvar_sub
      bound_type = upper
      bound_value = 1
    []
    [alpha_damagedvar_lower_bound]
      type = VariableConstantBounds
      variable = bounds_dummy
      bounded_variable = alpha_damagedvar_sub
      bound_type = lower
      bound_value = initial_damage_sub_aux
    []
    [B_damagedvar_upper_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = B_damagedvar_sub
      bound_type = upper
      bound_value = 1
    []
    [B_damagedvar_lower_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = B_damagedvar_sub
      bound_type = lower
      bound_value = 0
    []
[]

[AuxKernels]
    [get_Cd]
        type = MaterialRealAux
        variable = Cd_aux
        property = Cd
    []
[]

[Materials]
    # # damage
    [damage_mat]
        type = DiffusedDamageBreakageMaterialSubApp
        I2_aux = I2_sub_aux
        xi_aux = xi_sub_aux
        initial_damage_aux = initial_damage_sub_aux
        #use strain rate dependent Cd
        use_cd_strain_dependent = false
        strain_rate = deviatroic_strain_rate_sub_aux
        zero_cd_below_hat = true
    []
    # #add shear perturbation to the system
    # [damage_perturbation]
    #     type = PerturbationRadialSource
    #     nucl_center = '0 0 0'
    #     peak_value = 1.0
    #     thickness = 200
    #     length = 2000
    #     duration = 10
    #     perturbation_type = 'damage'
    #     sigma_divisor = 2.0
    #     # output_properties = 'shear_stress_perturbation damage_perturbation'
    #     # outputs = exodus
    # []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    start_time = -1e-12
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 30
    nl_abs_tol = 1e-10
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
    petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'
    verbose = true
[]
