#implicit continuum damage-breakage model dynamics
#solid properties
#-------------------------------------------------#
# solid_density = 2640
youngs_modulus = 48.5e9
poissons_ratio = 0.22
# solid_bulk_compliance = 3.46e-11
lambda_o = ${fparse youngs_modulus*poissons_ratio/(1+poissons_ratio)/(1-2*poissons_ratio)}
shear_modulus_o = ${fparse youngs_modulus/(2*(1+poissons_ratio))}
#-------------------------------------------------#

#damage-breakage properties
#-------------------------------------------------#
xi_o = -0.8073
xi_d = -0.8073
# Cg = 1e-12
Cd_constant = 80
CdCb_multiplier = 100
beta_width = 0.05
#-------------------------------------------------#

[Mesh]
    [./msh]
        type = FileMeshGenerator
        # file = '../meshfile/mesh_adaptive_test.msh'
        file = '../mesh/2dborehole_order1.msh'
    [] 
[]

[GlobalParams]
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = ${lambda_o}
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = ${shear_modulus_o}
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = ${xi_o}
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = ${xi_d}
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant #specify by auxiliary variable
    Cd_constant = ${Cd_constant}

    #strain rate dependent Cd options
    # m_exponent = 0.8
    # strain_rate_hat = 1e-8
    # cd_hat = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd #specify by auxiliary variable
    CdCb_multiplier = ${CdCb_multiplier}

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf" #specify by auxiliary variable
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = ${beta_width}

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
        block = '8'
    []
    [forcing_term_alpha]
        type = DamageEvolutionConditionalForcing
        variable = alpha_damagedvar_sub
        coupled = B_damagedvar_sub
        block = '8'
    []
    # [perturb_source_alpha]
    #     type = PerturbationSource
    #     variable = alpha_damagedvar_sub
    #     damage_source = 'damage_perturbation'
    #     block = '3'
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
        block = '8'
    []
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

[Materials]
    # # damage
    [damage_mat]
        type = DiffusedDamageBreakageMaterialSubApp
        I2_aux = I2_sub_aux
        xi_aux = xi_sub_aux
        initial_damage_aux = initial_damage_sub_aux
        strain_rate = deviatroic_strain_rate_sub_aux
    []
    #add shear perturbation to the system
    [damage_perturbation]
        type = PerturbationRadialSource
        nucl_center = '0 0 0'
        peak_value = 0.0
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'damage'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
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
    solve_type = 'NEWTON'
    start_time = -1e-12
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 10
    nl_abs_tol = 1e-10
    petsc_options_iname = '-snes_type'
    petsc_options_value = 'vinewtonrsls'
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-3
        cutback_factor_at_failure = 0.5
        optimal_iterations = 10
        growth_factor = 1.25
        max_time_step_bound = 100
    []
    [./TimeIntegrator]
        type = ImplicitEuler
        # type = BDF2
        # type = CrankNicolson
    [../]
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 50
        # show = 'Cd_aux'
    [../]
[]