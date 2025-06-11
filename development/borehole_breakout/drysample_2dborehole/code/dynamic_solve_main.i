#implicit continuum damage-breakage model dynamics
#solid properties
#-------------------------------------------------#
solid_density = 2640
youngs_modulus = 48.5e9
poissons_ratio = 0.22
lambda_o = ${fparse youngs_modulus*poissons_ratio/(1+poissons_ratio)/(1-2*poissons_ratio)}
shear_modulus_o = ${fparse youngs_modulus/(2*(1+poissons_ratio))}
length_scale = 1.3e-3
#-------------------------------------------------#

#damage-breakage properties
#-------------------------------------------------#
xi_o = -0.8073
xi_d = -0.8073
Cg = 1e-12
Cd_constant = 80
CdCb_multiplier = 100
beta_width = 0.05
CBH_constant = 0
C_1 = 0
C_2 = 0.05
chi = 0.8
m1 = 10
m2 = 1
#-------------------------------------------------#

#strain rate dependent Cd options
#-------------------------------------------------#
use_cd_strain_dependent = true
m_exponent = 0.8
strain_rate_hat = 1e-4
cd_hat = 1
#-------------------------------------------------#

#boundary conditions
#-------------------------------------------------#
outer_confinement_pressure = 17.2e6
loading = 'if (t > 1e-3, -2.6477e-5 - 3.3e-7 * t, -2.6477e-5)'
#-------------------------------------------------#

#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/2dborehole_order1.msh'
    [] 
[]

[GlobalParams]

    displacements = 'disp_x disp_y'
      
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = ${lambda_o}
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = ${shear_modulus_o}
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = ${xi_o}
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = ${xi_d}
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = ${Cg}
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = ${m1}
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = ${m2}
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = ${chi}

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
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    #
    [alpha_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    [B_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
    [I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [deviatroic_strain_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [structural_stress_coefficient_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [gradx_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    [grady_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    #spatial damage parameters
    [cg_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
    [nonlocal_xi]
        order = FIRST
        family = MONOMIAL
    []
    
[]

[AuxKernels]
    [vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    #
    [get_xi]
        type = MaterialRealAux
        variable = xi_aux
        property = strain_invariant_ratio
        block = '8'
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_aux
        property = second_elastic_strain_invariant
        block = '8'
    [] 
    [get_deviatroic_strain_rate]
        type = MaterialRealAux
        variable = deviatroic_strain_rate_aux
        property = deviatroic_strain_rate
        block = '8'
    []
    #
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
    []
[]

[Kernels]
    [inertia_x]
        type = InertialForce
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
    [inertia_y]
        type = InertialForce
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
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

[Materials]
    [strain]
        type = ComputeSmallStrain
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBMDiffused
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        output_properties = 'stress elastic_strain_tensor plastic_strain_tensor total_strain_tensor strain_invariant_ratio'
        outputs = exodus
        block = '8'
    []
    #elastic material
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = ${youngs_modulus}
        poissons_ratio = ${poissons_ratio}
    []
    [compute_stress]
        type = ComputeLinearElasticStress
        outputs = exodus
        block = '6 7'
    []
    #strain invariant ratio
    [comp_strain_invariant_ratio]
        type = ComputeXi 
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
        block = '6 7'
    []
    #nonlocal eqstrain
    [nonlocal_eqstrain]
        type = ElkNonlocalEqstrain
        average_UO = eqstrain_averaging
        output_properties = 'eqstrain_nonlocal'
        outputs = exodus
    []
    #shear stress perturbation
    [damage_perturbation]
        type = GenericConstantMaterial
        prop_names = 'shear_stress_perturbation damage_perturbation'
        prop_values = '0.0 0.0' #initial perturbation values
    []
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = ${solid_density}
    []  
[]

[UserObjects]
    [eqstrain_averaging] #length scale = radius = grain size 
        type = ElkRadialAverage
        length_scale = ${length_scale}
        prop_name = strain_invariant_ratio
        radius = ${length_scale}
        weights = BAZANT
        execute_on = TIMESTEP_END
    []
[]

#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = ${loading}
    []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */inertia_x */inertia_y'
    start_time = -1e-12
    end_time = 1e-3 # dt used in the simulation
  []
[../]
  
[Executioner]
    type = Transient
    # solve_type = 'NEWTON'
    solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 1e10
    # num_steps = 10
    l_max_its = 30
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 10
    nl_abs_tol = 1e-12
    # petsc_options_iname = '-ksp_type -ksp_max_it -ksp_gmres_restart -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    # petsc_options_value = 'gmres          1000      100       hypre  boomeramg True'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
    petsc_options_value = ' lu       mumps       100'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    # automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    # dt = 10
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-3
        cutback_factor_at_failure = 0.5
        optimal_iterations = 20
        growth_factor = 1.25
        max_time_step_bound = 10
    []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 10 ###
        # show = 'vel_x vel_y vel_z alpha_damagedvar_aux B_damagedvar_aux xi_aux deviatroic_strain_rate_aux nonlocal_xi stress_22 elastic_strain_tensor_22 plastic_strain_tensor_22 total_strain_tensor_22'
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
    [../]
[]

[BCs]
    #fix bottom boundary
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 3
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 3
        value = 0
    []
    #applied load on top boundary
    [applied_top_z_dispload]
        type = FunctionDirichletBC
        variable = disp_y
        boundary = 2
        function = applied_load_top
    [] 
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = '4 5'
          factor = ${outer_confinement_pressure}
          displacements = 'disp_x disp_y'
        [../]
    []
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'dynamic_solve_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
        cli_args = 'Cd_constant=${Cd_constant};CdCb_multiplier=${CdCb_multiplier};beta_width=${beta_width};lambda_o=${lambda_o};shear_modulus_o=${shear_modulus_o};xi_o=${xi_o};xi_d=${xi_d};CBH_constant=${CBH_constant};C_1=${C_1};C_2=${C_2};use_cd_strain_dependent=${use_cd_strain_dependent};m_exponent=${m_exponent};strain_rate_hat=${strain_rate_hat};cd_hat=${cd_hat}'
        clone_parent_mesh = true
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_damagedvar_sub B_damagedvar_sub structural_stress_coefficient_sub'
        variable = 'alpha_damagedvar_aux B_damagedvar_aux structural_stress_coefficient_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'I2_aux nonlocal_xi deviatroic_strain_rate_aux' 
        variable = 'I2_sub_aux xi_sub_aux deviatroic_strain_rate_sub_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]