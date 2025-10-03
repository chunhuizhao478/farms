confining_pressure = 100e6

#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = './meshfile/mesh_wohole_coarse.msh'
    [] 
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9
        
    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -1.0
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -1.0
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-12
    
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
    [disp_z]
        order = FIRST
        family = LAGRANGE    
    []
[]

[AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
    #
    [vel_finitediff_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_finitediff_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_finitediff_z]
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
    #outputs
    [pk2_stress_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [green_lagrange_elastic_strain_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [plastic_strain_01]
        order = CONSTANT
        family = MONOMIAL       
    []
    [total_lagrange_strain_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [fy]
    []
[]

[AuxKernels]
    # [accel_x]
    #     type = NewmarkAccelAux
    #     variable = accel_x
    #     displacement = disp_x
    #     velocity = vel_x
    #     beta = 0.25
    #     execute_on = 'TIMESTEP_END'
    # []
    [vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    # [accel_y]
    #     type = NewmarkAccelAux
    #     variable = accel_y
    #     displacement = disp_y
    #     velocity = vel_y
    #     beta = 0.25
    #     execute_on = 'TIMESTEP_END'
    # []
    [vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_END'
    []
    #
    [get_xi]
        type = MaterialRealAux
        variable = xi_aux
        property = strain_invariant_ratio
        block = '3'
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_aux
        property = second_elastic_strain_invariant
        block = '3'
    [] 
    [get_deviatroic_strain_rate]
        type = MaterialRealAux
        variable = deviatroic_strain_rate_aux
        property = deviatroic_strain_rate
        block = '3'
    []
    #
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
    []
    #get outputs
    [get_pk2_stress_01]
        type = MaterialRankTwoTensorAux
        variable = pk2_stress_01
        property = pk2_stress
        i = 0
        j = 1
        block = '3'
    []
    [get_green_lagrange_elastic_strain_01]
        type = MaterialRankTwoTensorAux
        variable = green_lagrange_elastic_strain_01
        property = green_lagrange_elastic_strain
        i = 0
        j = 1
        block = '3'
    []
    [get_plastic_strain_01]
        type = MaterialRankTwoTensorAux
        variable = plastic_strain_01
        property = plastic_strain
        i = 0
        j = 1
        block = '3'   
    []
    [get_total_lagrange_strain_01]
        type = MaterialRankTwoTensorAux
        variable = total_lagrange_strain_01
        property = total_lagrange_strain
        i = 0
        j = 1
        block = '3'
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

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-3.3e-7 * t - (0.105/8.01e10) * ${confining_pressure}'
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        # outputs = exodus
    []
    # # damage
    [damage_mat]
        type = DiffusedDamageBreakageMaterialMainApp
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        structural_stress_coefficient = structural_stress_coefficient_aux
        #build L matrix using velocity
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
        #use cg
        # use_spatial_cg = true
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Diffused
        large_kinematics = true
        # output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain strain_invariant_ratio'
        # outputs = exodus
        block = '3'
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
    []
    [define_shear_stress_perturbation]
        type = PerturbationRadialAdvanced
        nucl_center = '0 0 0'
        peak_value = 0
        thickness = 1000
        length = 2000
        duration = 1
        perturbation_type = mean_stress
        temporal_ramp = cosine
        thickness_taper = cosine
        hold_after_duration = true
        sigma_divisor = 2.0
        output_properties = 'mean_stress_perturbation'
        outputs = exodus
    []
    #elastic material
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        # output_properties = 'green_lagrange_strain pk2_stress'
        # outputs = exodus
        block = '1 2'
    []
    #strain invariant ratio
    [comp_strain_invariant_ratio]
        type = ComputeXi 
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
        block = '1 2'
    []
    #nonlocal eqstrain
    [nonlocal_eqstrain]
        type = ElkNonlocalEqstrain
        average_UO = eqstrain_averaging
        output_properties = 'eqstrain_nonlocal'
        outputs = exodus
    []
[] 

[UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 1e-2
        prop_name = strain_invariant_ratio
        radius = 1e-2
        weights = BAZANT3D
        execute_on = TIMESTEP_END
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
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 1e5
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 30
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
    # petsc_options_iname = '-ksp_type -pc_type'
    # petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'basic'
    dt = 10
    # verbose = true
    fixed_point_max_its = 10
    accept_on_max_fixed_point_iteration = false
    fixed_point_rel_tol = 1e-6
    fixed_point_abs_tol = 1e-8
    [./TimeIntegrator]
        type = ImplicitEuler
    []
[]

[Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 1
      show = 'vel_x vel_y alpha_damagedvar_aux B_damagedvar_aux xi_aux nonlocal_xi pk2_stress_01 green_lagrange_elastic_strain_01 plastic_strain_01 total_lagrange_strain_01 deviatroic_strain_rate_aux mean_stress_perturbation'
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
        execute_on = 'TIMESTEP_END'
    [../]
    [out]
        type = Checkpoint
        time_step_interval = 200
        num_files = 2
    []
[]

[BCs]
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 6
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 6
        value = 0
    []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = 6
        value = 0
    []
    [applied_top_z]
        type = FunctionDirichletBC
        variable = disp_z
        boundary = 5
        function = applied_load_top
    []
    [./Pressure]
        #assign pressure on inner surface
        [pressure_inner]
            boundary = 4 #confirm the boundary id
            factor = ${confining_pressure}
        []             
    []  
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'dynamic_solve_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
        sub_cycling = true
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

#compute the reaction force on the top boundary
[Postprocessors]
  [Fy]
    type = SidesetReaction
    boundary = 5
    direction = '0 0 -1'
    stress_tensor = pk1_stress        # matches the Total Lagrangian formulation
  []
[]
