#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_adaptive.msh'
    [] 

[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    porepressure = 'porepressure'
      
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 15.62e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 19.92e9

    # Water bulk modulus (2.2 GPa)
    fluid_bulk_modulus = 2.2e9      

     # Initial permeability (1 milli-darcy) 
    permeability_solid_o = 1e-20 

     # Initial porosity (15%)
    porosity_solid_o = 0.008   
     
    # Solid grains bulk modulus (36 GPa - typical for quartz)  
    solid_bulk_modulus_s = 50.00e9 

    initial_viscosity_fluid = 1e-3
    

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
    [porepressure]
        order = FIRST
        family = LAGRANGE
    []
    
[]

[AuxVariables]
    [initial_I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_xi_aux]
        order = FIRST
        family = MONOMIAL
   []
[]

[AuxKernels]
    [get_initial_I2]
        type = MaterialRealAux
        variable = initial_I2_aux
        property = I2_initial
    []
    [get_initial_xi]
        type = MaterialRealAux
        variable = initial_xi_aux
        property = xi_initial
    []
[]

[Kernels]
    [grad_stress_x]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
    [grad_stress_y]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_y
        component = 1
        large_kinematics = true
    []
    [grad_stress_z]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_z
        component = 2
        large_kinematics = true
    []
    [./mass1]
        type = FluidSolidCoupling
        variable = porepressure
    [../]
    [./mass2]
        type = PorePressureTimeDerivative
        variable = porepressure
    [../]
    [./darcy_flow]
        type = FluidDiffusion
        variable = porepressure
        large_kinematics = false
    []
[]

[Materials]
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        # outputs = exodus
    []
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 48.5e9
        poissons_ratio = 0.22
    []
    [compute_stress]
        type = ComputePoroStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus

    []
    [comp_strain_invariant_ratio]
        type = ComputeXi 
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
    []
    [porous_prop]
        type =   IntactPorousSolidProperties
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
   # solve_type = 'PJFNK'
    solve_type = 'NEWTON'
    start_time = -1e-12
    end_time = 3036
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 20
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
  #  automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'bt'
    # dt = 10
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1
        cutback_factor_at_failure = 0.5
        optimal_iterations = 6
        growth_factor = 1.25
        max_time_step_bound = 6
    []
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 1 ###
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
    [../]
    [checkpoint]
        type = Checkpoint
        time_step_interval = 20
        num_files = 2
    []
[]

[Functions]
  # 1. Axial Stress Function
  [axial_stress_function]
    type = ParsedFunction
    expression = 'if(t <= 1236, 1e6 + 16666.667*t, 21.6e6)'
  []

  # 2. Confining Stress Function
  [confining_stress_function]
    type = ParsedFunction
    expression = 'if(t <= 1236, 16666.667*t, 20.6e6)'
  []

  # 3. Pore Pressure Function
  [pore_pressure_function]
    type = ParsedFunction
    expression = 'if(t < 408, 0, if(t <= 1236, 4106.28*t, 3.4e6))'
  []
[]

[BCs]
    #fix bottom boundary
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 7
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 7
        value = 0
    []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = 7
        value = 0
    []
    #applied load on top boundary
    [applied_top_z_dispload]
        type = FunctionNeumannBC
        variable = disp_z
        boundary = 6
        function = axial_stress_function
    [] 
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = 4
          function = confining_stress_function
          displacements = 'disp_x disp_y'
        [../]
    []
    [./PorePressure]
          type = FunctionDirichletBC
          boundary = 5
          variable = porepressure
          function = pore_pressure_function
    []
[]