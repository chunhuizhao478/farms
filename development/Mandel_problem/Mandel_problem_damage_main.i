[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10 
  ny = 1
  nz = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
  zmin = 0
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  porepressure = 'porepressure'

      ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 0.5
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 0.75
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8073
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.8073
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8

    # Water bulk modulus (2.2 GPa)
    fluid_bulk_modulus = 8       

     # Initial permeability (1 milli-darcy) 
    permeability_solid_o = 1.5 

     # Initial porosity (15%)
    porosity_solid_o = 0.5   
     
     # Granular bulk modulus (15 GPa)      
    solid_bulk_modulus_g = 1 

    # Solid grains bulk modulus (36 GPa - typical for quartz)  
    solid_bulk_modulus_s = 2.5   

    permeability_evolution_with_damage = 0
    initial_grain_size = 1 
    ultimate_grain_size = 1
    initial_viscosity_fluid = 1

    anand_param_go_mat = 1
    anand_param_eta_cv_mat = 1
    anand_param_p_mat = 1
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [porepressure]
  []
[]

[BCs]
  [roller_xmin]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'
  []
  [roller_ymin]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'
  []
  [plane_strain]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'back front'
  []
  [xmax_drained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = right
  []
  [top_velocity]
    type = FunctionDirichletBC
    variable = disp_y
    function = top_velocity
    boundary = top
  []
[]

[Functions]
  [top_velocity]
    type = PiecewiseLinear
    x = '0 0.002 0.006 0.014 0.03 0.046 0.062 0.078 0.094 0.11 0.126 0.142 0.158 0.174 0.19 0.206 0.222 0.238 0.254 0.27 0.286 0.302 0.318 0.334 0.35 0.366 0.382 0.398 0.414 0.43 0.446 0.462 0.478 0.494 0.51 0.526 0.542 0.558 0.574 0.59 0.606 0.622 0.638 0.654 0.67 0.686 0.702'
    y = '-0.0002091242 -0.0002136514 -0.0002170636 -0.0002214434 -0.0002275459 -0.0002322983 -0.0002363412 -0.0002398737 -0.0002429855 -0.0002457335 -0.0002481619 -0.0002503085 -0.0002522060 -0.0002538834 -0.0002553662 -0.0002566770 -0.0002578358 -0.0002588601 -0.0002597656 -0.0002605661 -0.0002612738 -0.0002618993 -0.0002624523 -0.0002629412 -0.0002633733 -0.0002637553 -0.0002640930 -0.0002643916 -0.0002646555 -0.0002648888 -0.0002650950 -0.0002652773 -0.0002654385 -0.0002655809 -0.0002657069 -0.0002658182 -0.0002659166 -0.0002660036 -0.0002660805 -0.0002661485 -0.0002662086 -0.0002662618 -0.0002663087 -0.0002663500 -0.0002663869 -0.0002664194 -0.0002664481'
  []
[]

[ICs]
  # Initialize all auxiliary variables properly
  [alpha_damage_ic]
    type = ConstantIC
    variable = alpha_damagedvar_aux
    value = 0.0
  []
  [B_damage_ic]
    type = ConstantIC
    variable = B_damagedvar_aux
    value = 0.0
  []
  [vel_x_ic]
    type = ConstantIC
    variable = vel_x
    value = 0.0
  []
  [vel_y_ic]
    type = ConstantIC
    variable = vel_y
    value = 0.0
  []
  [vel_z_ic]
    type = ConstantIC
    variable = vel_z
    value = 0.0
  []
  [struct_stress_ic]
    type = ConstantIC
    variable = structural_stress_coefficient_aux
    value = 1.0
  []
  [cg_aux_ic]
    type = ConstantIC
    variable = cg_aux
    value = 1e-12
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
  [vel_z]
    order = FIRST
    family = LAGRANGE
  []
  [alpha_damagedvar_aux]
    order = FIRST
    family = LAGRANGE
  []
  [B_damagedvar_aux]
    order = FIRST
    family = LAGRANGE
  []
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
  [gradx_alpha_damagedvar]
    order = CONSTANT
    family = MONOMIAL
  []
  [grady_alpha_damagedvar]
    order = CONSTANT
    family = MONOMIAL
  []
  [cg_aux]
    order = FIRST
    family = LAGRANGE
  []
  [nonlocal_xi]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  
  # Keep damage variables at zero
  [alpha_damage_zero]
    type = ConstantAux
    variable = alpha_damagedvar_aux
    value = 0.0
  []
  [B_damage_zero]
    type = ConstantAux
    variable = B_damagedvar_aux
    value = 0.0
  []
  [struct_stress_one]
    type = ConstantAux
    variable = structural_stress_coefficient_aux
    value = 1.0
  []
  
  # Keep velocities at zero
  [vel_x_zero]
    type = ConstantAux
    variable = vel_x
    value = 0.0
  []
  [vel_y_zero]
    type = ConstantAux
    variable = vel_y
    value = 0.0
  []
  [vel_z_zero]
    type = ConstantAux
    variable = vel_z
    value = 0.0
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
  [../]
[]

[Materials]
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        # outputs = exodus
    []
    #
    [damage_mat]
        type = DiffusedDamageBreakageMaterialMainApp
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        structural_stress_coefficient = structural_stress_coefficient_aux
        #build L matrix using velocity
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Diffused
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain strain_invariant_ratio'
        outputs = exodus
    []

    #shear stress perturbation
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'shear_stress_perturbation damage_perturbation'
        prop_values = '0.0 0.0'
    []
[] 


[Postprocessors]
  [dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = if(0.15*t<0.01,0.15*t,0.01)
  []
[]




[Preconditioning]
  [basic]
    type = SMP
    full = true
  []
[]


[Executioner]
  type = Transient
  solve_type = PJFNK
  start_time = 0
  end_time = 0.7
  [TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 0.001
  []
[]

[Outputs]
  exodus = true
[]



