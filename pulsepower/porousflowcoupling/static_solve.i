confinement_pressure  = 1e6
initial_pore_pressure = 0.0965e6

[Mesh]
    [./msh]
      type = FileMeshGenerator
      file =  '../meshfile/debug_sample.msh'
    []
    [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = '-0.00570169 -0.00612895 0'
      new_boundary = corner_ptr
      input = msh
      use_closest_node=true
    []
    displacements = 'disp_x disp_y disp_z'
  []
  
  [GlobalParams]
    displacements = 'disp_x disp_y disp_z'
    PorousFlowDictator = dictator #All porous modules must contain
  []

  [Variables]
    #displacement components
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
    #pore pressure
    [pp]
        order = FIRST
        family = LAGRANGE 
    []
  []

  [Kernels]
    #effective stress tensor
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        use_displaced_mesh = false
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        use_displaced_mesh = false
    []
    #effective pressure coupling on stress tensor: _pf * biot_coefficient
    #this effective pressure _pf = saturation * _pp
    #it is declared in "PorousFlowEffectiveFluidPressure"
    [poro_x]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = 0.7
      variable = disp_x
      component = 0
    []
    [poro_y]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = 0.7
      variable = disp_y
      component = 1
    []
    [poro_z]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = 0.7
      component = 2
      variable = disp_z
    []
    #flux * grad(test)
    [flux]
      type = PorousFlowFullySaturatedDarcyBase
      variable = pp
      gravity = '0 0 0'
    []
  []
  
  [BCs]
    [./Pressure]
        #assign pressure on top surface
        [static_pressure_top]
          boundary = 2
          factor = ${confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
          use_displaced_mesh = false
        []
        #assign pressure on outer surface
        [static_pressure_outer]
          boundary = 3
          factor = ${confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
          use_displaced_mesh = false
        []  
        [static_pressure_bottom]
          boundary = 5
          factor = ${confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
          use_displaced_mesh = false
        []              
    []
    # fix ptr
    [./fix_cptr1_x]
      type = DirichletBC
      variable = disp_x
      boundary = corner_ptr
      value = 0
    []
    [./fix_cptr2_y]
      type = DirichletBC
      variable = disp_y
      boundary = corner_ptr
      value = 0
    []
    [./fix_cptr3_z]
      type = DirichletBC
      variable = disp_z
      boundary = corner_ptr
      value = 0
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 50e9
      poissons_ratio = 0.373
    [../]
    [./elastic_stress]
      type = ComputeFiniteStrainElasticStress
      outputs = exodus
    []
    [strain]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y disp_z'
      outputs = exodus
    []
    #
    [temperature]
      type = PorousFlowTemperature
    []
    [eff_fluid_pressure_qp]
      type = PorousFlowEffectiveFluidPressure
    []
    #compute volumetric strain and its rate
    [vol_strain]
      type = PorousFlowVolumetricStrain
      outputs = exodus
    []
    #This Material is used for the fully saturated single-phase situation "
    #"where porepressure is the primary variable", saturation = 1.0
    [ppss]
      type = PorousFlow1PhaseFullySaturated
      porepressure = pp
    []
    #List of variables that represent the mass fractions.
    #If no "variables are provided then num_phases=1=num_components."
    [massfrac]
      type = PorousFlowMassFraction
    []
    #compute porosity
    [porosity]
      type = PorousFlowPorosityConst # only the initial value of this is ever used
      porosity = 0.008
    []
    #comopute permeability
    [permeability]
      type = PorousFlowPermeabilityConst
      permeability = '5E-19 0 0 0 5E-19 0 0 0 5E-19'
    []
    #compute biot modulus
    [biot_modulus]
      type = PorousFlowConstantBiotModulus
      biot_coefficient = 0.7
      solid_bulk_compliance = 1.524e-11 #checked
      fluid_bulk_modulus = 2.24e+9
    []
    #Compute density and viscosity
    [simple_fluid_qp]
      type = PorousFlowSingleComponentFluid
      fp = the_simple_fluid
      phase = 0
    []
  []

  #provide fluid properties for porous flow 
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2.24e9
      density0 = 1000
      thermal_expansion = 0
      viscosity = 1e-3
    []
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []

  # this user object must contain for porous flow
  [UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'pp'
      number_fluid_phases = 1
      number_fluid_components = 1
    []
  []

  [ICs]
    [disp_x_ic]
      type = ConstantIC
      variable = disp_x
      value = 0
    []
    [disp_y_ic]
      type = ConstantIC
      variable = disp_y
      value = 0
    []
    [disp_z_ic]
      type = ConstantIC
      variable = disp_z
      value = 0
    []
    [pp_ic]
      type = ConstantIC
      variable = pp
      value = ${initial_pore_pressure}
    []
  []
  
  [Executioner]
    type = Steady
    solve_type = Newton
    petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    petsc_options_value = '101                asm      lu'
  
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    automatic_scaling = true
  []
  
  [Outputs]
    exodus = true
  []