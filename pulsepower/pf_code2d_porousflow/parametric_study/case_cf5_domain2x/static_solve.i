confinement_pressure  = 5000000.0
initial_pore_pressure = 96500.0

#solid properties
#----------------------------------------------------#
E = 50e9
nu = 0.3
K = '${fparse E/3.0/(1.0-2.0*nu)}' #bulk modulus of porous material
K_s = 50e9 #bulk modulus of solid grains
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
fluid_density = 1000
biot_coefficient = ${fparse 1 - K/K_s}
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-3
porosity = 0.008
solid_bulk_modulus_compliance = ${fparse 1.0/K} #this is the bulk modulus of the porous medium
permeability = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
#----------------------------------------------------#

#initial damage box 1
bottom_left1 = '-0.0025 -2e-4 0'
top_right1 = '0.0025 2e-4 0'

#initial damage box 2
bottom_left2 = '-2e-4 -0.0025 0'
top_right2 = '2e-4 0.0025 0'

[Mesh]
    [./msh]
      type = FileMeshGenerator
      file =  '../../../2dmeshfile/fieldscale_test1_2d_extend2x.msh'
    []
    [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = '0.02 0.02 0'
      new_boundary = corner_ptr
      input = msh
      use_closest_node=true
    []
    [./subdomain_id]
      type = SubdomainBoundingBoxGenerator
      bottom_left = ${bottom_left1}
      top_right = ${top_right1}
      location = INSIDE
      block_id = 1
      input = extranodeset1
    []
    [./subdomain_id2]
      type = SubdomainBoundingBoxGenerator
      bottom_left = ${bottom_left2}
      top_right = ${top_right2}
      location = INSIDE
      block_id = 1
      input = subdomain_id
    []
    displacements = 'disp_x disp_y'
  []

  [GlobalParams]
    displacements = 'disp_x disp_y'
    PorousFlowDictator = dictator #All porous modules must contain
  []

  [Variables]
    #displacement components
    [disp_x]
        order = FIRST
        family = LAGRANGE
        scaling = 1e-6
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
        scaling = 1e-6
    []
    #pore pressure
    [pp]
        order = FIRST
        family = LAGRANGE
    []
  []

  [AuxVariables]
    #reaction force
    [fx]
    []
    [fy]
    []
    [fz]
    []
  []

  [Kernels]
    #effective stress tensor
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = false
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
    #effective pressure coupling on stress tensor: _pf * biot_coefficient
    #this effective pressure _pf = saturation * _pp
    #it is declared in "PorousFlowEffectiveFluidPressure"
    [poro_x]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_x
      component = 0
    []
    [poro_y]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_y
      component = 1
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
        #assign pressure on outer surface
        [static_pressure_outer]
          boundary = 1
          factor = ${confinement_pressure}
          displacements = 'disp_x disp_y'
          use_displaced_mesh = false
          save_in_disp_x = fx
          save_in_disp_y = fy
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
  []

  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = ${E}
      poissons_ratio = ${nu}
    [../]
    [./elastic_stress]
      type = ComputeLinearElasticStress
      outputs = exodus
    []
    [strain]
      type = ComputeSmallStrain
      displacements = 'disp_x disp_y'
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
      porosity = ${porosity}
    []
    #comopute permeability
    [permeability]
      type = PorousFlowPermeabilityConst
      permeability = ${permeability}
    []
    #compute biot modulus
    [biot_modulus]
      type = PorousFlowConstantBiotModulus
      biot_coefficient = ${biot_coefficient}
      solid_bulk_compliance = ${solid_bulk_modulus_compliance}
      fluid_bulk_modulus = ${fluid_bulk_modulus}
    []
    #Compute density and viscosity
    [simple_fluid_qp]
      type = PorousFlowSingleComponentFluid
      fp = the_simple_fluid
      phase = 0
    []

    # Fluid driving energy density psi_f
    [fluid_driving_energy]
      type = ElkPorousFlowFluidDrivingEnergy
      biot_coefficient = ${biot_coefficient}
    []
  []

  #provide fluid properties for porous flow
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = ${fluid_bulk_modulus}
      density0 = ${fluid_density}
      thermal_expansion = 0
      viscosity = ${viscosity}
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
      porous_flow_vars = 'pp disp_x disp_y'
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
    [pp_ic]
      type = ConstantIC
      variable = pp
      value = ${initial_pore_pressure}
    []
  []

  [Executioner]
    type = Steady
    solve_type = Newton
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'

    line_search = 'bt'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-10
    l_tol = 1e-5
    # automatic_scaling = true
  []

  [Outputs]
    exodus = true
    [csv]
      type = CSV
      execute_on = 'initial timestep_end'
      time_step_interval = 1
      show = 'full_energy_static full_input_energy_static  solid_elastic_energy_static fluid_elastic_energy_total_static'
    []
  []

###############################Energy Calculation##############################
# input energy
###############################################################################
[Postprocessors]
  [external_work]
    type = FarmsExternalWork
    boundary = '1'
    forces = 'fx fy fz'
    use_displacement_work = true
  []
[]

[Postprocessors]
  [full_input_energy_static]
      type = ParsedPostprocessor
      expression = '-1 * external_work'
      pp_names = 'external_work'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# solid elastic energy
###############################################################################
[AuxVariables]
  [solid_elastic_energy]
      order = FIRST
      family = MONOMIAL
  []
[]

[AuxKernels]
  [solid_elastic_energy]
    type = ElasticEnergyAux
    variable = solid_elastic_energy
  []
[]

[Postprocessors]
  [solid_elastic_energy_static]
    type = ElementIntegralVariablePostprocessor
    variable = solid_elastic_energy
  []
[]

[Postprocessors]
  [solid_elastic_energy_total_static]
      type = ParsedPostprocessor
      expression = 'solid_elastic_energy_static'
      pp_names = 'solid_elastic_energy_static'
      execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# fluid elastic energy
###############################################################################
# [Postprocessors]
#   [fluid_elastic_energy_static]
#     type = ElementIntegralMaterialProperty
#     mat_prop = fluid_driving_energy_density
#   []
# []

# [Postprocessors]
#   [fluid_elastic_energy_total_static]
#       type = ParsedPostprocessor
#       expression = 'fluid_elastic_energy_static'
#       pp_names = 'fluid_elastic_energy_static'
#       execute_on = 'INITIAL TIMESTEP_END'
#   []
# []

[AuxVariables]
  [fluid_elastic_energy]
      order = CONSTANT
      family = MONOMIAL
  []
[]

[AuxKernels]
  [fluid_elastic_energy]
      type = ParsedAux
      variable = fluid_elastic_energy
      coupled_variables = 'elastic_strain_00 elastic_strain_11 elastic_strain_22 pp'
      expression = "0.5 * ${biot_coefficient} * -pp * (elastic_strain_00+elastic_strain_11+elastic_strain_22)"
  []
[]

[Postprocessors]
  [fluid_elastic_energy_total_static]
      type = ElementIntegralVariablePostprocessor
      variable = fluid_elastic_energy
  []
[]


# Full Energy
###############################################################################
[Postprocessors]
  [full_energy_static]
    type = ParsedPostprocessor
    expression = 'solid_elastic_energy_total_static + fluid_elastic_energy_total_static'
    pp_names = 'solid_elastic_energy_total_static fluid_elastic_energy_total_static'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
