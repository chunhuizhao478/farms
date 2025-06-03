#solid properties
#-------------------------------------------------#
# solid_density = 2640
youngs_modulus = 48.5e9
poissons_ratio = 0.22
solid_bulk_compliance = 3.46e-11
lambda_o = ${fparse youngs_modulus*poissons_ratio/(1+poissons_ratio)/(1-2*poissons_ratio)}
shear_modulus_o = ${fparse youngs_modulus/(2*(1+poissons_ratio))}
#-------------------------------------------------#

#damage-breakage properties
#-------------------------------------------------#
xi_o = -0.8073
xi_d = -0.8073
Cg = 1e-12
#-------------------------------------------------#

#fluid properties
#-------------------------------------------------#
porosity = 0.008
permeability = '1E-20 0 0 0 1E-20 0 0 0 1E-20'
fluid_density = 1000
viscosity = 1e-3
biot_coefficient = 0.5
fluid_bulk_modulus = 2.2e+9
#-------------------------------------------------#

#boundary conditions
#-------------------------------------------------#
outer_confinement_pressure = 20.6e6
inner_confinement_pressure = 3.4e6
#-------------------------------------------------#

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/mesh_adaptive_test.msh'
    []
    displacements = 'disp_x disp_y disp_z'
[]
    
[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
    PorousFlowDictator = dictator #All porous modules must contain 

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
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
[]
  
[Variables]
    #displacement components
    [disp_x]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    [disp_z]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    #pore pressure
    [pp]
        order = FIRST
        family = LAGRANGE  
        initial_condition = 3.4e6
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []  
    #
    [./biot_modulus_aux]
        order = CONSTANT
        family = MONOMIAL      
    []   
    [./effective_perm00_aux]
        order = CONSTANT
        family = MONOMIAL      
    [] 
    [./effective_perm11_aux]
        order = CONSTANT
        family = MONOMIAL      
    []
    [./effective_perm22_aux]
        order = CONSTANT
        family = MONOMIAL      
    [] 
    #
    [initial_I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_xi_aux]
        order = FIRST
        family = MONOMIAL
    [] 
    #
    [alpha_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.0
    []
    [B_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.0
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
    #pressure coupling on stress tensor
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
    [poro_z]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = ${biot_coefficient}
        variable = disp_z
        component = 2
    []
    #flux * grad(test)
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        multiply_by_density = false
        variable = pp
        gravity = '0 0 0'
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
        type = FunctionDirichletBC
        variable = disp_z
        boundary = 6
        function = applied_load_top
    []
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = 4
          factor = ${outer_confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
        [../]
        [./inner_boundary]
          boundary = 5
          factor = ${inner_confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
        [../]
    []
    #inner pressure
    # [inner_pressure]
    #     type = PorousFlowPiecewiseLinearSink
    #     variable = pp
    #     boundary = '5'
    #     pt_vals = '-1e9 1e9' # x coordinates defining g
    #     multipliers = '-1e9 1e9' # y coordinates defining g
    #     PT_shift = ${inner_confinement_pressure}
    #     flux_function = 1E-6 # Variable C
    #     fluid_phase = 0
    # []
[]

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-2.6477e-5'
    []
[]

[AuxKernels]
    [get_initial_I2]
        type = MaterialRealAux
        variable = initial_I2_aux
        property = I2_initial
        block = '1 2'
    []
    [get_initial_xi]
        type = MaterialRealAux
        variable = initial_xi_aux
        property = xi_initial
        block = '1 2'
    []
    [get_initial_I2_cdbm]
        type = MaterialRealAux
        variable = initial_I2_aux
        property = second_elastic_strain_invariant
        block = '3'
    []
    [get_initial_xi_cdbm]
        type = MaterialRealAux
        variable = initial_xi_aux
        property = strain_invariant_ratio
        block = '3'
    []
[]
    
[Materials]
    #shear stress perturbation
    [damage_perturbation]
        type = PerturbationRadial
        nucl_center = '0 0 0'
        peak_value = 0
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'shear_stress'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
        outputs = exodus
    []
    [getxi]
        type = ComputeXi
        outputs = exodus
        block = '1 2'
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBMDiffused
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        output_properties = 'stress elastic_strain_tensor plastic_strain_tensor total_strain_tensor strain_invariant_ratio'
        outputs = exodus
        block = '3'
    []
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = ${youngs_modulus}
      poissons_ratio = ${poissons_ratio}
    [../]
    [./elastic_stress]
      type = ComputeLinearElasticStress
      outputs = exodus
      block = '1 2'
    []
    [strain]
      type = ComputeSmallStrain
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
      solid_bulk_compliance = ${solid_bulk_compliance}
      fluid_bulk_modulus = ${fluid_bulk_modulus}
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
  
#this user object must contain for porous flow
[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pp'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[]

[Executioner]
    type = Steady
    solve_type = PJFNK
    # solve_type = Newton
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-10
    nl_abs_tol = 1e-12
    l_tol = 1e-5
[]
  
[Outputs]
    exodus = true
[]