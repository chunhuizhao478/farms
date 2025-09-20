# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #
# [Note]: This serves as a test file, to run the full problem, please extend the domain size by modifying nx, ny, xmin, xmax, ymin, ymax

#parameters

##mesh parameters

##element size
elem_size = 100 #!!! element size near the fault, need to be consistent with the mesh file

##main fault parameters
xmin_fault = -15000 #xmin of fault
xmax_fault = 15000 #xmax of fault

##-------------------------##
##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
Cs = '${fparse shear_modulus_o / density }' #shear wave speed
Cp = '${fparse (lambda_o + 2 * shear_modulus_o) / density }' #pressure wave speed
##-------------------------##

##Slip weakening parameters##
Dc = 0.4 #characteristic length (m)
q = 0.1 #damping ratio
mu_s = 0.677 #static friction coefficient
mu_d = 0.525 #dynamic friction coefficient
##-------------------------##

##CDB model parameters##
xi_0 = -0.8 #strain invariants ratio: onset of damage evolution
xi_d = -0.8 #strain invariants ratio: onset of breakage healing

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
use_strain_rate_dependent_Cd = true #use strain rate dependent Cd
m_exponent = 0.8 #strain rate dependent parameters
strain_rate_hat = 1e-4 #strain rate dependent parameters
cd_hat = 10 #strain rate dependent parameters
###

CdCb_multiplier = 100 #multiplier between Cd and Cb
CBH_constant = 0 #coefficient of healing for breakage evolution
C_1 = 0 #coefficient of healing for damage evolution
C_2 = 0.05 #coefficient of healing for damage evolution
beta_width = 0.05 #coefficient gives width of transitional region
C_g = 1e-10 #material parameter: compliance or fluidity of the fine grain granular material
m1 = 10 #coefficient of power law indexes
m2 = 1 #coefficient of power law indexes
chi = 0.8 #energy ratio
##-------------------------##

##nucleation parameters##
peak_shear_stress = 81.6e6 #peak shear stress for nucleation (Pa)
nucl_center = '0 0' #nucleation center (x z)
nucl_radius = 1500 #nucleation radius (m)
##-------------------------##

##model parameters##
dt = 0.001 #time step size

end_time = 100.0 #end time for simulation

# num_steps = 40 #end_time or num_steps only one of them is needed
exodus_time_step_interval = 40 #time step interval for output
sample_snapshots_time_step_interval = 400 #time step interval for sample snapshots output
# csv_time_step_interval = 2 #time step interval for csv output
checkpoint_time_step_interval = 40 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##------------------------------------------------------------------------##

[Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 2
      nx = 800
      ny = 400
      xmin = -40000
      xmax = 40000
      ymin = -20000
      ymax = 20000
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y>0 & x>${xmin_fault} & x<${xmax_fault}'
      block_id = 1
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'y<0 & x>${xmin_fault} & x<${xmax_fault}'
      block_id = 2
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block_2
      split_interface = true
      block_pairs = '1 2'
    []
  []

  [GlobalParams]
    ##------------slip weakening------------##
    #primary variables
    displacements = 'disp_x disp_y'
    #damping ratio
    q = ${q}
    #characteristic length (m)
    Dc = ${Dc}
    #static friction coefficient
    mu_s = ${mu_s}
    #dynamic friction coefficient
    mu_d = ${mu_d}
    #element edge length (m)
    len = ${elem_size}

    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = ${lambda_o}
    
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = ${shear_modulus_o}

    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = ${xi_0}

    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = ${xi_d}

    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8

    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = ${Cd_constant}

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = ${CdCb_multiplier}

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = ${CBH_constant}

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = ${C_1}

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = ${C_2}

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = ${beta_width}

    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = ${C_g}

    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = ${m1}

    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = ${m2}

    # energy ratio
    chi = ${chi}    

  []

  [AuxVariables]
    [./resid_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_x]
      order = FIRST
      family = LAGRANGE
    []
    [./vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    ###
    #output jump, jump rate, traction quantities
    [jump_x_aux]
      order = FIRST
      family = MONOMIAL
    []
    [jump_x_rate_aux]
      order = FIRST
      family = MONOMIAL
    []
    [traction_x_aux]
      order = FIRST
      family = MONOMIAL
    [] 
    [jump_y_aux]
      order = FIRST
      family = MONOMIAL
    []
    [jump_y_rate_aux]
      order = FIRST
      family = MONOMIAL
    []
    [traction_y_aux]
      order = FIRST
      family = MONOMIAL
    []
    ###
    #output CDB model properties
    [alpha_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [B_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    [] 
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    ###
    [deviatoric_strain_rate_aux]
      order = FIRST
      family = MONOMIAL
    []  
  []

  [Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
      boundary = 'Block1_Block2'
      strain = SMALL
      generate_output='traction_x traction_y jump_x jump_y'
    [../]
  []

  [Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          generate_output = 'stress_xx stress_yy stress_xy'
          extra_vector_tags = 'restore_tag'
        []
      []
    []
  []

  [Problem]
    extra_tag_vectors = 'restore_tag'
  []

  [AuxKernels]
    [Displacment_x]
      type = ProjectionAux
      variable = disp_slipweakening_x
      v = disp_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_y]
      type = ProjectionAux
      variable = disp_slipweakening_y
      v = disp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_x]
      type = CompVarRate
      variable = vel_slipweakening_x
      coupled = disp_x
      execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
      type = CompVarRate
      variable = vel_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_END'
    []
    [Residual_x]
      type = ProjectionAux
      variable = resid_slipweakening_x
      v = resid_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
      type = ProjectionAux
      variable = resid_slipweakening_y
      v = resid_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [restore_x]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_x'
      variable = 'resid_x'
    []
    [restore_y]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_y'
      variable = 'resid_y'
    []
    ### slip weakening strike direction
    [get_jump_x_aux]
      type = MaterialRealAux
      property = jump_x
      variable = jump_x_aux
      boundary = 'Block1_Block2'
      execute_on = 'TIMESTEP_END'
    []
    [get_jump_x_rate_aux]
      type = FDCompVarRate
      variable = jump_x_rate_aux
      coupled = jump_x
      execute_on = 'TIMESTEP_END'
      boundary = 'Block1_Block2'
    []
    [get_traction_x_aux]
      type = MaterialRealAux
      property = traction_x
      variable = traction_x_aux
      boundary = 'Block1_Block2'
      execute_on = 'TIMESTEP_END'
    []
    ### slip weakening normal direction
    [get_jump_y_aux]
      type = MaterialRealAux
      property = jump_y
      variable = jump_y_aux
      boundary = 'Block1_Block2'
      execute_on = 'TIMESTEP_END'
    []
    [get_jump_y_rate_aux]
      type = FDCompVarRate
      variable = jump_y_rate_aux
      coupled = jump_y
      execute_on = 'TIMESTEP_END'
      boundary = 'Block1_Block2'
    []
    [get_traction_y_aux]
      type = MaterialRealAux
      property = traction_y
      variable = traction_y_aux
      boundary = 'Block1_Block2'
      execute_on = 'TIMESTEP_END'
    []
    ### get CDB model properties
    [get_alpha_damagedvar]
      type = MaterialRealAux
      variable = alpha_damagedvar_aux
      property = alpha_damagedvar
      execute_on = 'TIMESTEP_END'
    []
    [get_B]
      type = MaterialRealAux
      variable = B_aux
      property = B
      execute_on = 'TIMESTEP_END'
    []
    [get_xi]
      type = MaterialRealAux
      variable = xi_aux
      property = xi
      execute_on = 'TIMESTEP_END'
    []
    [get_initial_damage]
      type = SolutionAux
      solution = init_sol_components
      variable = initial_damage_aux
      from_variable = 'alpha_damagedvar_aux'
      execute_on = initial
    [../]
    ###
    [get_deviatoric_strain_rate]
      type = MaterialRealAux
      variable = deviatoric_strain_rate_aux
      property = deviatoric_strain_rate
      execute_on = 'TIMESTEP_END'
    []
  []

  [Kernels]
    [./inertia_x]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_x
    []
    [./inertia_y]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_y
    []
    [./Reactionx]
      type = StiffPropDamping
      variable = 'disp_x'
      component = '0'
    []
    [./Reactiony]
      type = StiffPropDamping
      variable = 'disp_y'
      component = '1'
    []
  []

  [Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeDamageBreakageStress3DSlipWeakening
        output_properties = 'B alpha_damagedvar xi I1 I2 deviatoric_strain_rate'
        use_strain_rate_dependent_Cd = ${use_strain_rate_dependent_Cd}
        m_exponent = ${m_exponent}
        strain_rate_hat = ${strain_rate_hat}
        cd_hat = ${cd_hat}
        zero_Cd_below_threshold = true
        outputs = exodus
    []
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'initial_breakage damage_perturbation density'
        prop_values = '0 0 ${density}'
    []
    [initial_damage_surround]
      type = ParsedMaterial
      property_name = 'initial_damage'
      coupled_variables = initial_damage_aux
      expression = 'initial_damage_aux'
      outputs = exodus
    []
    [./czm_mat]
        type = SlipWeakeningFrictionczm2dCDBM
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        peak_shear_stress = ${peak_shear_stress}
        nucl_center = ${nucl_center}
        nucl_radius = ${nucl_radius}
        boundary = 'Block1_Block2'
    [../]
    [./static_initial_strain_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_strain_tensor
        tensor_functions = 'func_initial_strain_xx   func_initial_strain_xy      func_initial_strain_xz 
                            func_initial_strain_xy   func_initial_strain_yy      func_initial_strain_yz
                            func_initial_strain_xz   func_initial_strain_yz      func_initial_strain_zz'
        output_properties = 'static_initial_strain_tensor'
        outputs = exodus
    [../]
    [./static_initial_stress_tensor] #this is used in the SlipWeakeningFrictionczm3dCDBM
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
        output_properties = 'static_initial_stress_tensor'
        outputs = exodus
    [../]
  []

  [Functions]
    ###strain field###
    [./func_initial_strain_xx]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_00'
    []
    [./func_initial_strain_xy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_01'
    []
    [./func_initial_strain_xz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_02'
    []
    [./func_initial_strain_yy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_11'
    []
    [./func_initial_strain_yz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_12'
    []
    [./func_initial_strain_zz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'elastic_strain_22'
    []
    ###stress field###
    [./func_initial_stress_xx]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_00'
    []
    [./func_initial_stress_xy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_01'
    []
    [./func_initial_stress_xz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_02'
    []
    [./func_initial_stress_yy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_11'
    []
    [./func_initial_stress_yz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_12'
    []
    [./func_initial_stress_zz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = 'stress_22'
    []
  []

  [UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_readcsv_out.e'
      system_variables = 'elastic_strain_00 elastic_strain_01 elastic_strain_02
                          elastic_strain_11 elastic_strain_12 elastic_strain_22
                          stress_00 stress_01 stress_02 stress_11 stress_12 stress_22 alpha_damagedvar_aux'
      timestep = LATEST
      force_preaux = true
      execute_on = 'INITIAL'
    [../]
  []

  [Executioner]
    type = Transient
    dt = ${dt}
    end_time = ${end_time}
    # num_steps = ${num_steps}
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
      use_constant_mass = true
    []
  []

  [Outputs]
    [exodus]
      type = Exodus
      execute_on = 'timestep_end'
      show = 'vel_slipweakening_x vel_slipweakening_y disp_slipweakening_x disp_slipweakening_y  alpha_damagedvar_aux B_aux xi_aux stress_xx stress_yy stress_xy deviatoric_strain_rate_aux initial_damage_aux initial_damage'
      time_step_interval = ${exodus_time_step_interval}
    []
    [out]
      type = Checkpoint
      time_step_interval = ${checkpoint_time_step_interval}
      num_files = ${checkpoint_num_files}
    []
    [sample_snapshots]
      type = Exodus
      execute_on = 'timestep_end'
      show = 'vel_slipweakening_x vel_slipweakening_y  disp_slipweakening_x disp_slipweakening_y  alpha_damagedvar_aux B_aux xi_aux stress_xx stress_yy stress_xy deviatoric_strain_rate_aux'
      time_step_interval = ${sample_snapshots_time_step_interval}
    []
  []    

  [BCs]
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = right
    []
  []
