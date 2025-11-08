#parameters

##mesh parameters

##element size
elem_size = 100 #!!! element size near the fault, need to be consistent with the mesh file

##main fault parameters
xmin_fault = -15000 #xmin of fault
xmax_fault = 15000 #xmax of fault
ymin_domain = -2000 #ymin of domain
ymax_domain = 2000 #ymax of domain

##nonlocal averaging parameters
nonlocal_eqstrain_blocks = '100 200' #damage blocks
local_eqstrain_blocks = '1 2' #elastic block
nonlocal_eqstrain_blocks_1 = '100'
nonlocal_eqstrain_blocks_2 = '200'
nonlocal_averaging_length_scale = 200
nonlocal_averaging_radius = 400
##-------------------------##

##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
Cs = '${fparse sqrt(shear_modulus_o / density) }'
Cp = '${fparse sqrt((lambda_o + 2 * shear_modulus_o) / density) }'
##-------------------------##

##Slip weakening parameters##
Dc = 0.4 #characteristic length (m)
q = 0.2 #damping ratio
mu_s = 0.677 #static friction coefficient
mu_d = 0.525 #dynamic friction coefficient
##-------------------------##

##initial stress parameters##
sigmaxx = -120e6 #initial strike stress (Pa)
sigmayy = -120e6 #initial normal stress (Pa)
sigmaxy = 70e6 #initial shear stress (Pa)
##-------------------------##

##CDB model parameters##
xi_0 = -0.9 #strain invariants ratio: onset of damage evolution
xi_d = -0.9 #strain invariants ratio: onset of breakage healing
##-------------------------##

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
use_strain_rate_dependent_Cd = false #use strain rate dependent Cd
m_exponent = 0.8 #strain rate dependent parameters
strain_rate_hat = 5e-9 #strain rate dependent parameters
cd_hat = 10 #strain rate dependent parameters
##-------------------------##

CdCb_multiplier = 100 #multiplier between Cd and Cb
CBH_constant = 0 #coefficient of healing for breakage evolution
C_1 = 0 #coefficient of healing for damage evolution
C_2 = 0.05 #coefficient of healing for damage evolution
beta_width = 0.05 #coefficient gives width of transitional region
C_g = 1e-12 #material parameter: compliance or fluidity of the fine grain granular material
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
dt = 0.005 #time step size

end_time = 3600 #end time for simulation

# num_steps = 40 #end_time or num_steps only one of them is needed
exodus_time_step_interval = 40 #time step interval for output
# sample_snapshots_time_step_interval = 400 #time step interval for sample snapshots output
# csv_time_step_interval = 2 #time step interval for csv output
checkpoint_time_step_interval = 400 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##------------------------------------------------------------------------##

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/2dmesh_100m.msh'
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y>0 & x>${xmin_fault} & x<${xmax_fault} & y<${ymax_domain}'
      block_id = 100
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'y<0 & x>${xmin_fault} & x<${xmax_fault} & y>${ymin_domain}'
      block_id = 200
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block_2
      split_interface = true
      block_pairs = '100 200'
    []
    [./sidesets]
        input = split
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
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
    [eqstrain_nonlocal_aux]
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
      boundary = 'Block100_Block200'
      strain = SMALL
      generate_output='jump_x'
    [../]
  []

  [Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          #generate_output = 'stress_xx stress_yy stress_xy'
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
      boundary = 'Block100_Block200'
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
    [get_eqstrain_nonlocal]
      type = MaterialRealAux
      variable = eqstrain_nonlocal_aux
      property = eqstrain_nonlocal
      execute_on = 'TIMESTEP_END'
    []
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
    [stress_medium_nonlocal]
        type = ComputeDamageBreakageStress3DSlipWeakeningNonlocal
        output_properties = 'B alpha_damagedvar xi I1 I2 deviatoric_strain_rate'
        use_strain_rate_dependent_Cd = ${use_strain_rate_dependent_Cd}
        m_exponent = ${m_exponent}
        strain_rate_hat = ${strain_rate_hat}
        cd_hat = ${cd_hat}
        use_nonlocal_eqstrain = true
        nonlocal_eqstrain_blocks = ${nonlocal_eqstrain_blocks}
        zero_Cd_below_threshold = true
        static_solve_flag = false
        outputs = exodus
    []
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'eqstrain_nonlocal_initial initial_damage initial_breakage damage_perturbation density'
        prop_values = '0 0 0 0 ${density}'
    []
    #[initial_damage_surround]
    #  type = ParsedMaterial
    #  property_name = 'initial_damage'
    #  coupled_variables = initial_damage_aux
    #  expression = 'initial_damage_aux'
    #  #outputs = exodus
    #[]
    [./czm_mat]
        type = SlipWeakeningFrictionczm2dCDBM
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        peak_shear_stress = ${peak_shear_stress}
        nucl_center = ${nucl_center}
        nucl_radius = ${nucl_radius}
        boundary = 'Block100_Block200'
    [../]
    [./static_initial_strain_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_strain_tensor
        tensor_functions = '0   0      0
                            0   0      0
                            0   0      0'
        #output_properties = 'static_initial_strain_tensor'
        #outputs = exodus
    [../]
    [./static_initial_stress_tensor] #this is used in the SlipWeakeningFrictionczm3dCDBM
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
        #output_properties = 'static_initial_stress_tensor'
        #outputs = exodus
    [../]
    ##elastic domain
    #[elasticity]
    #    type = ComputeIsotropicElasticityTensor
    #    lambda = 32.04e9
    #    shear_modulus = 32.04e9
    #    use_displaced_mesh = false
    #[]
    #[stress]
    #    type = ComputeLinearElasticStress
    #    block = ${local_eqstrain_blocks}
    #[]
    #the ComputeDamageBreakageStress3DSlipWeakeningNonlocal takes old value for updating damage/breakage
    [nonlocal_eqstrain_block1]
        type = ElkNonlocalEqstrainUpdated
        average_UO = eqstrain_averaging_block1
        block = ${nonlocal_eqstrain_blocks_1}
    []
    [nonlocal_eqstrain_block2]
        type = ElkNonlocalEqstrainUpdated
        average_UO = eqstrain_averaging_block2
        block = ${nonlocal_eqstrain_blocks_2}
    []
    #for the block outside the region, nonlocal strain is equal to the local strain
    [nonlocal_eqstrain_block]
        type = ParsedMaterial
        property_name = eqstrain_nonlocal
        coupled_variables = 'xi_aux'
        expression = 'xi_aux'
        block = ${local_eqstrain_blocks}
    []
  []

[Functions]
  ###stress field###
  [./func_initial_stress_xx]
    type = ConstantFunction
    value = ${sigmaxx}
  []
  [./func_initial_stress_xy]
    type = ConstantFunction
    value = ${sigmaxy}
  []
  [./func_initial_stress_xz]
    type = ConstantFunction
    value = 0
  []
  [./func_initial_stress_yy]
    type = ConstantFunction
    value = ${sigmayy}
  []
  [./func_initial_stress_yz]
    type = ConstantFunction
    value = 0
  []
  [./func_initial_stress_zz]
    type = ConstantFunction
    value = 0
  []
[]

  [UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [eqstrain_averaging_block1]
        type = ElkRadialAverageUpdated
        length_scale = ${nonlocal_averaging_length_scale}
        prop_name = xi
        radius = ${nonlocal_averaging_radius}
        weights = BAZANT
        execute_on = TIMESTEP_END
        block = ${nonlocal_eqstrain_blocks_1}
    []
    [eqstrain_averaging_block2]
        type = ElkRadialAverageUpdated
        length_scale = ${nonlocal_averaging_length_scale}
        prop_name = xi
        radius = ${nonlocal_averaging_radius}
        weights = BAZANT
        execute_on = TIMESTEP_END
        block = ${nonlocal_eqstrain_blocks_2}
    []
  []

  [Executioner]
    type = Transient
    dt = ${dt}
    end_time = ${end_time}
    #num_steps = 10
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
      show = 'vel_slipweakening_x vel_slipweakening_y disp_slipweakening_x disp_slipweakening_y  alpha_damagedvar_aux B_aux xi_aux eqstrain_nonlocal_aux deviatoric_strain_rate_aux'
      time_step_interval = ${exodus_time_step_interval}
    []
    [out]
      type = Checkpoint
      time_step_interval = ${checkpoint_time_step_interval}
      num_files = ${checkpoint_num_files}
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
