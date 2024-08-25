[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/cdbm_tpv2053d_small.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0
                    0 0 -1
                    0 0 1'
        new_boundary = 'left right bottom top back front'
    []    
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '-8000  -8000  -8000;
                  8000  -8000  -8000;
                 -8000  -8000   8000;
                  8000  -8000   8000'
        new_boundary = corner_ptr
        input = sidesets
    [] 
[]
  
[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
  
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
      
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
  
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
  
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.9
  
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
  
    ##Compute gamma_damaged_r, xi_1
    #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
    #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
    #check struct_param.m 
  
    #coefficient of damage solid modulus
    gamma_damaged_r = 34.785e9
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.825
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.8
    a0 = 7.42e9
    a1 = -21.341e9
    a2 = 19.028e9
    a3 = -4.924e9

    #
    D = 0
    
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [vel_x]
    []  
    [vel_y]
    []
    [vel_z]
    []
[]

[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_z
      coupled = disp_z
      execute_on = 'TIMESTEP_END'
    []
[]
  
[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
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
    [./inertia_z]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_z
    []
[]
  
[Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeDamageBreakageStress3Doverstress
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi'
        outputs = exodus
    []
    [strain]
        type = ComputeSmallStrain
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2700
    []
    [nonADdensity]
        type = GenericConstantMaterial
        prop_names = nonADdensity
        prop_values = 2700
    []
    [./static_initial_stress_tensor]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_stress_xx     func_stress_xy      func_stress_xz 
                            func_stress_xy     func_stress_yy      func_stress_yz
                            func_stress_xz     func_stress_yz      func_stress_zz'
        # outputs = exodus
    [../]
    # [./initial_damage_B1]
    #     type = GenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     block = '1'
    #     prop_values = '0.8'
    #     outputs = exodus
    # [] 
    # [./initial_damage_B2]
    #     type = GenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     block = '2'
    #     prop_values = '0'
    #     outputs = exodus
    # [] 
    # [./initial_damage_B3]
    #     type = GenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     block = '3'
    #     prop_values = '0.7'
    #     outputs = exodus
    # [] 
    # [./initial_damage_B4]
    #     type = GenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     block = '4'
    #     prop_values = '0'
    #     outputs = exodus
    # [] 
    # [./initial_damage_B5]
    #     type = GenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     block = '5'
    #     prop_values = '0'
    #     outputs = exodus
    # [] 
    [initial_damage]
        type = InitialDamageBenchmark
        nucl_center = '0 -2500 0'
        fault_plane = '-5000 5000 -5000 0 -500 500'
        nucl_distance = 400
        nucl_thickness = 400
        nucl_damage = 0.85
        e_damage = 0.6
        e_sigma = 1e3
        outputs = exodus
    []
    # [initial_damage]
    #     type = InitialDamageBenchmarkSperical
    #     nucl_center = '-2666.7 -4000 0'
    #     e_damage = 0.85
    #     e_sigma = 1e3
    #     outputs = exodus
    # []
[]
  
[Functions]
    [func_stress_xx]
        type = ConstantFunction
        value = -135e6
    [../]
    [func_stress_xy]
        type = ConstantFunction
        value = 0
    [../]
    [func_stress_xz]
        type = InitialShearStress3D
        nucl_center = '0 -2500 0'
        e_sigma = 1e3
        min_val = 60e6
        max_val = 85e6
    [../]
    [func_stress_yy]
        type = ConstantFunction
        value = -58.75e6
    [../]
    [func_stress_yz]
        type = ConstantFunction
        value = 0
    [../]
    [func_stress_zz]
        type = ConstantFunction
        value = -100e6
    [../]  
[]

#0.4/5773
[Executioner]
    type = Transient
    dt = 1e-4
    end_time = 10.0
    # num_steps = 10000
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]
  
[Outputs]
    exodus = true
    time_step_interval = 100
    [sample_snapshots]
        type = Exodus
        time_step_interval = 200
    []
    [snapshots]
        type = Exodus
        time_step_interval = 100
        overwrite = true
    []
[]