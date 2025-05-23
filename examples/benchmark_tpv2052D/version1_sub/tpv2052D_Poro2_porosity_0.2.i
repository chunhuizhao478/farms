# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 200
        ny = 200
        xmin = -15000
        xmax = 15000
        ymin = -15000
        ymax = 15000
        elem_type = QUAD9
        
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y>0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
        add_interface_on_two_sides = true
    []   
[]

[GlobalParams]
    displacements = 'disp_x disp_y' 
    q = 0.1
[]

[Variables]
    [./disp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_vel_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_vel_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    [./vel_x]
        order = SECOND
        family = LAGRANGE
    []
    [./vel_y]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_x]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_y]
        order = SECOND
        family = LAGRANGE
    []
    [./resid_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./traction_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./traction_y]
        order = SECOND
        family = LAGRANGE
    [../]
[]

[AuxKernels]
    [velocity_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
    []
    [velocity_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
    [] 
    #obtain system residuals by tagging
  #  [restore_x]
  #      type = TagVectorAux
 #       vector_tag = 'restore_tag_x'
  #      v = 'disp_x'
 #       variable = 'resid_x'
 #       execute_on = 'TIMESTEP_END'
 #   []
   # [restore_y]
  #      type = TagVectorAux
  ##      vector_tag = 'restore_tag_y'
  #      v = 'disp_y'
 #       variable = 'resid_y'
  #      execute_on = 'TIMESTEP_END'
  #  []
 #   [restore_pressurex]
 #       type = TagVectorAux
 #       vector_tag = 'restore_pressx_tag'
 #       v = 'disp_x'
  #      variable = 'resid_pressure_x'
  #      execute_on = 'TIMESTEP_END'
  #  []
 #   [restore_pressurey]
 #       type = TagVectorAux
  #      vector_tag = 'restore_pressy_tag'
  #      v = 'disp_y'
  #      variable = 'resid_pressure_y'
  #      execute_on = 'TIMESTEP_END'
 #   []
  #  [restore_dampx]
  #      type = TagVectorAux
 #       vector_tag = 'restore_dampx_tag'
 #       v = 'disp_x'
 #       variable = 'resid_damp_x'
 #       execute_on = 'TIMESTEP_END'
  #  []
   # [restore_dampy]
  #      type = TagVectorAux
  #      vector_tag = 'restore_dampy_tag'
  #      v = 'disp_y'
  #      variable = 'resid_damp_y'
  #      execute_on = 'TIMESTEP_END'
  #  []
[]

[Problem]
#    extra_tag_vectors = 'restore_tag_x restore_tag_y restore_pressx_tag restore_pressy_tag restore_dampx_tag restore_dampy_tag'
[]

[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false   
      #  extra_vector_tags = 'restore_tag_x' 
        save_in = 'resid_x'
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
      #   extra_vector_tags = 'restore_tag_y' 
      save_in = 'resid_y'
    [../]
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_y
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_x]
        type = CoupledFluidInertialForce
        variable = disp_x
        fluid_vel = fluid_vel_x
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_y]
        type = CoupledFluidInertialForce
        variable = disp_y
        fluid_vel = fluid_vel_y
        use_displaced_mesh = false
    [../]
    [./darcyflow_x]
        type = DynamicDarcyFlow2
        variable = fluid_vel_x
        skeleton_acceleration = disp_x
    [../]
    [./darcyflow_y]
        type = DynamicDarcyFlow2
        variable = fluid_vel_y
        skeleton_acceleration = disp_y
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_x
        porepressure = p
        component = 0
       #  extra_vector_tags = 'restore_pressx_tag' 
       save_in = 'resid_pressure_x'
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
       #  extra_vector_tags = 'restore_pressy_tag' 
       save_in = 'resid_pressure_y'
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_x
        porepressure = p
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_y
        porepressure = p
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p
       #  displacements = 'disp_x disp_y'
       
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p
        u = fluid_vel_x
        v = fluid_vel_y
        pressure = p
    [../]
    [./Reactionx]
        type = StiffPropDamping
        variable = 'disp_x'
        component = '0'
       #  extra_vector_tags = 'restore_dampx_tag' 
       save_in = 'resid_damp_x'
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = 'disp_y'
        component = '1'
       # extra_vector_tags = 'restore_dampy_tag' 
       save_in = 'resid_damp_y'
    []
[]

[Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 15.46e9
        shear_modulus = 13.86e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [Strain]
        type = ComputeSmallStrain
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2320
    []
    [./rhof]
        type = GenericConstantMaterial
        prop_names = rhof
        prop_values = 1000
    [../]
    [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 2.24
    [../]
    [./porosity]
        type = GenericConstantMaterial
        prop_names = porosity
        prop_values = 0.2
    [../]
    [./hydconductivity]
        type = GenericConstantMaterial
        prop_names = hydconductivity
        prop_values = 1.1280533319e-9
    [../]
    [./hydconductivity_layer]
        type = GenericConstantMaterial
        prop_names = hydconductivity_layer
        prop_values = 1.1280533319e-9
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.567
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 1.0084e10
    [../]
    [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
    [../]
[]

[BCs]
    [./trac_fault_main_x]
        type = CoupledVarNeumannBC
        variable = disp_x
        boundary = 'Block0_Block1'
        v = traction_x
        coef = -1 
    [../]
    [./trac_fault_main_y]
        type = CoupledVarNeumannBC
        variable = disp_y
        boundary = 'Block0_Block1'
        v = traction_y
        coef = -1 
    [../]
    [./trac_fault_sec_x]
        type = CoupledVarNeumannBC
        variable = disp_x
        boundary = 'Block1_Block0'
        v = traction_x
       #coef = -1 
    [../]
    [./trac_fault_sec_y]
        type = CoupledVarNeumannBC
        variable = disp_y
        boundary = 'Block1_Block0'
        v = traction_y
       # coef = -1 
    [../]
    [./fault_p]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = Block0_Block1
        value = 0.0
    [../]
    [./fault_n]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = Block1_Block0
        value = 0.0
    [../]
    [./left_vf]
        type = DirichletBC
        variable = fluid_vel_x
        boundary = left
        value = 0.0
    [../]
    [./right_vf]
        type = DirichletBC
        variable = fluid_vel_x
        boundary = right
        value = 0.0
    [../]
    [./top_vf]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = top
        value = 0.0
    [../]
    [./bot_vf]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = bottom
        value = 0.0
    [../]
  ##non-reflecting bc
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = right
    []
[]


[UserObjects]
   #evalute system residual after system solve before auxkernels
 #    [recompute_residual_tag_x]
 #        type = ResidualEvaluationUserObject
 #        vector_tag = 'restore_tag_x'
 #        force_preaux = true
  #       execute_on = 'TIMESTEP_END'
  #   []
  #   [recompute_residual_tag_y]
  #       type = ResidualEvaluationUserObject
  #       vector_tag = 'restore_tag_y'
  #       force_preaux = true
  #       execute_on = 'TIMESTEP_END'
  #   []
 ##     [recompute_residual_pressure_tag_x]
 #        type = ResidualEvaluationUserObject
 #        vector_tag = 'restore_pressx_tag'
 #        force_preaux = true
 #        execute_on = 'TIMESTEP_END'
 #    []
 #    [recompute_residual_pressure_tag_y]
#         type = ResidualEvaluationUserObject
 #        vector_tag = 'restore_pressy_tag'
 #        force_preaux = true
 #        execute_on = 'TIMESTEP_END'
 #    []
 #    [recompute_residual_tag_dampx]
  #       type = ResidualEvaluationUserObject
 #        vector_tag = 'restore_dampx_tag'
 #        force_preaux = true
 #        execute_on = 'TIMESTEP_END'
 #    []
 #    [recompute_residual_tag_dampy]
 #        type = ResidualEvaluationUserObject
  #       vector_tag = 'restore_dampy_tag'
  #       force_preaux = true
 #        execute_on = 'TIMESTEP_END'
  #   []

[]


[Preconditioning]
    [./smp]
        type = SMP
        full = true
       # petsc_options_iname = '-snes_test_jacobian  -pc_factor_mat_solver_type -pc_factor_shift_type'
        #petsc_options_value = 'lu umfpack NONZERO'
    [../]
[]

[Executioner]
    type = Transient
    dt = 0.002
    end_time = 3.6
    automatic_scaling = true
    [TimeIntegrator]
         type = CentralDifference
    []
    
[]

[Outputs]
    exodus = true
    time_step_interval = 20
     #print_linear_residuals = true
[]

[MultiApps]
    #allocate transfer from mainApp to subApp
    [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'sub.i'
      execute_on = 'INITIAL TIMESTEP_END'
    [../]
[]

[Transfers]
    #get displacement residuals from subApp to mainApp
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'traction_sub_x traction_sub_y'
        variable = 'traction_x traction_y'
        execute_on = 'INITIAL TIMESTEP_END'
    []
    #push system residual vector from mainApp to subApp
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'disp_x disp_y resid_x resid_y resid_pressure_x resid_pressure_y resid_damp_x resid_damp_y  p'
        variable = 'disp_czm_x disp_czm_y resid_sub_x resid_sub_y resid_pressure_sub_x resid_pressure_sub_y resid_damp_sub_x resid_damp_sub_y p_sub_main'
        execute_on = 'INITIAL TIMESTEP_END'
    []
[]

