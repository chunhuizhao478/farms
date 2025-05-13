[GlobalParams]
    order = FIRST
    family = LAGRANGE
    displacements = 'disp_x disp_y'
    out_of_plane_strain = strain_zz
  []
  
  [Mesh]
    [./square]
      type = GeneratedMeshGenerator
      dim = 2
      nx = 2
      ny = 2
    [../]
  []
  
  [Variables]
    [./disp_x]
    [../]
    [./disp_y]
    [../]
  
    [./strain_zz]
    [../]
  []
  
  [Kernels]
    [./disp_x]
      type = StressDivergenceTensors
      variable = disp_x
      component = 0
    [../]
    [./disp_y]
      type = StressDivergenceTensors
      variable = disp_y
      component = 1
    [../]

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
  
    [./solid_z]
      type = WeakPlaneStress
      variable = strain_zz
    [../]
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      poissons_ratio = 0.0
      youngs_modulus = 1
    [../]
    [./strain]
      type = ComputePlaneSmallStrain
    [../]
    [./stress]
      type = ComputeLinearElasticStress
    [../]
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 1180 #kg/m^3
    []
  []
  
  [Preconditioning]
    [./SMP]
      type = SMP
      full = true
    [../]
  []
  
  [Executioner]
    type = Transient
    num_steps = 1
    petsc_options = '-ksp_view_pmat'
    [TimeIntegrator]
        type = FarmsCentralDifference
        solve_type = CONSISTENT
    []
  []