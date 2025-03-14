[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        xmin = -1000
        xmax = 1000
        ymin = -1000
        ymax = 1000
        nx = 10
        ny = 10
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
    out_of_plane_strain = strain_zz
[]

[Variables]
    [./disp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./strain_zz]
    []
[]

[Physics/SolidMechanics/QuasiStatic]
    [plane_stress]
      planar_formulation = WEAK_PLANE_STRESS
      strain = SMALL
      generate_output = 'stress_xx stress_yy stress_xy'
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
[]

[Materials]
    #take PMMA: poisson's ratio 0.37, shear modulus 2.17GPa
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e9
        shear_modulus = 1e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 1180 #kg/m^3
    []
[]
    
[Executioner]
    #dt = prefactor * dx / pressure wave speed = 0.0009 / 2985
    type = Transient
    dt = 1e-4 
    num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]
    
[Outputs]
    exodus = true
    time_step_interval = 1
[]
