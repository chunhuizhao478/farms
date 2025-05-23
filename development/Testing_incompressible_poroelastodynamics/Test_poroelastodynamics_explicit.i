[Mesh]
type = FileMesh
file = test.msh
block_id = '8'
block_name = 'domain'
boundary_id = '9 10 11 12 13'
boundary_name = 'bottom right topright topleft left'
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
[./v_x]
order = SECOND
family = LAGRANGE
[../]
[./v_y]
order = SECOND
family = LAGRANGE
[../]
[./a_x]
order = SECOND
family = LAGRANGE
[../]
[./a_y]
order = SECOND
family = LAGRANGE
[../]
[./af_x]
order = SECOND
family = LAGRANGE
[../]
[./af_y]
order = SECOND
family = LAGRANGE
[../]
[]


[Kernels]
   [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false     
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
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
    [./darcyflofluid_vel_x]
        type = DynamicDarcyFlow2
        variable = fluid_vel_x
        skeleton_acceleration = disp_x
    [../]
    [./darcyflofluid_vel_y]
        type = DynamicDarcyFlow2
        variable = fluid_vel_y
        skeleton_acceleration = disp_y
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_x
        porepressure = p
        component = 0
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling
        variable = fluid_vel_x
        porepressure = p
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling
        variable = fluid_vel_y
        porepressure = p
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p
        displacements = 'disp_x disp_y'
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p
        u = fluid_vel_x
        v = fluid_vel_y
        pressure = p
    [../]
[]


[AuxKernels]
[./accel_x]
type = NewmarkAccelAux
variable = a_x
displacement = disp_x
velocity = v_x
beta = 0.25
execute_on = timestep_end
[../]
[./accel_y]
type = NewmarkAccelAux
variable = a_y
displacement = disp_y
velocity = v_y
beta = 0.25
execute_on = timestep_end
[../]
[./vel_x]
type = NewmarkVelAux
variable = v_x
acceleration = a_x
gamma = 0.5
execute_on = timestep_end
[../]
[./vel_y]
type = NewmarkVelAux
variable = v_y
acceleration = a_y
gamma = 0.5
execute_on = timestep_end
[../]
[./fluidaccel_x]
type = NewmarkPoreFluidAccelAux
variable = af_x
darcyvel = fluid_vel_x
gamma = 0.5
execute_on = timestep_end
[../]
[./fluidaccel_y]
type = NewmarkPoreFluidAccelAux
variable = af_y
darcyvel = fluid_vel_y
gamma = 0.5
execute_on = timestep_end
[../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 14.5e6
    poissons_ratio = 0.3
    block = 'domain'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    block = 'domain'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 'domain'
  [../]
  [./density]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = density
    prop_values = 1986
  [../]
  [./rhof]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = rhof
    prop_values = 1000
  [../]
  [./porosity]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = porosity
    prop_values = 0.42
  [../]
  [./hydconductivity]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = hydconductivity
    prop_values = 1.01937e-8
  [../]
  [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 1
  [../]
  [./biotcoeff]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = biot_coefficient
    prop_values = 1.0
  [../]
[./constants]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = 'rho mu'
    prop_values = '1  1'
[../]
[]

[Functions]
active = 'bc_func bc_func2'
[./bc_func]
type = ParsedFunction
#value = 'if(x<5.0,0.0,15.0)*if(t<0.1,10*t,1.0)'
expression = 'if(t<0.1,-10*t*15e3,-1.0*15e3)'
[../]
[./bc_func2]
type = ParsedFunction
expression = 0
[../]
[]


[BCs]
[./bottom_y]
type = FunctionDirichletBC
variable = disp_y
boundary = 'bottom'
function = bc_func2
[../]
[./topright_y]
type = FunctionNeumannBC
  variable = disp_y
  boundary = 'topright'
  function = bc_func
[../]
[./topleft_x]
type = FunctionNeumannBC
variable = disp_y
boundary = 'topleft'
function = bc_func2
[../]
[./left_x]
type = FunctionDirichletBC
variable = disp_x
boundary = 'left'
function = bc_func2
[../]
[./right_x]
type = FunctionDirichletBC
variable = disp_x
boundary = 'right'
function = bc_func2
[../]
[./fluidbottom_y]
type = FunctionDirichletBC 
variable = fluid_vel_y
boundary = 'bottom'
function = bc_func2
[../]
[./fluidleft_x]
type = FunctionDirichletBC 
variable = fluid_vel_x
boundary = 'left'
function = bc_func2
[../]
[./fluidright_x]
type = FunctionDirichletBC 
variable = fluid_vel_x
boundary = 'right'
function = bc_func2
[../]
[./fluidtopright_y]
type = FunctionDirichletBC 
variable = fluid_vel_y
boundary = 'topright'
function = bc_func2
[../]
[./porepressure]
type = FunctionDirichletBC
variable = p
boundary = 'topleft'
function = bc_func2
[../]
[]


[Preconditioning]
[./smp]
type = SMP
full = true
#petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -pc_factor_shift_type'
#petsc_options_value = 'lu umfpack NONZERO'
[../]
[]
[Postprocessors]
[./rightcornerdisp]
type = PointValue
point = '10 10 0'
variable = disp_y
[../]
[./leftcornerdisp]
type = PointValue
point = '0 10 0'
variable = disp_y
[../]
[./rightcornerdarcy]
type = PointValue
point = '10 10 0'
variable = fluid_vel_y
[../]
[./leftcornerdarcy]
type = PointValue
point = '0 10 0'
variable = fluid_vel_y
[../]
[]
[Executioner]
type = Transient
solve_type = 'PJFNK'
start_time = 0
end_time = 10.0
dtmax = 0.002
dtmin =  0.002
[./TimeStepper]
type = ConstantDT
dt =  0.002
[../]
    [TimeIntegrator]
        type = NewmarkBeta
    []
[]
[Outputs]
exodus = true
[]