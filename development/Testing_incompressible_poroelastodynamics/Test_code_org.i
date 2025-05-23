[Mesh]
type = FileMesh
file = test.msh
block_id = '8'
block_name = 'domain'
boundary_id = '9 10 11 12 13'
boundary_name = 'bottom right topright topleft left'
[]
[Variables]
[./u_x]
order = SECOND
family = LAGRANGE
[../]
[./u_y]
order = SECOND
family = LAGRANGE
[../]
[./w_x]
order = SECOND
family = LAGRANGE
[../]
[./w_y]
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
variable = u_x
component = 0
displacements = 'u_x u_y'
# use_displaced_mesh = false
[../]
[./stressdiv_y]
type = StressDivergenceTensors
variable = u_y
component = 1
displacements = 'u_x u_y'
# use_displaced_mesh = false
[../]
[./skeletoninertia_x]
type = InertialForce
variable = u_x
velocity = v_x
acceleration = a_x
beta = 0.25
gamma = 0.5
use_displaced_mesh = false
[../]
[./skeletoninertia_y]
type = InertialForce
variable = u_y
velocity = v_y
acceleration = a_y
beta = 0.25
gamma = 0.5
use_displaced_mesh = false
[../]
[./porefluidIFcoupling_x]
type = PoreFluidInertialForceCoupling
variable = u_x
fluidaccel = af_x
darcyvel = w_x
gamma = 0.5
[../]
[./porefluidIFcoupling_y]
type = PoreFluidInertialForceCoupling
variable = u_y
fluidaccel = af_y
darcyvel = w_y
gamma = 0.5
[../]
[./darcyflow_x]
type = DynamicDarcyFlow
variable = w_x
skeletondisp = u_x
skeletonvel = v_x
skeletonaccel = a_x
fluidaccel = af_x
gravity = 9.81
beta = 0.25
gamma = 0.5
[../]
[./darcyflow_y]
type = DynamicDarcyFlow
variable = w_y
skeletondisp = u_y
skeletonvel = v_y
skeletonaccel = a_y
fluidaccel = af_y
gravity = 9.81
beta = 0.25
gamma = 0.5
[../]
[./poromechskeletoncoupling_x]
type = PoroMechanicsCoupling
variable = u_x
porepressure = p
component = 0
[../]
[./poromechskeletoncoupling_y]
type = PoroMechanicsCoupling
variable = u_y
porepressure = p
component = 1
[../]
[./poromechfluidcoupling_x]
type = PoroMechanicsCoupling
variable = w_x
porepressure = p
component = 0
[../]
[./poromechfluidcoupling_y]
type = PoroMechanicsCoupling
variable = w_y
porepressure = p
component = 1
[../]
[./massconservationskeleton]
type = MassConservationNewmark
variable = p
displacements = 'u_x u_y'
velocities = 'v_x v_y'
accelerations = 'a_x a_y'
beta = 0.25
gamma = 0.5
[../]
[./massconservationfluid]
type = INSMass
variable = p
u = w_x
v = w_y
pressure = p
[../]
[]
[AuxKernels]
[./accel_x]
type = NewmarkAccelAux
variable = a_x
displacement = u_x
velocity = v_x
beta = 0.25
execute_on = timestep_end
[../]
[./accel_y]
type = NewmarkAccelAux
variable = a_y
displacement = u_y
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
darcyvel = w_x
gamma = 0.5
execute_on = timestep_end
[../]
[./fluidaccel_y]
type = NewmarkPoreFluidAccelAux
variable = af_y
darcyvel = w_y
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
displacements = 'u_x u_y'
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
prop_values = 0.0001
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
type = PresetDisplacement
variable = u_y
acceleration = a_y
velocity = v_y
boundary = 'bottom'
#function = bc_func2
beta = 0.25
[../]
[./topright_y]
type = FunctionNeumannBC
  variable = u_y
  boundary = 'topright'
  function = bc_func
[../]
[./topleft_y]
type =  FunctionNeumannBC
variable = u_y
boundary = 'topleft'
function = bc_func2
[../]
[./left_x]
type = PresetDisplacement
variable = u_x
acceleration = a_x
velocity = v_x
boundary = 'left'
function = bc_func2
beta = 0.25
[../]
[./right_x]
type = PresetDisplacement
variable = u_x
acceleration = a_x
velocity = v_x
boundary = 'right'
function = bc_func2
beta = 0.25
[../]
[./fluidbottom_y]
type = FunctionDirichletBC 
variable = w_y
boundary = 'bottom'
function = bc_func2
[../]
[./fluidleft_x]
type = FunctionDirichletBC 
variable = w_x
boundary = 'left'
function = bc_func2
[../]
[./fluidright_x]
type = FunctionDirichletBC 
variable = w_x
boundary = 'right'
function = bc_func2
[../]
[./fluidtopright_y]
type = FunctionDirichletBC 
variable = w_y
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
variable = u_y
[../]
[./leftcornerdisp]
type = PointValue
point = '0 10 0'
variable = u_y
[../]
[./rightcornerdarcy]
type = PointValue
point = '10 10 0'
variable = w_y
[../]
[./leftcornerdarcy]
type = PointValue
point = '0 10 0'
variable = w_y
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
[]
[Outputs]
exodus = true
[]