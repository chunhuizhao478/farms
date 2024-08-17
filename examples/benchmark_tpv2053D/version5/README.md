Use *FarmsSlipWeakeningBase::FarmsSlipWeakeningCZM* compute traction, feed the values into interface kernel *FarmsCZM*

The coordinate system follows (Day,2005), for planar fault system in global coordinate:

x, y defines the fault surface, z is the normal to the fault surface

global coordinate is the same as local coordinate after transformation

x - strike
y - dip
z - normal

tpv2053d_test_hex_uniform.i : uniform hex elements
tpv2053d_test_tet.i : adaptive tet elements