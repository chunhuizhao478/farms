/*
Material Description of Slip Weakening Friction
*/

#include "PoroSlipWeakeningFriction2dNoInertia.h"
#include "InterfaceKernel.h"

//User-Specific App name
registerMooseObject("farmsApp", PoroSlipWeakeningFriction2dNoInertia);

InputParameters
PoroSlipWeakeningFriction2dNoInertia::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("pressure_plus","Pressure at postive side of the fault");
  params.addRequiredCoupledVar("pressure_minus","Pressure at minus side of the fault");
  return params;
}

PoroSlipWeakeningFriction2dNoInertia::PoroSlipWeakeningFriction2dNoInertia(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _nodal_area(coupledValue("nodal_area")),
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_pressure_plus(coupledNeighborValue("pressure_plus")),
    _interface_pressure_minus(coupledValue("pressure_minus")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density"))
     
{
}

//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
double PoroNoInertiacomputeMusDistribution2D(Real x_coord)
{
  double mu_s = 0.0;
  if (x_coord>=-15.0e3 && x_coord<=15.0e3)
  {
    mu_s = 0.677;
  }
  else
  {
    mu_s = 10000.0;
  }
  return mu_s;
} 

//Compute Spatial Distribution of Background Shear Stress (T1_o)
//Problem-Specific: TPV205-2D
double PoroNoInertiacomputeT1oDistribution2D(Real x_coord)
{
  double T1_o = 0.0;
  if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3)))
  {
      T1_o = 81.6e6;
  }
  else
  {
      T1_o = 70.0e6;
  }
  return T1_o;
}

void
PoroSlipWeakeningFriction2dNoInertia::computeInterfaceTractionAndDerivatives()
{   
  //Incremental displacement jump 
  RealVectorValue displacement_jump_rate = (_interface_displacement_jump[_qp] - _interface_displacement_jump_older[_qp])*(1/_dt);

  //Stress state 
  //RealVectorValue stress_x = _stress[_qp].row(0); 
  RealVectorValue stress_y = _stress[_qp].row(1);  
  RealVectorValue dstress_y = _dstress[_qp].row(1);  
  //RealVectorValue stress_z = _stress[_qp].row(2); 
  Real shear = stress_y(0); 
  Real dshear = dstress_y(0); 
  //Real normal = stress_y(1);


  //Parameter initialization
  Real mu_s = 0; 
  Real mu_d = _mu_d; 
  Real Dc = _Dc; 
  Real tau_f =0;
  Real dtau_f =0;
  Real T1_o = 0;
  Real dT1d_disp =0;
  Real dfd_disp =0;
  Real T2_o = _T2_o;
  Real area = _nodal_area[_qp];


  //pressure state 
   // Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
  //  Real p_plus = _interface_pressure_plus[_qp];
  //  Real p_minus = _interface_pressure_minus[_qp];
   Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
   // Real p = p_minus;
   Real p_plus = _interface_pressure_plus[_qp];
   Real p_minus = _interface_pressure_minus[_qp];
  //  Real p = 0;

  //Obtain x_coord of current qp
  Real x_coord =_q_point[_qp](0);

  //Compute node mass
  Real M = _density[_qp] * area * area/2; 

  //Compute mu_s for current qp
  mu_s = PoroNoInertiacomputeMusDistribution2D(x_coord);

  //Compute T1_o for current qp
  T1_o = PoroNoInertiacomputeT1oDistribution2D(x_coord);

  //Compute sticking stress
  Real T1  = -(1/_dt)*M*displacement_jump_rate(1)/(2*area) + shear + T1_o;
  Real T2  = - T2_o + p;

//Compute friction strength
  if (std::abs(_interface_displacement_jump[_qp](1)) < Dc)
  {
    tau_f = (mu_s - (mu_s - mu_d)*std::abs(_interface_displacement_jump[_qp](1))/Dc)*(-T2); // square for shear component
    dtau_f = -(mu_s - mu_d)/Dc*(-T2)*(_interface_displacement_jump[_qp](1))/(std::abs(_interface_displacement_jump[_qp](1)));
  }
  else
  {
    tau_f = mu_d * (-T2);
    dtau_f = 0;
  }

  //Compute fault traction
  if (std::abs(T1)<tau_f)
  {
    dT1d_disp =  -(1/_dt)*( M )* (1/_dt) / (2*area) + dshear ; 
  }else{
     T1 = tau_f*T1/std::abs(T1);
     dfd_disp =  -(1/_dt)*( M )* (1/_dt) / (2*area) +  dshear ; 
     dT1d_disp =  dtau_f *T1/std::abs(T1) + tau_f  * (dfd_disp*std::abs(T1)-T1*std::copysign(1.0, T1)*dfd_disp)/std::pow(std::abs(T1), 2);  
  }

  //Assign back traction as BC in CZM
  RealVectorValue traction;

  traction(0) = 0; 
  traction(1) = -(T1-T1_o);
  traction(2) = 0;

  //Assign back dtraction as BC in CZM
  RankTwoTensor dtraction;

  dtraction(0,0) = 0; 
  dtraction(0,1) = 0;
  dtraction(0,2) = 0; 

  dtraction(1,0) = 0; 
  dtraction(1,1) = -dT1d_disp;
  dtraction(1,2) = 0;

  dtraction(2,0) = 0; 
  dtraction(2,1) = 0;
  dtraction(2,2) = 0; 

  // std::cout << "T1: " << -T1 << std::endl;
  //  std::cout << "dshear: " <<   dshear << std::endl;
  // std::cout << "shear " << shear << std::endl;
  
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = dtraction;
}