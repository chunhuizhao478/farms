/*
Material Description of Slip Weakening Friction
*/

#include "SlipWeakeningFriction2d.h"
#include "InterfaceKernel.h"

//User-Specific App name
registerMooseObject("farmsApp", SlipWeakeningFriction2d);

InputParameters
SlipWeakeningFriction2d::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addParam<Real>("area", 1.0, "Value of fault edge mesh size");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  return params;
}

SlipWeakeningFriction2d::SlipWeakeningFriction2d(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _area(getParam<Real>("area")),
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),

    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density"))
{
}

//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
double computeMusDistribution2D(Real x_coord)
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
double computeT1oDistribution2D(Real x_coord)
{
  double T1_o = 0.0;
  if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3)))
  {
      T1_o = 81.6e6;
  }
  else if ((x_coord<=(-7.5e3+1.5e3))&&(x_coord>=(-7.5e3-1.5e3)))
  {
      T1_o = 78.0e6;
  }
  else if ((x_coord<=(7.5e3+1.5e3))&&(x_coord>=(7.5e3-1.5e3)))
  {
      T1_o = 62.0e6;
  }
  else
  {
      T1_o = 70.0e6;
  }
  return T1_o;
}

void
SlipWeakeningFriction2d::computeInterfaceTractionAndDerivatives()
{   
  //Incremental displacement jump 
  RealVectorValue displacement_jump_rate = (_interface_displacement_jump[_qp] - _interface_displacement_jump_older[_qp])*(1/_dt);

  //Stress state 
  //RealVectorValue stress_x = _stress[_qp].row(0); 
  RealVectorValue stress_y = _stress[_qp].row(1);  
  //RealVectorValue stress_z = _stress[_qp].row(2); 
  Real shear = stress_y(0); 
  //Real normal = stress_y(1);
  
  //Parameter initialization
  Real mu_s = 0; 
  Real mu_d = _mu_d; 
  Real Dc = _Dc; 
  Real tau_f =0;

  Real T1_o = 0;
  Real T2_o = _T2_o;
  Real area = _area;

  //Obtain x_coord of current qp
  Real x_coord =_q_point[_qp](0);

  //Compute node mass
  Real M = _density[_qp] * area * area / 2 ; 

  //Compute mu_s for current qp
  mu_s = computeMusDistribution2D(x_coord);

  //Compute T1_o for current qp
  T1_o = computeT1oDistribution2D(x_coord);

  //Compute sticking stress
  Real T1  = -(1/_dt)*M*displacement_jump_rate(1)/(2*area) + shear + T1_o;
  Real T2  = - T2_o;

  //Compute friction strength
  if (std::abs(_interface_displacement_jump[_qp](1)) < Dc)
  {
    tau_f = (mu_s - (mu_s - mu_d)*std::abs(_interface_displacement_jump[_qp](1))/Dc)*(-T2); // square for shear component
  }
  else
  {
    tau_f = mu_d * (-T2);
  }

  //Compute fault traction
  if (std::abs(T1)<tau_f)
  {
   
  }else{
     T1 = tau_f*T1/std::abs(T1);
  }

  //Assign back traction as BC in CZM
  RealVectorValue traction;

  traction(0) = 0; 
  traction(1) = -(T1-T1_o);
  traction(2) = 0;

  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;
}