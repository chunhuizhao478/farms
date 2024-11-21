/*
Material Description of Slip Weakening Friction
*/

#include "PoroSlipWeakening2dImpermeable2.h"
#include "InterfaceKernel.h"

//User-Specific App name
registerMooseObject("farmsApp", PoroSlipWeakening2dImpermeable2);

InputParameters
PoroSlipWeakening2dImpermeable2::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addCoupledVar("pressure_plus","Pressure at postive side of the fault");
  params.addCoupledVar("pressure_minus","Pressure at minus side of the fault");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  return params;
}

PoroSlipWeakening2dImpermeable2::PoroSlipWeakening2dImpermeable2(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _interface_pressure_plus(coupledValue("pressure_plus")),
    _interface_pressure_minus(coupledValue("pressure_minus")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _nodal_area(coupledValue("nodal_area")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density"))
{
}

//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
double computeMusDistribution2Dporo2(Real x_coord)
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
double computeT1oDistribution2Dporo2(Real x_coord)
{
  double T1_o = 0.0;
  if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3)))
  {
      T1_o = 81.6e6;
  }
  else if ((x_coord<=(-7.5e3+1.5e3))&&(x_coord>=(-7.5e3-1.5e3)))
  {
      T1_o = 70.0e6;
  }
  else if ((x_coord<=(7.5e3+1.5e3))&&(x_coord>=(7.5e3-1.5e3)))
  {
      T1_o = 70.0e6;
  }
  else
  {
      T1_o = 70.0e6;
  }
  return T1_o;
}

void
PoroSlipWeakening2dImpermeable2::computeInterfaceTractionAndDerivatives()
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


  //pressure state 
  //  Real p = std::max(std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]), 0.0);
  Real p = _interface_pressure_minus[_qp];


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


  //Obtain x_coord of current qp
  Real x_coord =_q_point[_qp](0);

  //Compute node mass
  Real Ms = _density[_qp] * area * area / 2 ; 

  //Compute mu_s for current qp
  mu_s = computeMusDistribution2Dporo2(x_coord);

  //Compute T1_o for current qp
  T1_o = computeT1oDistribution2Dporo2(x_coord);

  //Compute sticking stress
  Real T1  = -(1/_dt)*Ms*displacement_jump_rate(1)/(2*area) + shear + T1_o;
  Real T2  = - (T2_o-p);

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
     dT1d_disp =  -(1/_dt)*( Ms ) * (1/_dt)/ (2*area) + dshear ; 
  }else{
     T1 = tau_f*T1/std::abs(T1);
     dfd_disp =  -(1/_dt)*( Ms ) * (1/_dt)/ (2*area) + dshear ; 
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
  // dtraction(0,3) = 0; 
  // dtraction(0,4) = 0;
  // dtraction(0,5) = 0; 

  dtraction(1,0) = 0; 
  dtraction(1,1) = -dT1d_disp;
  dtraction(1,2) = 0; 
  // dtraction(1,3) = 0; 
  // dtraction(1,4) = dT1d_fluid_vel;
  // dtraction(1,5) = 0; 

   dtraction(2,0) = 0; 
   dtraction(2,1) = 0;
   dtraction(2,2) = 0; 
  // dtraction(2,3) = 0; 
  // dtraction(2,4) = 0;
  // dtraction(2,5) = 0; 

  // std::cout << "dshear: " <<   dshear << std::endl;
  // std::cout << "dT1d_disp: " <<   dT1d_disp << std::endl;

  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] =  -dT1d_disp;
  // _dinterface_traction_djump_vf[_qp] = 0;
  // _dinterface_traction_dpressure[_qp] = 0;

  // _interface_pressure[_qp] = 0;
  // _dinterface_pressure_djump[_qp] = 0;
  // _dinterface_pressure_djump_vf[_qp] = 0;
}