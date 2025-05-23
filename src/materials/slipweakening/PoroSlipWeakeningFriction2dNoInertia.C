/*
Material Description of Slip Weakening Friction
*/
#include "FEProblem.h"
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
  params.addParam<Real>("elem_size", 1.0, "Value of element size");
  params.addRequiredCoupledVar("pressure_plus","Pressure at postive side of the fault");
  params.addRequiredCoupledVar("pressure_minus","Pressure at minus side of the fault");
  params.addRequiredCoupledVar("react_x","reaction in x dir");
  params.addRequiredCoupledVar("react_y","reaction in y dir");
  params.addRequiredCoupledVar("react_damp_x","reaction in x dir");
  params.addRequiredCoupledVar("react_damp_y","reaction in y dir");
  params.addRequiredCoupledVar("react_pressure_x","reaction in x dir");
  params.addRequiredCoupledVar("react_pressure_y","reaction in y dir");
  return params;
}

PoroSlipWeakeningFriction2dNoInertia::PoroSlipWeakeningFriction2dNoInertia(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _nodal_area(coupledValue("nodal_area")),
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _nodal_area2(getParam<Real>("elem_size")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_pressure_plus(coupledNeighborValue("pressure_plus")),
    _interface_pressure_minus(coupledValue("pressure_minus")),
    _reaction_x(coupledValue("react_x")),
    _reaction_neighbor_x(coupledNeighborValue("react_x")),
    _reaction_y(coupledValue("react_y")),
    _reaction_neighbor_y(coupledNeighborValue("react_y")),
    _reaction_damp_x(coupledValue("react_damp_x")),
    _reaction_neighbor_damp_x(coupledNeighborValue("react_damp_x")),
    _reaction_damp_y(coupledValue("react_damp_y")),
    _reaction_neighbor_damp_y(coupledNeighborValue("react_damp_y")),
    _reaction_pressure_x(coupledValue("react_pressure_x")),
    _reaction_neighbor_pressure_x(coupledNeighborValue("react_pressure_x")),
    _reaction_pressure_y(coupledValue("react_pressure_y")),
    _reaction_neighbor_pressure_y(coupledNeighborValue("react_pressure_y")),
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
      T1_o = 73.0e6;
  }
  return T1_o;
}

void
PoroSlipWeakeningFriction2dNoInertia::computeInterfaceTractionAndDerivatives()
{   
  if(_fe_problem.getCurrentExecuteOnFlag()=="LINEAR")
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
  Real area2 = _nodal_area2/2;

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

  //Compute reaction on plus/minus faces in local coordinate
   RealVectorValue elem_react_plus(_reaction_x[_qp],_reaction_y[_qp]);
   RealVectorValue local_elem_react_plus = _rot[_qp].transpose() * elem_react_plus;

   RealVectorValue elem_react_minus(-_reaction_neighbor_x[_qp],-_reaction_neighbor_y[_qp]);
   RealVectorValue local_elem_react_minus = _rot[_qp].transpose() * elem_react_minus;
   
   Real R_plus_local_y  =    local_elem_react_plus(0);
   Real Rx_plus  =  - local_elem_react_plus(1);
   Real R_minus_local_y =    local_elem_react_minus(0);
   Real Rx_minus =  - local_elem_react_minus(1);

   //Compute damping on plus/minus faces in local coordinate
   RealVectorValue elem_react_damp_plus(_reaction_damp_x[_qp],_reaction_damp_y[_qp]);
   RealVectorValue local_elem_react_damp_plus = _rot[_qp].transpose() * elem_react_damp_plus;

   RealVectorValue elem_react_damp_minus(-_reaction_neighbor_damp_x[_qp],-_reaction_neighbor_damp_y[_qp]);
   RealVectorValue local_elem_react_damp_minus = _rot[_qp].transpose() * elem_react_damp_minus;
   
   Real R_plus_damp_local_y  =    local_elem_react_damp_plus(0);
   Real Rx_plus_damp  =  - local_elem_react_damp_plus(1);
   Real R_minus_damp_local_y =    local_elem_react_damp_minus(0);
   Real Rx_minus_damp =  - local_elem_react_damp_minus(1);

  //Compute pressure on plus/minus faces in local coordinate
   RealVectorValue elem_react_pressure_plus(_reaction_pressure_x[_qp],_reaction_pressure_y[_qp]);
   RealVectorValue local_elem_react_pressure_plus = _rot[_qp].transpose() * elem_react_pressure_plus;

   RealVectorValue elem_react_pressure_minus(-_reaction_neighbor_pressure_x[_qp],-_reaction_neighbor_pressure_y[_qp]);
   RealVectorValue local_elem_react_pressure_minus = _rot[_qp].transpose() * elem_react_pressure_minus;
   
   Real R_plus_pressure_local_y  =    local_elem_react_pressure_plus(0);
   Real Rx_plus_pressure  =  - local_elem_react_pressure_plus(1);
   Real R_minus_pressure_local_y =    local_elem_react_pressure_minus(0);
   Real Rx_minus_pressure =  - local_elem_react_pressure_minus(1);

  //Compute node mass
  Real M = _density[_qp] * area * area2; 

  //Compute mu_s for current qp
  mu_s = PoroNoInertiacomputeMusDistribution2D(x_coord);

  //Compute T1_o for current qp
  T1_o = PoroNoInertiacomputeT1oDistribution2D(x_coord);

  //Compute sticking stress
    Real T1 =  -(1/_dt)*M*displacement_jump_rate(1)/(2*area)
               +  ( ( Rx_plus + Rx_plus_damp - Rx_minus  - Rx_minus_damp ) / ( 2*area) )  + T1_o;
  
  
    Real T2 =  -T2_o 
          + ( (R_plus_local_y + R_plus_damp_local_y + R_plus_pressure_local_y - R_minus_local_y - R_minus_damp_local_y - R_minus_pressure_local_y) / ( 2*area) );

   Real T2_SW =  T2_o - p;
//Compute friction strength
  if (std::abs(_interface_displacement_jump[_qp](1)) < Dc)
  {
    tau_f = (mu_s - (mu_s - mu_d)*std::abs(_interface_displacement_jump[_qp](1))/Dc)*(T2_SW); // square for shear component
    dtau_f = -(mu_s - mu_d)/Dc*(T2_SW)*(_interface_displacement_jump[_qp](1))/(std::abs(_interface_displacement_jump[_qp](1)));
  }
  else
  {
    tau_f = mu_d * (T2_SW);
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

   //  traction(0) =  -(T2+T2_o); 
  traction(0) = T2+T2_o;
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
}