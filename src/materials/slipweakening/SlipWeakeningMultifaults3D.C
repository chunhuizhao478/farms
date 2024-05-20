/* Debug */

#include "SlipWeakeningMultifaults3D.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningMultifaults3D);

InputParameters
SlipWeakeningMultifaults3D::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z","displacement in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z","reaction in z dir");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  params.addRequiredCoupledVar("tria_area","area of triangle element along the faults");
  params.addRequiredCoupledVar("cohesion","cohesion in shear stress");
  params.addRequiredCoupledVar("forced_rupture_time","time of forced rupture");
  return params;
}

SlipWeakeningMultifaults3D::SlipWeakeningMultifaults3D(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _Dc(getParam<Real>("Dc")),
    _nodal_area(coupledValue("nodal_area")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _disp_slipweakening_x(coupledValue("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x(coupledNeighborValue("disp_slipweakening_x")),
    _disp_slipweakening_y(coupledValue("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y(coupledNeighborValue("disp_slipweakening_y")),
    _disp_slipweakening_z(coupledValue("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z(coupledNeighborValue("disp_slipweakening_z")),
    _reaction_slipweakening_x(coupledValue("reaction_slipweakening_x")),
    _reaction_slipweakening_neighbor_x(coupledNeighborValue("reaction_slipweakening_x")),
    _reaction_slipweakening_y(coupledValue("reaction_slipweakening_y")),
    _reaction_slipweakening_neighbor_y(coupledNeighborValue("reaction_slipweakening_y")),
    _reaction_slipweakening_z(coupledValue("reaction_slipweakening_z")),
    _reaction_slipweakening_neighbor_z(coupledNeighborValue("reaction_slipweakening_z")),
    _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
    _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
    _disp_slipweakening_z_old(coupledValueOld("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z")),
    _mu_s(coupledValue("mu_s")),
    _mu_d(coupledValue("mu_d")),
    _sts_init(getMaterialPropertyByName<RankTwoTensor>(_base_name + "static_initial_stress_tensor_slipweakening")),
    _tria_area(coupledValue("tria_area")),
    _tria_area_neighbor(coupledNeighborValue("tria_area")),
    _Co(coupledValue("cohesion")),
    _T(coupledValue("forced_rupture_time"))
{
}

void
SlipWeakeningMultifaults3D::computeInterfaceTractionAndDerivatives()
{   

   //std::cout<<"pass"<<std::endl;

   //Global Displacement Jump
   RealVectorValue displacement_jump_global(_disp_slipweakening_x[_qp]-_disp_slipweakening_neighbor_x[_qp],_disp_slipweakening_y[_qp]-_disp_slipweakening_neighbor_y[_qp],_disp_slipweakening_z[_qp]-_disp_slipweakening_neighbor_z[_qp]);
   
   //Global Displacement Jump Old
   RealVectorValue displacement_jump_old_global(_disp_slipweakening_x_old[_qp]-_disp_slipweakening_neighbor_x_old[_qp],_disp_slipweakening_y_old[_qp]-_disp_slipweakening_neighbor_y_old[_qp],_disp_slipweakening_z_old[_qp]-_disp_slipweakening_neighbor_z_old[_qp]);

   //Global Displacement Jump Rate
   RealVectorValue displacement_jump_rate_global = (displacement_jump_global - displacement_jump_old_global)*(1/_dt);

   //Local Displacement Jump / Displacement Jump Rate
   RealVectorValue displacement_jump      = _rot[_qp].transpose() * displacement_jump_global;
   RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;

   //Parameter initialization
   Real mu_s = _mu_s[_qp]; 
   Real mu_d = _mu_d[_qp]; 
   Real Dc = _Dc; 
   Real tau_f = 0;
   
    //Involve Background Stress Projection
    //Local Init Stress
    RankTwoTensor sts_init_local = _rot[_qp].transpose() * _sts_init[_qp] * _rot[_qp];
    RealVectorValue local_normal(1.0,0.0,0.0);

    //Local Traction
    RealVectorValue traction_local =  sts_init_local * local_normal;

    Real T1_o = -traction_local(1); 
    Real T2_o = -traction_local(0); 
    Real T3_o = -traction_local(2); 

   Real area = _nodal_area[_qp];

   //Reaction force in local coordinate
   RealVectorValue R_plus_global(-_reaction_slipweakening_x[_qp],-_reaction_slipweakening_y[_qp], -_reaction_slipweakening_z[_qp]);
   RealVectorValue R_minus_global(-_reaction_slipweakening_neighbor_x[_qp],-_reaction_slipweakening_neighbor_y[_qp], -_reaction_slipweakening_neighbor_z[_qp]);

   RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
   RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

   Real R_plus_local_x  = R_plus_local(1);
   Real R_plus_local_y  = R_plus_local(0);
   Real R_plus_local_z  = R_plus_local(2);
   Real R_minus_local_x = R_minus_local(1);
   Real R_minus_local_y = R_minus_local(0);
   Real R_minus_local_z = R_minus_local(2);

    //Compute node mass //equal length tetrahedron
    Real M = _density[_qp] * sqrt(3) / 8 * area * area * area / 3;

    //Compute sticking stress
    Real T1 =   (1/_dt)*M*displacement_jump_rate(1)/(2*area*area) + (R_plus_local_x - R_minus_local_x)/(2*area*area) + T1_o;
    Real T3 =   (1/_dt)*M*displacement_jump_rate(2)/(2*area*area) + (R_plus_local_z - R_minus_local_z)/(2*area*area) + T3_o;
    Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area*area) + ( (R_minus_local_y - R_plus_local_y) / ( 2*area*area ) ) - T2_o ;

    //Compute fault traction
   if (T2<0)
   {
   }else{
     T2 = 0;
   }

   //Compute friction strength
  //  if (std::norm(displacement_jump) < Dc)
  //  {
  //    tau_f = (mu_s - (mu_s - mu_d)*std::norm(displacement_jump)/Dc)*(-T2); // square for shear component
  //  } 
  //  else
  //  {
  //    tau_f = mu_d * (-T2);
  //  }

  //parameter f1
  Real f1 = 0.0;
  if ( std::norm(displacement_jump) < Dc ){
    f1 = ( 1.0 * std::norm(displacement_jump) ) / ( 1.0 * Dc );
  }
  else{
    f1 = 1;
  }

  //parameter f2
  Real f2 = 0.0;
  Real t0 = 0.05; //0.5; //s //reduce the nucleation time MODIFIED TPV24
  Real T = _T[_qp];
  if ( _t < T ){
    f2 = 0.0;
  }
  else if ( _t > T && _t < T + t0 ){
    f2 = ( _t - T ) / t0;
  }
  else{
    f2 = 1;
  }

  Real mu = mu_s + ( mu_d - mu_s ) * std::max(f1,f2);

  //Pf
  Real z_coord = _q_point[_qp](2);
  Real fluid_density = 1000; //kg/m^3 fluid density
  Real gravity = 9.8; //m/s^2
  Real Pf = fluid_density * gravity * abs(z_coord);

  //tau_f
  //T2: total normal stress acting on the fault, taken to be positive in compression: abs(T2)
  Real effective_stress = abs(T2) - Pf;
  tau_f = _Co[_qp] + mu * std::max(effective_stress,0.0);

  //Compute fault traction
  if (std::sqrt(T1*T1 + T3*T3)<tau_f)
  {
  
  }else{
    T1 = tau_f*T1/std::sqrt(T1*T1 + T3*T3);
    T3 = tau_f*T3/std::sqrt(T1*T1 + T3*T3);
  }

  //Assign back traction in CZM
  RealVectorValue traction;

  traction(0) = T2+T2_o; 
  traction(1) = -T1+T1_o; 
  traction(2) = -T3+T3_o;

  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;

}