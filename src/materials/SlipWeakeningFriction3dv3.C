/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#include "SlipWeakeningFriction3dv3.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningFriction3dv3);

InputParameters
SlipWeakeningFriction3dv3::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addParam<Real>("T3_o", 1.0, "Background shear traction in z dir");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z","displacement in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z","reaction in z dir");
  return params;
}

SlipWeakeningFriction3dv3::SlipWeakeningFriction3dv3(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _T3_o(getParam<Real>("T3_o")),
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _nodal_area(coupledValue("nodal_area")),
    _nodal_area_neighbor(coupledValue("nodal_area")),
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
    _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z"))
{
}

//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-3D
double computeMusDistribution3Dv3(Real x_coord, Real z_coord)
{
  double mu_s = 0.0;
  if ((x_coord<=(0.0+15e3))&&(x_coord>=(0.0-15e3)) && (z_coord<=(0.0+15e3))&&(z_coord>=(0.0-0.0)))
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
//Problem-Specific: TPV205-3D
double computeT1oDistribution3Dv3(Real x_coord, Real z_coord)
{
  double T1_o = 0.0;
  if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3))&& (z_coord<=(7.5e3+1.5e3))&&(z_coord>=(7.5e3-1.5e3)))
  {
      T1_o = 81.6e6;
  }
  else if ((x_coord<=(-7.5e3+1.5e3))&&(x_coord>=(-7.5e3-1.5e3)) && (z_coord<=(7.5e3+1.5e3))&&(z_coord>=(7.5e3-1.5e3)))
  {
      T1_o = 78.0e6;
  }
  else if ((x_coord<=(7.5e3+1.5e3))&&(x_coord>=(7.5e3-1.5e3)) && (z_coord<=(7.5e3+1.5e3))&&(z_coord>=(7.5e3-1.5e3)))
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
SlipWeakeningFriction3dv3::computeInterfaceTractionAndDerivatives()
{   
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
   Real mu_s = 0; 
   Real mu_d = _mu_d; 
   Real Dc = _Dc; 
   Real tau_f = 0;
   
   Real T1_o = 0;
   Real T2_o = _T2_o;
   Real T3_o = _T3_o;

   Real area = _nodal_area[_qp];
   Real len  = std::sqrt(area);

   Real x_coord =_q_point[_qp](0);
   Real z_coord =_q_point[_qp](2);

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

   //Compute node mass
   Real M = _density[_qp] * len * len * len / 2; 

   //Compute mu_s for current qp
   mu_s = computeMusDistribution3Dv3(x_coord, z_coord);

   //Compute T1_o for current qp
   T1_o = computeT1oDistribution3Dv3(x_coord, z_coord);

   //Compute sticking stress
   Real T1 =  (1/_dt)*M*displacement_jump_rate(1)/(2*len*len) + (R_plus_local_x - R_minus_local_x)/(2*len*len) + T1_o;
   Real T3 =  (1/_dt)*M*displacement_jump_rate(2)/(2*len*len) + (R_plus_local_z - R_minus_local_z)/(2*len*len) + T3_o;
   Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*len*len) + ( (R_minus_local_y - R_plus_local_y) / ( 2*len*len ) ) - T2_o ;

   //Compute fault traction
   if (T2<0)
   {
   }else{
     T2 = 0;
   }

   //Compute friction strength
   if (std::norm(displacement_jump) < Dc)
   {
     tau_f = (mu_s - (mu_s - mu_d)*std::norm(displacement_jump)/Dc)*(-T2); // square for shear component
   } 
   else
   {
     tau_f = mu_d * (-T2);
   }

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