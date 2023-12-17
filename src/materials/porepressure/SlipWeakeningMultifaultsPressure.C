/* Debug */

#include "SlipWeakeningMultifaultsPressure.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningMultifaultsPressure);

InputParameters
SlipWeakeningMultifaultsPressure::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addRequiredParam<Real>("Dc", "Value of characteristic length");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  params.addRequiredCoupledVar("tria_area","area of triangle element along the faults");
  params.addRequiredParam<Real>("effec_sts_coeff", "biot coefficient");
  params.addRequiredCoupledVar("pressure","pressure");
  return params;
}

SlipWeakeningMultifaultsPressure::SlipWeakeningMultifaultsPressure(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _Dc(getParam<Real>("Dc")),
    _nodal_area(coupledValue("nodal_area")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _disp_slipweakening_x(coupledValue("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x(coupledNeighborValue("disp_slipweakening_x")),
    _disp_slipweakening_y(coupledValue("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y(coupledNeighborValue("disp_slipweakening_y")),
    _reaction_slipweakening_x(coupledValue("reaction_slipweakening_x")),
    _reaction_slipweakening_neighbor_x(coupledNeighborValue("reaction_slipweakening_x")),
    _reaction_slipweakening_y(coupledValue("reaction_slipweakening_y")),
    _reaction_slipweakening_neighbor_y(coupledNeighborValue("reaction_slipweakening_y")),
    _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
    _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
    _mu_s(coupledValue("mu_s")),
    _mu_d(coupledValue("mu_d")),
    _sts_init(getMaterialPropertyByName<RankTwoTensor>(_base_name + "static_initial_stress_tensor_slipweakening")),
    _tria_area(coupledValue("tria_area")),
    _tria_area_neighbor(coupledNeighborValue("tria_area")),
    _alpha(getParam<Real>("effec_sts_coeff")),
    _pressure(coupledValue("pressure"))
{
}

void
SlipWeakeningMultifaultsPressure::computeInterfaceTractionAndDerivatives()
{   
   //Global Displacement Jump
   RealVectorValue displacement_jump_global(_disp_slipweakening_x[_qp]-_disp_slipweakening_neighbor_x[_qp],_disp_slipweakening_y[_qp]-_disp_slipweakening_neighbor_y[_qp]);
   
   //Global Displacement Jump Old
   RealVectorValue displacement_jump_old_global(_disp_slipweakening_x_old[_qp]-_disp_slipweakening_neighbor_x_old[_qp],_disp_slipweakening_y_old[_qp]-_disp_slipweakening_neighbor_y_old[_qp]);

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

   Real area = _nodal_area[_qp];

   //Reaction force in local coordinate
   RealVectorValue R_plus_global(-_reaction_slipweakening_x[_qp],-_reaction_slipweakening_y[_qp], 0);
   RealVectorValue R_minus_global(-_reaction_slipweakening_neighbor_x[_qp],-_reaction_slipweakening_neighbor_y[_qp], 0);

   RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
   RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

   Real R_plus_local_x  = R_plus_local(1);
   Real R_plus_local_y  = R_plus_local(0);
   Real R_minus_local_x = R_minus_local(1);
   Real R_minus_local_y = R_minus_local(0);

    //Compute node mass
    Real M = _density[_qp] * sqrt(3) / 4 * area * area / 3;

    //Compute sticking stress
    Real T1 =  (1/_dt)*M*displacement_jump_rate(1)/(2*area) + (R_plus_local_x - R_minus_local_x)/(2*area) + T1_o;
    
    Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area) + ( (R_minus_local_y - R_plus_local_y) / ( 2*area) ) - T2_o ;

   //Compute fault traction
   // add limits -2MPa
    if (T2<0)
    {
    }else{
    T2 = 0;
    }

   //Compute friction strength
   //normal stress is decreased because of fluid injection
   if ( T1 > 0 ){
      if (std::abs(displacement_jump(1)) < Dc)
      {
        tau_f = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-T2-_alpha*_pressure[_qp]); // square for shear component
      } 
      else
      {
        tau_f = mu_d * (-T2-_alpha*_pressure[_qp]);
      }
   }
   else{
      if (std::abs(displacement_jump(1)) < Dc)
      {
        tau_f = (-mu_s + (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-T2-_alpha*_pressure[_qp]); // square for shear component
      } 
      else
      {
        tau_f = -mu_d * (-T2-_alpha*_pressure[_qp]);
      }
   }
   

   //Compute fault traction
   if ( (T1 < 0 && T1 > tau_f) || (T1 > 0 && T1 < tau_f) ) //stuck
   {

   }else{ //slip //pass diff in traction vector // std::abs(T1) may cause problem from neg to pos 
     if (T1 > 0){ 
        T1 =  1 * tau_f*T1/std::abs(T1);
     }
     else{
        T1 = -1 * tau_f*T1/std::abs(T1);
     }
   }

    // Real x_coord = _q_point[_qp](0);
    // Real y_coord = _q_point[_qp](1);
    // Real origin_x = 0;
    // Real origin_y = -100;
    // if (x_coord >= origin_x-30 and x_coord <= origin_x+30 and y_coord >= origin_y-30 and y_coord <= origin_y+30)
    // {
    //   std::cout<<-T2-_alpha*_pressure[_qp]<<" "<<_alpha*_pressure[_qp]<<std::endl;
    // }

   //Assign back traction in CZM
   RealVectorValue traction;

   traction(0) = T2+T2_o; 
   traction(1) = -T1+T1_o; 
   traction(2) = 0;

   _interface_traction[_qp] = traction;
   _dinterface_traction_djump[_qp] = 0;

}