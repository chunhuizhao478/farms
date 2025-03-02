/* 
This file aims to take care of special cases: fault opening/rejoining, reversal in the direction of slip in slip weakening simulation
Created by Chunhui Zhao, Dec 28, 2023
*/

#include "SlipWeakeningMultifaultsGeneralizedPropCSV.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningMultifaultsGeneralizedPropCSV);

InputParameters
SlipWeakeningMultifaultsGeneralizedPropCSV::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBaseSWF2D::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addRequiredCoupledVar("D_c","characteristic length spatial distribution");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("ini_shear_sts","initial shear stress spatial distribution");
  params.addRequiredCoupledVar("ini_normal_sts","initial normal stress spatial distribution");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  return params;
}

SlipWeakeningMultifaultsGeneralizedPropCSV::SlipWeakeningMultifaultsGeneralizedPropCSV(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseSWF2D(parameters),
    _Dc(coupledValue("D_c")),
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
    _ini_shear_sts(coupledValue("ini_shear_sts")),
    _ini_normal_sts(coupledValue("ini_normal_sts")),
    _flag_track_opening_old(getMaterialPropertyOldByName<Real>("flag_track_opening")),
    _flag_track_activecase_old(getMaterialPropertyOldByName<Real>("flag_track_activecase")),
    _jump_track_opening_old(getMaterialPropertyOldByName<Real>("jump_track_opening")),
    _jump_track_reversal_old(getMaterialPropertyOldByName<Real>("jump_track_reversal")),
    _T1_old(getMaterialPropertyOldByName<Real>("T1"))
{
}

void
SlipWeakeningMultifaultsGeneralizedPropCSV::computeInterfaceTractionAndDerivatives()
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
   Real Dc = _Dc[_qp]; 
   Real tau_f = 0;

   //Extract from Function
   Real T1_o = _ini_shear_sts[_qp];
   Real T2_o = -1 * _ini_normal_sts[_qp];

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
    //Real M = _density[_qp] * area * area / 2; 
    Real M = _density[_qp] * sqrt(3) / 4 * area * area / 3;

    //Compute sticking stress
    Real T1 =  (1/_dt)*M*displacement_jump_rate(1)/(2*area) + (R_plus_local_x - R_minus_local_x)/(2*area) + T1_o;
    
    Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area) + ( (R_minus_local_y - R_plus_local_y) / ( 2*area) ) - T2_o ;

   //Compute fault traction
    if (T2<0)
    {
    }else{
    T2 = 0;
    }

   //special cases:
   //initialize
   Real flag_track_opening = _flag_track_opening_old[_qp]; //flag tracks open (0 close, 1 open)
   Real flag_track_activecase = _flag_track_activecase_old[_qp]; //flag tracks which case is on (0 none, 1 opening/rejoining, 2 slip reversal)

   Real jump_track_opening = _jump_track_opening_old[_qp]; //jump at the instant of rejoin
   Real jump_track_reversal = _jump_track_reversal_old[_qp]; //jump at the instant of shear stress reversal

   Real jump_effective = 0.0; //effective jump after considering special cases

   //fault opening
   if ( T2 == 0 ){ //if the fault opens
      //label flag as 1
      flag_track_opening = 1;
   }
   else if ( T2 < 0 && flag_track_opening == 1 ){ //if the fault rejoins
      //label flag as 0
      flag_track_opening = 0;
      //label flag activecase as 1
      flag_track_activecase = 1;
      //save current jump
      jump_track_opening = displacement_jump(1);
   }
   else{}
   
   //slip reversal
   if ( T1 * _T1_old[_qp] < 0 ){ //if the shear stress flip sign
      //label flag activecase as 2
      flag_track_activecase = 2;
      //save current jump
      jump_track_reversal = displacement_jump(1);
   }
   else{}

   //compute effective jump
   if ( flag_track_activecase == 1 ){ //rejoin case happens
      //compute effective jump due to fault open
      jump_effective = displacement_jump(1) - jump_track_opening;
   }
   else if ( flag_track_activecase == 2 ){ //slip reversal case happens
      //compute effective jump due to slip reversal
      jump_effective = displacement_jump(1) - jump_track_reversal;
   }
   else{ //not special cases
      jump_effective = displacement_jump(1);
   }

   //update material props
   _flag_track_opening[_qp] = flag_track_opening;
   _flag_track_activecase[_qp] = flag_track_activecase;
   _jump_track_opening[_qp] = jump_track_opening;
   _jump_track_reversal[_qp] = jump_track_reversal;
   
   //Compute friction strength
   if ( T1 > 0 ){
      if (std::abs(jump_effective) < Dc)
      {
        tau_f = (mu_s - (mu_s - mu_d)*std::abs(jump_effective)/Dc)*(-T2); // square for shear component
      } 
      else
      {
        tau_f = mu_d * (-T2);
      }
   }
   else{
      if (std::abs(jump_effective) < Dc)
      {
        tau_f = (-mu_s + (mu_s - mu_d)*std::abs(jump_effective)/Dc)*(-T2); // square for shear component
      } 
      else
      {
        tau_f = -mu_d * (-T2);
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

   //save T1, T2
   _T1[_qp] = T1;
   _T2[_qp] = T2;

   //save jump_effective
   _jump_effective[_qp] = jump_effective;

   //Assign back traction in CZM
   RealVectorValue traction;

   traction(0) = T2+T2_o; 
   traction(1) = -T1+T1_o; 
   traction(2) = 0;

   _interface_traction[_qp] = traction;
   _dinterface_traction_djump[_qp] = 0;

}