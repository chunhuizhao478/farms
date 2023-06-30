/* Debug */

#include "SlipWeakeningTurkeyBranchv4.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningTurkeyBranchv4);

InputParameters
SlipWeakeningTurkeyBranchv4::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addRequiredCoupledVar("D_c","characteristic length spatial distribution");
  params.addParam<Real>("area",200,"fault area");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("ini_shear_sts","initial shear stress spatial distribution");
  params.addRequiredCoupledVar("ini_normal_sts","initial normal stress spatial distribution");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  params.addRequiredCoupledVar("tria_area","area of triangle element along the faults");
  return params;
}

SlipWeakeningTurkeyBranchv4::SlipWeakeningTurkeyBranchv4(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _Dc(coupledValue("D_c")),
    _area(getParam<Real>("area")),
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
    // _sts_init(getMaterialPropertyByName<RankTwoTensor>(_base_name + "static_initial_stress_tensor")),
    _tria_area(coupledValue("tria_area")),
    _tria_area_neighbor(coupledNeighborValue("tria_area")),
    _ini_shear_sts(coupledValue("ini_shear_sts")),
    _ini_normal_sts(coupledValue("ini_normal_sts"))
{
}

void
SlipWeakeningTurkeyBranchv4::computeInterfaceTractionAndDerivatives()
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
   
    // //Involve Background Stress Projection
    // //Local Init Stress
    // RankTwoTensor sts_init_local = _rot[_qp].transpose() * _sts_init[_qp] * _rot[_qp];
    // RealVectorValue local_normal(1.0,0.0,0.0);

    // //Local Traction
    // RealVectorValue traction_local =  sts_init_local * local_normal;

    // Real T1_o = -traction_local(1); 
    // Real T2_o = -traction_local(0); 

   //Predefined Local Background Stress
   //Real T1_o = _sts_init[_qp](0,1);
   //Real T2_o = -_sts_init[_qp](1,1);

   //Extract from Function
   Real T1_o = _ini_shear_sts[_qp];
   Real T2_o = -1 * _ini_normal_sts[_qp];

  //  Real area = _nodal_area[_qp];
  //  Real area = 50; 
   //Real area = 12.5;
   Real area = _area;

   Real x_coord = _q_point[_qp](0);
   Real y_coord = _q_point[_qp](1);

   //Reaction force in local coordinate
   RealVectorValue R_plus_global(-_reaction_slipweakening_x[_qp],-_reaction_slipweakening_y[_qp], 0);
   RealVectorValue R_minus_global(-_reaction_slipweakening_neighbor_x[_qp],-_reaction_slipweakening_neighbor_y[_qp], 0);

   RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
   RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

   Real R_plus_local_x  = R_plus_local(1);
   Real R_plus_local_y  = R_plus_local(0);
   Real R_minus_local_x = R_minus_local(1);
   Real R_minus_local_y = R_minus_local(0);

    //  Explicit Compute node mass
    //  Real M = _density[_qp] * area * area / 2; 
    //  Real M = _density[_qp] * sqrt(3) / 4 * area * area / 3;
    //  Real M_plus = _density[_qp] * _tria_area[_qp];
    //  Real M_minus = _density[_qp] * _tria_area_neighbor[_qp];

    //  //Compute sticking stress
    //  Real T1 =  (1/_dt)*M_plus*M_minus*displacement_jump_rate(1)/(area*(M_plus + M_minus)) + (M_minus * R_plus_local_x - M_plus * R_minus_local_x)/(area*(M_plus + M_minus)) + T1_o;
    
    //  Real T2 =  -(1/_dt)*M_plus*M_minus*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(area*(M_plus + M_minus)) + ( (M_plus*R_minus_local_y - M_minus*R_plus_local_y) / (area*(M_plus + M_minus)) ) - T2_o ;

    //Compute node mass
    //Real M = _density[_qp] * area * area / 2; 
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
   if ( T1 > 0 ){
      if (std::abs(displacement_jump(1)) < Dc)
      {
        tau_f = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-T2); // square for shear component
      } 
      else
      {
        tau_f = mu_d * (-T2);
      }
   }
   else{
      if (std::abs(displacement_jump(1)) < Dc)
      {
        tau_f = (-mu_s + (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-T2); // square for shear component
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

   //output - main fault
  //  if ( (x_coord > 4591) && (x_coord < 4656) && (y_coord > -560) && (y_coord < -508) ){
   if ( (x_coord > 400) && (x_coord < 1000) && (y_coord > -1) && (y_coord < 1) ){
        //std::cout<<"main fault point"<<std::endl;
        std::cout<<"xcoord: "<<x_coord<<" "<<"ycoord: "<<y_coord<<std::endl;
        std::cout<<"mu_s: "<<mu_s<<std::endl;
        std::cout<<"mu_d: "<<mu_d<<std::endl;
        std::cout<<"Dc: "<<Dc<<std::endl;
        std::cout<<"ini_shear_stress: "<<T1_o<<std::endl;
        std::cout<<"ini_normal_stress: "<<T2_o<<std::endl;
   }

  //  //output - branch
  //  if ( (x_coord > 10000-25) && (x_coord < 10000+25) && (y_coord < 0) ){
  //       std::cout<<"branch point"<<std::endl;
  //       std::cout<<"xcoord: "<<x_coord<<" "<<"ycoord: "<<y_coord<<std::endl;
  //       std::cout<<"T1_o: "<<T1_o<<" "<<"T2_o: "<<T2_o<<std::endl;
  //       std::cout<<"T1: "<<T1<<" "<<"T2: "<<T2<<std::endl;
  //       std::cout<<"_rot: "<<_rot[_qp]<<std::endl;
  //       std::cout<<"_normals: "<<_normals[_qp]<<std::endl;
  //       //std::cout<<"sts_init_local: "<<sts_init_local<<std::endl;
  //  }

   //Assign back traction in CZM
   RealVectorValue traction;

   traction(0) = T2+T2_o; 
   traction(1) = -T1+T1_o; 
   traction(2) = 0;

   _interface_traction[_qp] = traction;
   _dinterface_traction_djump[_qp] = 0;

}