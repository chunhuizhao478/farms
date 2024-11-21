/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#include "PermeableSlipWeakening2dv2.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", PermeableSlipWeakening2dv2);

InputParameters
PermeableSlipWeakening2dv2::validParams()
{
  InputParameters params = PoroCZMComputeLocalTractionBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("react_x","reaction in x dir");
  params.addRequiredCoupledVar("react_neighbor_x","reaction in x dir neighbor");
  params.addRequiredCoupledVar("react_y","reaction in y dir");
  params.addRequiredCoupledVar("react_neighbor_y","reaction in y dir neighbor");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("fluid_disp_x","fluid displacement in x dir");
  params.addRequiredCoupledVar("fluid_disp_y","fluid displacement in y dir");
  params.addCoupledVar("pressure_plus","Pressure at postive side of the fault");
  params.addCoupledVar("pressure_minus","Pressure at minus side of the fault");
  return params;
}

PermeableSlipWeakening2dv2::PermeableSlipWeakening2dv2(const InputParameters & parameters)
  : PoroCZMComputeLocalTractionBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _biot_coefficient(getMaterialPropertyByName<Real>(_base_name + "biot_coefficient")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _nodal_area(coupledValue("nodal_area")),
    _rhof(getMaterialPropertyByName<Real>(_base_name + "rhof")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _interface_pressure_plus(coupledValue("pressure_plus")),
    _interface_pressure_minus(coupledValue("pressure_minus")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_fluid_vel_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_fluid_vel_jump")),
    _fluid_disp_x(coupledValue("fluid_disp_x")),
    _fluid_disp_neighbor_x(coupledNeighborValue("fluid_disp_x")),
    _fluid_disp_y(coupledValue("fluid_disp_y")),
    _fluid_disp_neighbor_y(coupledNeighborValue("fluid_disp_y")),
    _reaction_x(coupledValue("react_x")),
    _reaction_neighbor_x(coupledValue("react_neighbor_x")),
    _reaction_y(coupledValue("react_y")),
    _reaction_neighbor_y(coupledValue("react_neighbor_y"))
    
{
}

//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
double porocomputeMusDistribution2DReactv2(Real x_coord)
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
double porocomputeT1oDistribution2DReactv2(Real x_coord)
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
PermeableSlipWeakening2dv2::computeInterfaceTractionAndDerivatives()
{   
  //Local displacement jump rate
   RealVectorValue displacement_jump_rate = (_interface_displacement_jump[_qp] - _interface_displacement_jump_older[_qp])*(1/_dt);
  
  //Local displacement jump
   RealVectorValue displacement_jump = _interface_displacement_jump[_qp];

  //Global Fluid Displacement Jump 
   RealVectorValue fluid_disp_jump_global = (_fluid_disp_x[_qp]-_fluid_disp_neighbor_x[_qp],_fluid_disp_y[_qp]-_fluid_disp_neighbor_y[_qp]);
 
  //Local Fluid velocity Jump 
   RealVectorValue fluid_vel_jump =  _interface_fluid_vel_jump_old[_qp];

   //Local Fluid Displacement Jump 
   RealVectorValue fluid_disp_jump = _rot[_qp].transpose() * fluid_disp_jump_global;

  //Stress state 
  //RealVectorValue stress_x = _stress[_qp].row(0); 
  RealVectorValue stress_y = _stress[_qp].row(1);  
  RealVectorValue dstress_y = _dstress[_qp].row(1);  
  //RealVectorValue stress_z = _stress[_qp].row(2); 
  Real shear = stress_y(0); 
  Real dshear = dstress_y(0); 
  //Real normal = stress_y(1);


   //pressure state 
   // Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
  //  Real p_plus = _interface_pressure_plus[_qp];
  //  Real p_minus = _interface_pressure_minus[_qp];
    Real alpha = _biot_coefficient[_qp];
  //  Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
   // Real p = p_minus;
   Real p_plus = 0;
   Real p_minus = 0;
   Real p = 0;


 //Parameter initialization
   Real mu_s = 0; 
   Real mu_d = _mu_d; 
   Real Dc = _Dc; 
   Real tau_f =0;
   Real T1_o = 0;
   Real T2_o = _T2_o;

   Real area = _nodal_area[_qp];

   Real x_coord =_q_point[_qp](0);

   //Compute reaction on plus/minus faces in local coordinate
   RealVectorValue elem_react_plus(_reaction_x[_qp],_reaction_y[_qp]);
   RealVectorValue local_elem_react_plus = _rot[_qp].transpose() * elem_react_plus;

   RealVectorValue elem_react_minus(-_reaction_neighbor_x[_qp],-_reaction_neighbor_y[_qp]);
   RealVectorValue local_elem_react_minus = _rot[_qp].transpose() * elem_react_minus;
   
   Real R_plus_local_y  =    local_elem_react_plus(0);
   //Real Rx_plus  =  - local_elem_react_plus(1);
   Real R_minus_local_y =    local_elem_react_minus(0);
   //Real Rx_minus =  - local_elem_react_minus(1);

   //Compute node mass
   Real M = _density[_qp] * area * area / 2; 
   //Compute fluid mass
   Real Mf = _rhof[_qp] * area * area / 2;

   //Compute mu_s for current qp
   mu_s = porocomputeMusDistribution2DReactv2(x_coord);

   //Compute T1_o for current qp
   T1_o = porocomputeT1oDistribution2DReactv2(x_coord);

   //Compute sticking stress
   Real T1 =  -(1/_dt)*M*displacement_jump_rate(1)/(2*area) - (1/_dt)*Mf*fluid_vel_jump(1)/(2*area)
              + shear + T1_o;
  
   Real T2 =  (1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area) 
              + (1/_dt)*Mf*(fluid_vel_jump(0)+(1/_dt)*fluid_disp_jump(0))/(2*area)  
              +( (R_minus_local_y - R_plus_local_y) / ( 2*area) - alpha * (p_plus - p_minus) ) - T2_o + p;
    // Real T2 = - T2_o;

    //Compute sticking stress derivatives
   Real dT1d_ddispx =  -(1/_dt)*M*(1/_dt)/(2*area); // + d(R_plus-R_minus)/ddispx
   Real dT1d_ddispy =  0; // + d(R_plus-R_minus)/ddispy
   Real dT1d_ddispz =  0; // + d(R_plus-R_minus)/ddispz
                     
   Real dT2d_ddispx =  0; // + d(R_plus-R_minus)/ddispx - d(p_plus - p_minus)/ddispx
   Real dT2d_ddispy =  (1/_dt)*M*(2/_dt)/(2*area) ; // d(R_plus-R_minus)/ddispy - d(p_plus - p_minus)/ddispy
   Real dT2d_ddispz =  0; // d(R_plus-R_minus)/ddispz - d(p_plus - p_minus)/ddispz

   Real dT1d_dfluidvelx = - (1/_dt)*Mf/(2*area);
   Real dT1d_dfluidvely =  0; 
   Real dT1d_dfluidvelz =  0; 
                     
   Real dT2d_dfluidvelx =  0; 
   Real dT2d_dfluidvely =  (1/_dt)*Mf*(2)/(2*area); 
   Real dT2d_dfluidvelz =  0; 

   Real dT1d_dp =  0;
   Real dT2d_dp =  1- alpha; 


   //Compute fault traction
   if (T2<0)
   {
   }else{
     T2 = 0;
   }

     Real dtau_f_ddipsx = 0;
     Real dtau_f_ddipsy = 0;
     Real dtau_f_ddipsz = 0;
     Real dtau_f_dfluidvelx  = 0;
     Real dtau_f_dfluidvely  = 0;
     Real dtau_f_dfluidvelz  = 0;
     Real dtau_f_dp  = 0;

   //Compute friction strength
   if (std::abs(displacement_jump(1)) < Dc)
   {
     tau_f = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-T2); // square for shear component

     dtau_f_ddipsx = -(mu_s - mu_d)/Dc*(-T2)*(displacement_jump(1))/(std::abs(displacement_jump(1)));
     dtau_f_ddipsy = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-dT2d_ddispy);
     dtau_f_ddipsz = 0;

     dtau_f_dfluidvelx  = 0;
     dtau_f_dfluidvely  = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-dT2d_dfluidvely);
     dtau_f_dfluidvelz  = 0;

     Real dtau_f_dp  = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(-dT2d_dp);;
   } 
   else
   {
     tau_f = mu_d * (-T2);

     dtau_f_ddipsx = 0;
     dtau_f_ddipsy = mu_d * (-dT2d_ddispy);
     dtau_f_ddipsz = 0;

     dtau_f_dfluidvelx  = 0;
     dtau_f_dfluidvely  = mu_d * (-dT2d_dfluidvely);
     dtau_f_dfluidvelz  = 0;

     dtau_f_dp  = mu_d * (-dT2d_dp);;
   }
   

   //Compute fault traction
   if (std::abs(T1)<tau_f)
   {
   
   }else{
     T1 = tau_f*T1/std::abs(T1);

     dT1d_ddispx = dtau_f_ddipsx*T1/std::abs(T1);
     dT1d_ddispy = dtau_f_ddipsy*T1/std::abs(T1);
     dT1d_ddispz = dtau_f_ddipsz*T1/std::abs(T1);

     dT1d_dfluidvelx = dtau_f_dfluidvelx*T1/std::abs(T1);
     dT1d_dfluidvely = dtau_f_dfluidvely*T1/std::abs(T1);
     dT1d_dfluidvelz = dtau_f_dfluidvelz*T1/std::abs(T1);

     dT1d_dp = dtau_f_dp*T1/std::abs(T1);
 
   }

   //Assign back traction in CZM
   RealVectorValue traction;

   traction(0) =  T2+T2_o;
   traction(1) = -T1+T1_o; 
   traction(2) = 0;

   RankTwoTensor dtraction_ddisp;

   dtraction_ddisp(0,0) = dT2d_ddispy; 
   dtraction_ddisp(0,1) = dT2d_ddispx;
   dtraction_ddisp(0,2) = 0; 

   dtraction_ddisp(1,0) = -dT1d_ddispy; 
   dtraction_ddisp(1,1) = -dT1d_ddispx;
   dtraction_ddisp(1,2) = 0; 

   dtraction_ddisp(2,0) = 0; 
   dtraction_ddisp(2,1) = 0;
   dtraction_ddisp(2,2) = 0; 

   RankTwoTensor dtraction_dfluidvel;

   dtraction_dfluidvel(0,0) = dT2d_dfluidvely; 
   dtraction_dfluidvel(0,1) = dT2d_dfluidvelx;
   dtraction_dfluidvel(0,2) = 0; 

   dtraction_dfluidvel(1,0) = -dT1d_dfluidvely; 
   dtraction_dfluidvel(1,1) = -dT1d_dfluidvelx;
   dtraction_dfluidvel(1,2) = 0; 

   dtraction_dfluidvel(2,0) = 0; 
   dtraction_dfluidvel(2,1) = 0;
   dtraction_dfluidvel(2,2) = 0; 

   RealVectorValue dtraction_dp;
     
   dtraction_dp(0) = dT2d_dp; 
   dtraction_dp(1) = -dT1d_dp; 
   dtraction_dp(2) = 0;
  
   std::cout << "T2 " << T2+T2_o << std::endl;

   _interface_traction[_qp] = traction;
   _dinterface_traction_djump[_qp] = dtraction_ddisp;
   _dinterface_traction_djump_vf[_qp] = dtraction_dfluidvel;
   _dinterface_traction_dpressure[_qp] = dtraction_dp;
  
   _interface_pressure[_qp] = p;
   _dinterface_pressure_djump[_qp] = 0;
   _dinterface_pressure_djump_vf[_qp] = 0;

}