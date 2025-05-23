/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#include "PoroSlipWeakening2d_sub.h"
#include "InterfaceKernel.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", PoroSlipWeakening2d_sub);

InputParameters
PoroSlipWeakening2d_sub::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("T2_o", 1.0, "Background normal traction");
  params.addParam<Real>("mu_d", 1.0, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addParam<Real>("elem_size", 1.0, "Value of element size");
  params.addRequiredCoupledVar("interface_disp_x","displacement in x dir");
  params.addRequiredCoupledVar("interface_disp_y","displacement in y dir");
  params.addRequiredCoupledVar("react_x","reaction in x dir");
  params.addRequiredCoupledVar("react_y","reaction in y dir");
  params.addRequiredCoupledVar("react_damp_x","reaction in x dir");
  params.addRequiredCoupledVar("react_damp_y","reaction in y dir");
  params.addRequiredCoupledVar("react_pressure_x","reaction in x dir");
  params.addRequiredCoupledVar("react_pressure_y","reaction in y dir");
  params.addRequiredCoupledVar("fluid_disp_x","reaction in x dir");
  params.addRequiredCoupledVar("fluid_disp_y","reaction in y dir");
  params.addRequiredCoupledVar("fluid_vel_x","reaction in x dir");
  params.addRequiredCoupledVar("fluid_vel_y","reaction in y dir");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("pressure_plus","Pressure at postive side of the fault");
  params.addRequiredCoupledVar("pressure_minus","Pressure at minus side of the fault");
  params.addRequiredParam<std::string>("permeability_type", 
    "Type of permeability condition (permeable/impermeable/semi_permeable)");
  return params;
}

PoroSlipWeakening2d_sub::PoroSlipWeakening2d_sub(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _T2_o(getParam<Real>("T2_o")), 
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _nodal_area2(getParam<Real>("elem_size")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _nodal_area(coupledValue("nodal_area")),
    _rhof(getMaterialPropertyByName<Real>(_base_name + "rhof")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _permeability_type(getParam<std::string>("permeability_type")), 
    _disp_x(coupledValue("interface_disp_x")),
    _disp_neighbor_x(coupledNeighborValue("interface_disp_x")),
    _disp_y(coupledValue("interface_disp_y")),
    _disp_neighbor_y(coupledNeighborValue("interface_disp_y")),
    _disp_x_older(coupledValueOld("interface_disp_x")),
    _disp_neighbor_x_older(coupledNeighborValueOld("interface_disp_x")),
    _disp_y_older(coupledValueOld("interface_disp_y")),
    _disp_neighbor_y_older(coupledNeighborValueOld("interface_disp_y")),
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
    _fluid_vel_x(coupledValue("fluid_disp_x")),
    _fluid_vel_neighbor_x(coupledNeighborValue("fluid_disp_x")),
    _fluid_vel_y(coupledValue("fluid_disp_y")),
    _fluid_vel_neighbor_y(coupledNeighborValue("fluid_disp_y")),
    _fluid_disp_x(coupledValue("fluid_disp_x")),
    _fluid_disp_neighbor_x(coupledNeighborValue("fluid_disp_x")),
    _fluid_disp_y(coupledValue("fluid_disp_y")),
    _fluid_disp_neighbor_y(coupledNeighborValue("fluid_disp_y"))
       
{
}


//Compute Spatial Distribution of Static Friction Parameter (mu_s)
//Problem-Specific: TPV205-2D
double porocomputeMusDistribution2DReactv10(Real x_coord)
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
double porocomputeT1oDistribution2DReactv10(Real x_coord)
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
PoroSlipWeakening2d_sub::computeInterfaceTractionAndDerivatives()
{   

  // if(_fe_problem.getCurrentExecuteOnFlag()=="LINEAR")
  //  {
  //Global displacement jump
   RealVectorValue displacement_jump_global = (_disp_x[_qp]-_disp_neighbor_x[_qp],_disp_y[_qp]-_disp_neighbor_y[_qp]);
   RealVectorValue displacement_jump_global_older = (_disp_x_older[_qp]-_disp_neighbor_x_older[_qp],_disp_y_older[_qp]-_disp_neighbor_y_older[_qp]);

  //Global displacement jump rate
   RealVectorValue displacement_jump_rate_global = (displacement_jump_global- displacement_jump_global_older)*(1/_dt);

  //Global Fluid Displacement Jump 
   RealVectorValue fluid_disp_jump_global = (_fluid_disp_x[_qp]-_fluid_disp_neighbor_x[_qp],_fluid_disp_y[_qp]-_fluid_disp_neighbor_y[_qp]);
 
  //Local Fluid velocity Jump 
   RealVectorValue fluid_vel_jump_global = (_fluid_vel_x[_qp]-_fluid_vel_neighbor_x[_qp],_fluid_vel_y[_qp]-_fluid_vel_neighbor_y[_qp]);

   //Local Fluid Displacement Jump 
   RealVectorValue displacement_jump = _rot[_qp].transpose() * displacement_jump_global;
   RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;
   RealVectorValue fluid_vel_jump = _rot[_qp].transpose() * fluid_vel_jump_global;
   RealVectorValue fluid_disp_jump = _rot[_qp].transpose() * fluid_disp_jump_global;

  //Stress state 
  //RealVectorValue stress_x = _stress[_qp].row(0); 
  RealVectorValue stress_y = _stress[_qp].row(1);  
  //RealVectorValue stress_z = _stress[_qp].row(2); 
  Real shear = stress_y(0); 
  Real normal = stress_y(1); 
  
  //pressure state 
  Real p_plus = _interface_pressure_plus[_qp];
  Real p_minus = _interface_pressure_minus[_qp];
  Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
  
  // Define pressure based on fault permeability condition
  if (_permeability_type == "permeable")
  {
     Real p = 0.0;  // For permeable condition
  }

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

    Real area2;
    if (_current_elem->default_order() == SECOND)  // QUAD9
        area2 = _nodal_area2/6.0;
    else  // QUAD4
        area2 = _nodal_area2/2.0;


   //Compute node mass
   Real M = _density[_qp] * area * area2; 
   //Compute fluid mass
   Real Mf = _rhof[_qp] * area * area2; 
   
   //Compute mu_s for current qp
   mu_s = porocomputeMusDistribution2DReactv10(x_coord);

   //Compute T1_o for current qp
   T1_o = porocomputeT1oDistribution2DReactv10(x_coord);

   //Compute sticking stress

    Real T1 =  -(1/_dt)*M*displacement_jump_rate(1)/(2*area) - (1/_dt)*Mf*fluid_vel_jump(1)/(2*area)
              +  ( ( Rx_plus + Rx_plus_damp + Rx_plus_pressure - Rx_minus  - Rx_minus_damp - Rx_minus_pressure) / ( 2*area) )  + T1_o;
  
    Real T2 =   //(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area) 
               0+ (1/_dt)*Mf*(fluid_vel_jump(0)+(1/_dt)*fluid_disp_jump(0))/(2*area)  
               +( (R_plus_local_y + R_plus_damp_local_y + R_plus_pressure_local_y - R_minus_local_y - R_minus_damp_local_y - R_minus_pressure_local_y) / ( 2*area) )
               - T2_o;


   //Compute fault traction
   if (T2<0)
   {
   }else{
     T2 = 0;
   }

   //Compute friction strength
   Real T2_SW = T2_o - p ;

   if (std::abs(displacement_jump(1)) < Dc)
   {
     tau_f = (mu_s - (mu_s - mu_d)*std::abs(displacement_jump(1))/Dc)*(T2_SW); // square for shear component
   } 
   else
   {
     tau_f = mu_d * (T2_SW);
   }
   

   //Compute fault traction
   if (std::abs(T1)<tau_f)
   {
   
   }else{
     T1 = tau_f*T1/std::abs(T1);
 
   }

   //Assign back traction in CZM
   RealVectorValue traction;

   traction(0) =  T2+T2_o;
   traction(1) = -T1+T1_o; 
   traction(2) = 0;

   _interface_traction[_qp] = traction;
   _dinterface_traction_djump[_qp] = 0;

}
// }