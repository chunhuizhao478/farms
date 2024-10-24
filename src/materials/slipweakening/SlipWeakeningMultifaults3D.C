/* Debug */

#include "SlipWeakeningMultifaults3D.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SlipWeakeningMultifaults3D);

InputParameters
SlipWeakeningMultifaults3D::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBaseLSW3D::validParams();
  params.addClassDescription("Linear Slip Weakening Traction Separation Law.");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z","displacement in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z","reaction in z dir");
  params.addRequiredCoupledVar("reaction_damp_x","reaction damping in x dir");
  params.addRequiredCoupledVar("reaction_damp_y","reaction damping in y dir");
  params.addRequiredCoupledVar("reaction_damp_z","reaction damping in z dir");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  params.addCoupledVar("cohesion","cohesion in shear stress");
  params.addCoupledVar("forced_rupture_time","time of forced rupture");
  return params;
}

SlipWeakeningMultifaults3D::SlipWeakeningMultifaults3D(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseLSW3D(parameters),
    _Dc(getParam<Real>("Dc")),
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
    _reaction_damp_x(coupledValue("reaction_damp_x")),
    _reaction_damp_neighbor_x(coupledNeighborValue("reaction_damp_x")),
    _reaction_damp_y(coupledValue("reaction_damp_y")),
    _reaction_damp_neighbor_y(coupledNeighborValue("reaction_damp_y")),
    _reaction_damp_z(coupledValue("reaction_damp_z")),
    _reaction_damp_neighbor_z(coupledNeighborValue("reaction_damp_z")),
    _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
    _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
    _disp_slipweakening_z_old(coupledValueOld("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z")),
    _mu_s(coupledValue("mu_s")),
    _mu_d(coupledValue("mu_d")),
    _sts_init(getMaterialPropertyByName<RankTwoTensor>(_base_name + "static_initial_stress_tensor_slipweakening")),
    _Co_coupled(isCoupled("cohesion")),
    _T_coupled(isCoupled("forced_rupture_time")),
    _Co(_Co_coupled ? &coupledValue("cohesion") : nullptr),
    _T(_T_coupled ? &coupledValue("forced_rupture_time") : nullptr),
    _accumulated_slip_along_normal_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_normal")),
    _accumulated_slip_along_strike_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_strike")),
    _accumulated_slip_along_dip_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_dip")),
    _slip_along_normal_old(getMaterialPropertyOldByName<Real>("slip_along_normal")),
    _slip_along_strike_old(getMaterialPropertyOldByName<Real>("slip_along_strike")),
    _slip_along_dip_old(getMaterialPropertyOldByName<Real>("slip_along_dip")),
    _current_elem_volume(_assembly.elemVolume()),
    _neighbor_elem_volume(_assembly.neighborVolume()),
    _current_side_volume(_assembly.sideElemVolume()),
    _jump_track_dip_old(getMaterialPropertyOldByName<Real>("jump_track_dip")),
    _T3_old(getMaterialPropertyOldByName<Real>("T3"))
{
}

void
SlipWeakeningMultifaults3D::computeInterfaceTractionAndDerivatives()
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
    Real mu_s = _mu_s[_qp]; 
    Real mu_d = _mu_d[_qp]; 
    Real Dc = _Dc; 
    Real tau_f = 0;
   
    //Involve Background Stress Projection
    //Local Init Stress
    RankTwoTensor _rot_static;
    _rot_static(0,0) = 0; _rot_static(0,1) = 1; _rot_static(0,2) = 0;
    _rot_static(1,0) = 1; _rot_static(1,1) = 0; _rot_static(1,2) = 0;
    _rot_static(2,0) = 0; _rot_static(2,1) = 0; _rot_static(2,2) = -1;

    RankTwoTensor sts_init_local = _rot_static.transpose() * _sts_init[_qp] * _rot_static;
    RealVectorValue local_normal(1.0,0.0,0.0);

    std::cout<<_rot[_qp]<<std::endl;

    //Local Traction
    RealVectorValue traction_local =  sts_init_local * local_normal;

    Real T1_o = -traction_local(1); 
    Real T2_o = -traction_local(0); 
    Real T3_o = -traction_local(2); 

    //for benchmarks, we don't rotate stress
    // Real T1_o = -1*_sts_init[_qp](0,1); 
    // Real T2_o = -1*_sts_init[_qp](1,1); 
    // Real T3_o = -1*_sts_init[_qp](2,2); 

    // std::cout<<"T1: "<<T1_o<<" "<<T1_o_old<<std::endl;
    // std::cout<<"T2: "<<T2_o<<" "<<T2_o_old<<std::endl;
    // std::cout<<"T3: "<<T3_o<<" "<<T3_o_old<<std::endl;

    //*Restoration Force*

    //Stress Divergence Components (label as stsdivcomp)

    //--------------------------------------------------------------------------------------------------//

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_stsdivcomp(-_reaction_slipweakening_x[_qp],-_reaction_slipweakening_y[_qp], -_reaction_slipweakening_z[_qp]);
    RealVectorValue R_minus_global_stsdivcomp(-_reaction_slipweakening_neighbor_x[_qp],-_reaction_slipweakening_neighbor_y[_qp], -_reaction_slipweakening_neighbor_z[_qp]);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_stsdivcomp = _rot[_qp].transpose() * R_plus_global_stsdivcomp;
    RealVectorValue R_minus_local_stsdivcomp = _rot[_qp].transpose() * R_minus_global_stsdivcomp;

    ///Get Components
    //current time step
    Real R_plus_local_normal_stsdivcomp  = R_plus_local_stsdivcomp(0);
    Real R_plus_local_strike_stsdivcomp  = R_plus_local_stsdivcomp(1);
    Real R_plus_local_dip_stsdivcomp     = R_plus_local_stsdivcomp(2);

    Real R_minus_local_normal_stsdivcomp = R_minus_local_stsdivcomp(0);
    Real R_minus_local_strike_stsdivcomp = R_minus_local_stsdivcomp(1);
    Real R_minus_local_dip_stsdivcomp    = R_minus_local_stsdivcomp(2);

    //--------------------------------------------------------------------------------------------------//

    //Damping Components Contribution (label as dampingcomp)

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_dampingcomp(-_reaction_damp_x[_qp],-_reaction_damp_y[_qp], -_reaction_damp_z[_qp]);
    RealVectorValue R_minus_global_dampingcomp(-_reaction_damp_neighbor_x[_qp],-_reaction_damp_neighbor_y[_qp], -_reaction_damp_neighbor_z[_qp]);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_dampingcomp = _rot[_qp].transpose() * R_plus_global_dampingcomp;
    RealVectorValue R_minus_local_dampingcomp = _rot[_qp].transpose() * R_minus_global_dampingcomp;

    ///Get Components
    //current time step
    Real R_plus_local_normal_dampingcomp  = R_plus_local_dampingcomp(0);
    Real R_plus_local_strike_dampingcomp  = R_plus_local_dampingcomp(1);
    Real R_plus_local_dip_dampingcomp     = R_plus_local_dampingcomp(2);
    
    Real R_minus_local_normal_dampingcomp = R_minus_local_dampingcomp(0);
    Real R_minus_local_strike_dampingcomp = R_minus_local_dampingcomp(1);
    Real R_minus_local_dip_dampingcomp    = R_minus_local_dampingcomp(2);

    //--------------------------------------------------------------------------------------------------//

    //Reaction force in local coordinate
    Real R_plus_local_x  = R_plus_local_strike_stsdivcomp + R_plus_local_strike_dampingcomp;
    Real R_plus_local_y  = R_plus_local_normal_stsdivcomp + R_plus_local_normal_dampingcomp;
    Real R_plus_local_z  = R_plus_local_dip_stsdivcomp + R_plus_local_dip_dampingcomp;
    Real R_minus_local_x = R_minus_local_strike_stsdivcomp + R_minus_local_strike_dampingcomp;
    Real R_minus_local_y = R_minus_local_normal_stsdivcomp + R_minus_local_normal_dampingcomp;
    Real R_minus_local_z = R_minus_local_dip_stsdivcomp + R_minus_local_dip_dampingcomp;

    //
    // std::cout<<_current_elem_volume<<" "<<_neighbor_elem_volume<<" "<<_current_side_volume<<std::endl;

    //Compute node mass //equal length tetrahedron
    //Real M_plus  = _density[_qp] * _current_elem_volume  / 2;
    //Real M_minus = M_plus;
    //Real area    = _current_side_volume;

    //Compute sticking stress
    //Real T1 =   (1/_dt)*M_plus*M_minus*displacement_jump_rate(1)/(area*(M_plus+M_minus)) + (M_minus*R_plus_local_x - M_plus*R_minus_local_x)/(area*(M_plus+M_minus)) + T1_o;
    //Real T3 =   (1/_dt)*M_plus*M_minus*displacement_jump_rate(2)/(area*(M_plus+M_minus)) + (M_minus*R_plus_local_z - M_plus*R_minus_local_z)/(area*(M_plus+M_minus)) + T3_o;
    //Real T2 =  -(1/_dt)*M_plus*M_minus*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(area*(M_plus+M_minus)) + ( (M_plus*R_minus_local_y - M_minus*R_plus_local_y) / (area*(M_plus+M_minus)) ) - T2_o ;

    Real area = 200;

    // //Compute node mass //equal length tetrahedron
    // //Real M = _density[_qp] * sqrt(3) / 8 * area * area * area / 3;
    // Real M = _density[_qp] * area * area * area / ( 6 * sqrt(2) ) * 2;

    // Real area_of_triangle = sqrt(3) / 4 * area * area * 2;

    // //Compute sticking stress
    // Real T1 =   (1/_dt)*M*displacement_jump_rate(1)/(2*area_of_triangle) + (R_plus_local_x - R_minus_local_x)/(2*area_of_triangle) + T1_o;
    // Real T3 =   (1/_dt)*M*displacement_jump_rate(2)/(2*area_of_triangle) + (R_plus_local_z - R_minus_local_z)/(2*area_of_triangle) + T3_o;
    // Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area_of_triangle) + ( (R_minus_local_y - R_plus_local_y) / ( 2*area_of_triangle ) ) - T2_o ;

    //Compute node mass //equal length tetrahedron
    Real M = _density[_qp] * area * area * area / 8;

    //Compute sticking stress
    Real T1 =   (1/_dt)*M*displacement_jump_rate(1)/(2*area*area) + (R_plus_local_x - R_minus_local_x)/(2*area*area) + T1_o;
    Real T3 =   (1/_dt)*M*displacement_jump_rate(2)/(2*area*area) - (R_plus_local_z - R_minus_local_z)/(2*area*area) + T3_o;
    Real T2 =  -(1/_dt)*M*(displacement_jump_rate(0)+(1/_dt)*displacement_jump(0))/(2*area*area) + ( (R_minus_local_y - R_plus_local_y) / ( 2*area*area ) ) - T2_o ;

    // Check reversal of T3
    Real jump_track_dip = _jump_track_dip_old[_qp]; //get old dip jump at traction direction switch
    Real jump_effective_dip = 0.0; //initialize effective jump along dip direction
    Real jump_absolute = 0.0; //initialize absolute jump which will be feed into slip weakening law
    if ( T3 * _T3_old[_qp] < 0 ){ //if there is traction direction switch, update track dip
      jump_track_dip = displacement_jump(2);
    }

    //calculate effective jump for dip direction
    jump_effective_dip = displacement_jump(2) - jump_track_dip;

    //calculate absolute jump
    jump_absolute = std::sqrt(displacement_jump(1)*displacement_jump(1)+jump_effective_dip*jump_effective_dip);

    //save material properties
    _T3[_qp] = T3;
    _jump_track_dip[_qp] = jump_track_dip;

    // //Note: The distance that the node has slipped is path-integrated. For example, if the node slips 0.4 m in one
    // //direction and then 0.1 m in the opposite direction, the value of is 0.5 m (and not 0.3 m).
    // //get current slip components
    // //Real slip_along_normal = displacement_jump(0);
    // Real slip_along_strike = displacement_jump(1);
    // Real slip_along_dip    = displacement_jump(2);
    // //get slip absolute difference value compared with previous step
    // //Real slip_diff_along_normal = abs(slip_along_normal-_slip_along_normal_old[_qp]);
    // Real slip_diff_along_strike = abs(slip_along_strike-_slip_along_strike_old[_qp]);
    // Real slip_diff_along_dip    = abs(slip_along_dip-_slip_along_dip_old[_qp]);
    // //update accumulated slip
    // //Real accumulated_slip_along_normal = _accumulated_slip_along_normal_old[_qp] + slip_diff_along_normal;
    // Real accumulated_slip_along_strike = _accumulated_slip_along_strike_old[_qp] + slip_diff_along_strike;
    // Real accumulated_slip_along_dip = _accumulated_slip_along_dip_old[_qp] + slip_diff_along_dip;
    // //compute total distance using accumulated slip
    // Real total_distance = sqrt(accumulated_slip_along_strike*accumulated_slip_along_strike+accumulated_slip_along_dip*accumulated_slip_along_dip);
    // //update accumulated slip components
    // //_accumulated_slip_along_normal[_qp] = accumulated_slip_along_normal;
    // _accumulated_slip_along_strike[_qp] = accumulated_slip_along_strike;
    // _accumulated_slip_along_dip[_qp] = accumulated_slip_along_dip; 
    // //update slip conponents
    // //_slip_along_normal[_qp] = slip_along_normal;
    // _slip_along_strike[_qp] = slip_along_strike;
    // _slip_along_dip[_qp] = slip_along_dip;

    //region overstress nuleation, same as tpv205, tpv14
    if ( !_T_coupled ){
      
      //Compute fault traction
      //min(0,sigma_N)
      if (T2<0)
      {
      }else{
        T2 = 0;
      }

      //Compute friction strength
      //In Seisol documentation states: tau = -C - min(0,sigma_N) * (mu_s - (mu_s - mu_d)/dc * min(total distance,dc))
      //tau_f = (mu_s - (mu_s - mu_d)*std::min(total_distance,Dc)/Dc)*(-T2);

      //Compute friction strength
      // if (std::norm(displacement_jump) < Dc)
      // {
      //   tau_f = (mu_s - (mu_s - mu_d)*std::norm(displacement_jump)/Dc)*(-T2); // square for shear component
      // }
      // else
      // {
      //   tau_f = mu_d * (-T2);
      // }

      if (jump_absolute < Dc)
      {
        tau_f = (mu_s - (mu_s - mu_d)*jump_absolute/Dc)*(-T2); // square for shear component
      }
      else
      {
        tau_f = mu_d * (-T2);
      }

    }
    //region with forced rupture time nucleation, same as tpv24
    // else{
      
    //   Real f1 = 0.0;
    //   if ( total_distance < Dc ){
    //     f1 = ( 1.0 * total_distance ) / ( 1.0 * Dc );
    //   }
    //   else{
    //     f1 = 1;
    //   }

    //   //parameter f2
    //   Real f2 = 0.0;
    //   Real t0 = 0.5; //0.5;
    //   Real T = (*_T)[_qp];
    //   if ( _t < T ){
    //     f2 = 0.0;
    //   }
    //   else if ( _t > T && _t < T + t0 ){
    //     f2 = ( _t - T ) / t0;
    //   }
    //   else{
    //     f2 = 1;
    //   }

    //   //if mu_s > 0.18, must be boundary location, set high mus without degradation
    //   Real mu = 0.0;
    //   if ( mu_s > 10 ){ //must be boundary elements, set high mus values
    //     mu = mu_s * 1;
    //   }
    //   else{
    //     mu = mu_s + ( mu_d - mu_s ) * std::max(f1,f2);
    //   }

    //   //Pf
    //   Real z_coord = _q_point[_qp](2);
    //   Real fluid_density = 1000; //kg/m^3 fluid density
    //   Real gravity = 9.8; //m/s^2
    //   Real Pf = fluid_density * gravity * abs(z_coord);

    //   //tau_f
    //   //T2: total normal stress acting on the fault, taken to be "positive" in compression: -T2
    //   //treat tension on the fault the same as if the effective normal stress equals zero.
    //   Real effective_stress = (-T2) - Pf;
    //   tau_f = (*_Co)[_qp] + mu * std::max(effective_stress,0.0);

    // }

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