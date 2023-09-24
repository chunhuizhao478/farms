//Implementation of Rate-and-State Friction Law in 3D

#include "RateStateFrictionLaw3DAsBC.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", RateStateFrictionLaw3DAsBC);

InputParameters
RateStateFrictionLaw3DAsBC::validParams()
{ 
  //Tn_o,Ts_o,Vini,statevarini are defined in "CZMComputeLocalTractionTotalBaseRSF3D"
  InputParameters params = CZMComputeLocalTractionTotalBaseRSF3D::validParams();
  params.addClassDescription("Rate-and-State Frictional Law.");
  params.addRequiredParam<Real>("len","element length");
  params.addRequiredParam<Real>("f_o","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_a","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_b","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_L","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("delta_o","slip rate parameter");
  params.addRequiredCoupledVar("reaction_rsf_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_rsf_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_rsf_z","reaction in z dir");
  params.addRequiredCoupledVar("reaction_damp_x","reaction damping in x dir");
  params.addRequiredCoupledVar("reaction_damp_y","reaction damping in y dir");
  params.addRequiredCoupledVar("reaction_damp_z","reaction damping in z dir");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir");
  return params;
}

RateStateFrictionLaw3DAsBC::RateStateFrictionLaw3DAsBC(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseRSF3D(parameters),
    _len(getParam<Real>("len")),
    _f_o(getParam<Real>("f_o")),
    _rsf_a(getParam<Real>("rsf_a")),
    _rsf_b(getParam<Real>("rsf_b")),
    _rsf_L(getParam<Real>("rsf_L")),
    _delta_o(getParam<Real>("delta_o")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _reaction_rsf_x(coupledValue("reaction_rsf_x")),
    _reaction_rsf_y(coupledValue("reaction_rsf_y")),
    _reaction_rsf_z(coupledValue("reaction_rsf_z")),
    _reaction_rsf_neighbor_x(coupledNeighborValue("reaction_rsf_x")),
    _reaction_rsf_neighbor_y(coupledNeighborValue("reaction_rsf_y")),
    _reaction_rsf_neighbor_z(coupledNeighborValue("reaction_rsf_z")),
    _reaction_damp_x(coupledValue("reaction_damp_x")),
    _reaction_damp_y(coupledValue("reaction_damp_y")),
    _reaction_damp_z(coupledValue("reaction_damp_z")),
    _reaction_damp_neighbor_x(coupledNeighborValue("reaction_damp_x")),
    _reaction_damp_neighbor_y(coupledNeighborValue("reaction_damp_y")),
    _reaction_damp_neighbor_z(coupledNeighborValue("reaction_damp_z")),
    _Ts_perturb(coupledValue("Ts_perturb")),
    _Ts_perturb_old(coupledValueOld("Ts_perturb")),
    _alongfaultvel_strike_plus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_strike_plus")),
    _alongfaultvel_strike_minus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_strike_minus")),
    _alongfaultvel_normal_plus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_normal_plus")),
    _alongfaultvel_normal_minus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_normal_minus")),
    _alongfaultvel_dip_plus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_dip_plus")),
    _alongfaultvel_dip_minus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_dip_minus")),
    _alongfaultdisp_strike_plus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_strike_plus")),
    _alongfaultdisp_strike_minus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_plus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_normal_minus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_normal_minus")),
    _alongfaultdisp_dip_plus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_dip_plus")),
    _alongfaultdisp_dip_minus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_dip_minus")),
    _alongfaultdisp_strike_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_strike_plus")),
    _alongfaultdisp_normal_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_dip_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_dip_plus")),
    _alongfaultdisp_strike_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_minus")),
    _alongfaultdisp_dip_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_dip_minus")),
    _sliprate_strike_old(getMaterialPropertyOldByName<Real>("sliprate_strike")),
    _sliprate_normal_old(getMaterialPropertyOldByName<Real>("sliprate_normal")),
    _sliprate_dip_old(getMaterialPropertyOldByName<Real>("sliprate_dip")),
    _sliprate_mag_old(getMaterialPropertyOldByName<Real>("sliprate_mag")),
    _sliprate_predict(getMaterialPropertyOldByName<Real>("sliprate_predict")),
    _slip_strike_old(getMaterialPropertyOldByName<Real>("slip_strike")),
    _slip_normal_old(getMaterialPropertyOldByName<Real>("slip_normal")),
    _slip_dip_old(getMaterialPropertyOldByName<Real>("slip_dip")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _statevar_older(getMaterialPropertyOlderByName<Real>("statevar")),
    _traction_strike_old(getMaterialPropertyOldByName<Real>("traction_strike")),
    _traction_normal_old(getMaterialPropertyOldByName<Real>("traction_normal")),
    _traction_dip_old(getMaterialPropertyOldByName<Real>("traction_dip"))
{
}

void
RateStateFrictionLaw3DAsBC::computeInterfaceTractionAndDerivatives()
{   
    //Define Parameters
    Real len = _len;
    Real f_o = _f_o;
    Real rsf_a = _rsf_a;
    Real rsf_b = _rsf_b;
    Real rsf_L = _rsf_L;
    Real delta_o = _delta_o;
    Real Tn_o = _Tn_o;
    Real Ts_o = _Ts_o;
    Real Td_o = _Td_o;
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal
    Real Td = 0.0; //dip
    
    //*Restoration Force*

    //Stress Divergence Components (label as stsdivcomp)

    //--------------------------------------------------------------------------------------------------//

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global(-_reaction_rsf_x[_qp],-_reaction_rsf_y[_qp], -_reaction_rsf_z[_qp]);
    RealVectorValue R_minus_global(-_reaction_rsf_neighbor_x[_qp],-_reaction_rsf_neighbor_y[_qp], -_reaction_rsf_neighbor_z[_qp]);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
    RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

    ///Get Components
    //current time step
    Real R_plus_local_normal_stsdivcomp  = R_plus_local(0);
    Real R_plus_local_strike_stsdivcomp  = R_plus_local(1);
    Real R_plus_local_dip_stsdivcomp     = R_plus_local(2);

    Real R_minus_local_normal_stsdivcomp = R_minus_local(0);
    Real R_minus_local_strike_stsdivcomp = R_minus_local(1);
    Real R_minus_local_dip_stsdivcomp    = R_minus_local(2);

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

    //Add restoration forces from two contributions
    Real R_plus_local_normal  = R_plus_local_normal_stsdivcomp  + R_plus_local_normal_dampingcomp;
    Real R_plus_local_strike  = R_plus_local_strike_stsdivcomp  + R_plus_local_strike_dampingcomp;
    Real R_plus_local_dip     = R_plus_local_dip_stsdivcomp     + R_plus_local_dip_dampingcomp;
    Real R_minus_local_normal = R_minus_local_normal_stsdivcomp + R_minus_local_normal_dampingcomp;
    Real R_minus_local_strike = R_minus_local_strike_stsdivcomp + R_minus_local_strike_dampingcomp;  
    Real R_minus_local_dip    = R_minus_local_dip_stsdivcomp    + R_minus_local_dip_dampingcomp;

    //--------------------------------------------------------------------------------------------------//

    //*Nodal Mass*
    ///HEX8 Element
    Real M = _density[_qp] * len * len * len * 0.5;

    //fault area
    Real len_len = ( 1.0 * len ) * ( 1.0 * len );
    
    //Assign shear perturbation at t
    Real Ts_perturb = _Ts_perturb[_qp];

    //Old Value
    ///Disp Plus Side
    Real alongfaultdisp_strike_plus_t = _alongfaultdisp_strike_plus_old[_qp];
    Real alongfaultdisp_normal_plus_t = _alongfaultdisp_normal_plus_old[_qp];
    Real alongfaultdisp_dip_plus_t    = _alongfaultdisp_dip_plus_old[_qp];
    ///Disp Minus Side
    Real alongfaultdisp_strike_minus_t = _alongfaultdisp_strike_minus_old[_qp];
    Real alongfaultdisp_normal_minus_t = _alongfaultdisp_normal_minus_old[_qp];
    Real alongfaultdisp_dip_minus_t = _alongfaultdisp_dip_minus_old[_qp];

    //Older Value
    //Plus Side
    Real alongfaultdisp_strike_plus_tminust = _alongfaultdisp_strike_plus_older[_qp];
    Real alongfaultdisp_normal_plus_tminust = _alongfaultdisp_normal_plus_older[_qp];
    Real alongfaultdisp_dip_plus_tminust = _alongfaultdisp_dip_plus_older[_qp];
    //Minus Side
    Real alongfaultdisp_strike_minus_tminust = _alongfaultdisp_strike_minus_older[_qp];
    Real alongfaultdisp_normal_minus_tminust = _alongfaultdisp_normal_minus_older[_qp];
    Real alongfaultdisp_dip_minus_tminust = _alongfaultdisp_dip_minus_older[_qp];

    //*Get slip rate and slip from time t-dt/2 and t
    ///strike direction
    Real sliprate_strike_tminusdtover2 = _sliprate_strike_old[_qp];
    Real slip_strike_t = _slip_strike_old[_qp];

    ///dip direction
    Real sliprate_dip_tminusdtover2 = _sliprate_dip_old[_qp];
    Real slip_dip_t = _slip_dip_old[_qp];

    ///normal direction
    Real sliprate_normal_tminusdtover2 = _sliprate_normal_old[_qp];
    Real slip_normal_t = _slip_normal_old[_qp];

    ///mag
    Real sliprate_mag_tminusdtover2 = _sliprate_mag_old[_qp];

    //Get State Variable at Current Time Step
    Real statevar_t = _statevar_old[_qp];

    //*Compute Trial Normal Traction and Normal Traction*
    Real Tn_trial = ( -1.0 * (1.0/_dt) * M * M * ( sliprate_normal_tminusdtover2 + (1.0/_dt) * slip_normal_t) ) / ( len_len * (M + M) ) + ( M * R_minus_local_normal - M * R_plus_local_normal ) / ( len_len * (M + M) ) - Tn_o;
    if (Tn_trial<0)
    {
        Tn = Tn_trial;
    }else{
        Tn = 0;
    }

    ///Make Tn positive
    Tn = abs(Tn);

    //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
    Real Ts_trial = ( M * M * sliprate_strike_tminusdtover2 )/( len_len * _dt * (M + M) ) + (M * R_plus_local_strike - M * R_minus_local_strike) / ( len_len * ( M + M ) ) + Ts_o + Ts_perturb;
    Real Td_trial = ( M * M * sliprate_dip_tminusdtover2    )/( len_len * _dt * (M + M) ) + (M * R_plus_local_dip    - M * R_minus_local_dip   ) / ( len_len * ( M + M ) ) + Td_o;
    Real Tmag_trial = sqrt(Ts_trial*Ts_trial+Td_trial*Td_trial);

    //const
    Real c = len_len * _dt * ( M + M ) / (M * M);
    Real Z = 0.5 / delta_o * exp((f_o + rsf_b * log(delta_o * statevar_t/rsf_L))/rsf_a);    
    
    //Setup while loop
    Real iterr = 1;
    Real max_iter = 10000;
    Real er = 1;
    Real solution; 
    Real guess_i = abs(sliprate_mag_tminusdtover2); //slip rate at time t-dt/2
    Real residual;
    Real jacobian;
    Real guess_j;
    while ( er > 1e-10 && iterr < max_iter ){  
        
        //Compute Residual
        residual = guess_i + c * Tn * rsf_a * asinh( 0.5*(guess_i+sliprate_mag_tminusdtover2) * Z ) - c * Tmag_trial;

        //Compute Jacobian
        jacobian = 1.0 + c * Tn * rsf_a * 0.5 * Z / sqrt( 1.0 + 0.5 * 0.5 * (guess_i+sliprate_mag_tminusdtover2) * (guess_i+sliprate_mag_tminusdtover2) * Z * Z );

        //Compute New guess
        guess_j = guess_i - residual / jacobian;

        //save
        solution = guess_j;

        //Compute err
        er = abs(guess_j - guess_i)/abs(guess_j);

        //Update Old guess
        guess_i = guess_j;

        if (iterr == max_iter){
            std::cout<<"NOT CONVERGED!"<<std::endl;
        }

        //update iterr
        iterr = iterr + 1;
    
    }

    Real sliprate_mag_tplusdtover2 = abs(solution); 
    
    //update state variable
    // theta_ref = ( theta_pre + _dt ) / ( 1.0 + _dt * dv_pre / rsf_L ); //febealg2dot1
    //theta_ref = ( _statevar_older[_qp] + 2 * _dt ) / ( 1.0 + 2 * _dt * dv_pre / rsf_L );
    Real coeffD = exp(-sliprate_mag_tplusdtover2*_dt/rsf_L);
    Real statevar_tplusdt = statevar_t * coeffD + (rsf_L/sliprate_mag_tplusdtover2) * (1-coeffD);

    //*Compute shear traction at time t*
    Real T_mag = Tn * rsf_a * asinh( 0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2) * Z );

    ///Get Components
    Ts = T_mag * ( Ts_trial / Tmag_trial );
    Td = T_mag * ( Td_trial / Tmag_trial );

    //*Compute quantities at time t + dt/2 or t + dt
    ////DISP
    Real du_strike_plus_t  =  alongfaultdisp_strike_plus_t  - alongfaultdisp_strike_plus_tminust  + _dt * _dt / M * (R_plus_local_strike  - len_len * (Ts - Ts_o - Ts_perturb));
    Real du_normal_plus_t  =  alongfaultdisp_normal_plus_t  - alongfaultdisp_normal_plus_tminust  + _dt * _dt / M * (R_plus_local_normal  - len_len * (Tn - Tn_o));
    Real du_dip_plus_t     =  alongfaultdisp_dip_plus_t     - alongfaultdisp_dip_plus_tminust     + _dt * _dt / M * (R_plus_local_dip     - len_len * (Td - Td_o));
    
    Real du_strike_minus_t =  alongfaultdisp_strike_minus_t  - alongfaultdisp_strike_minus_tminust + _dt * _dt / M * (R_minus_local_strike + len_len * (Ts - Ts_o - Ts_perturb));
    Real du_normal_minus_t =  alongfaultdisp_normal_minus_t  - alongfaultdisp_normal_minus_tminust + _dt * _dt / M * (R_minus_local_normal + len_len * (Tn - Tn_o));
    Real du_dip_minus_t    =  alongfaultdisp_dip_minus_t     - alongfaultdisp_dip_minus_tminust    + _dt * _dt / M * (R_minus_local_dip    + len_len * (Td - Td_o));

    Real alongfaultdisp_strike_plus_tplusdt  = alongfaultdisp_strike_plus_t  + du_strike_plus_t;
    Real alongfaultdisp_normal_plus_tplusdt  = alongfaultdisp_normal_plus_t  + du_normal_plus_t;
    Real alongfaultdisp_dip_plus_tplusdt     = alongfaultdisp_dip_plus_t     + du_dip_plus_t;
    
    Real alongfaultdisp_strike_minus_tplusdt = alongfaultdisp_strike_minus_t + du_strike_minus_t;
    Real alongfaultdisp_normal_minus_tplusdt = alongfaultdisp_normal_minus_t + du_normal_minus_t;
    Real alongfaultdisp_dip_minus_tplusdt    = alongfaultdisp_dip_minus_t    + du_dip_minus_t;

    //Update Slip Rate and Slip at t + dt/2 or t + dt
    ///Slip Rate
    Real sliprate_strike_tplusdtover2 = sliprate_mag_tplusdtover2 * ( Ts_trial / Tmag_trial );
    Real sliprate_normal_tplusdtover2 = ( alongfaultdisp_normal_plus_tplusdt - alongfaultdisp_normal_plus_t ) / _dt - ( alongfaultdisp_normal_minus_tplusdt - alongfaultdisp_normal_minus_t ) / _dt;
    Real sliprate_dip_tplusdtover2    = sliprate_mag_tplusdtover2 * ( Td_trial / Tmag_trial );

    ///Slip
    Real slip_strike_tplusdtover2 = slip_strike_t + sliprate_strike_tplusdtover2 * _dt;
    Real slip_normal_tplusdtover2 = alongfaultdisp_normal_plus_tplusdt - alongfaultdisp_normal_minus_tplusdt;
    Real slip_dip_tplusdtover2    = slip_dip_t + sliprate_dip_tplusdtover2 * _dt;

    //Save Results
    ///Traction (t)
    _traction_strike[_qp] = Ts;
    _traction_normal[_qp] = Tn;
    _traction_dip[_qp] = Td;

    ///Disp (t+dt)
    _alongfaultdisp_strike_plus[_qp] = alongfaultdisp_strike_plus_tplusdt;
    _alongfaultdisp_normal_plus[_qp] = alongfaultdisp_normal_plus_tplusdt;
    _alongfaultdisp_dip_plus[_qp]    = alongfaultdisp_dip_plus_tplusdt;

    _alongfaultdisp_strike_minus[_qp] = alongfaultdisp_strike_minus_tplusdt;
    _alongfaultdisp_normal_minus[_qp] = alongfaultdisp_normal_minus_tplusdt;
    _alongfaultdisp_dip_minus[_qp]    = alongfaultdisp_dip_minus_tplusdt;

    ///Slip Rate (t+dt/2)
    _sliprate_strike[_qp] = sliprate_strike_tplusdtover2;
    _sliprate_normal[_qp] = sliprate_normal_tplusdtover2;
    _sliprate_dip[_qp]    = sliprate_dip_tplusdtover2;
    _sliprate_mag[_qp]    = sliprate_mag_tplusdtover2;

    ///Slip (t+dt)
    _slip_strike[_qp] = slip_strike_tplusdtover2;
    _slip_normal[_qp] = slip_normal_tplusdtover2;
    _slip_dip[_qp]    = slip_dip_tplusdtover2;
    _slip_mag[_qp]    = sqrt(slip_strike_tplusdtover2*slip_strike_tplusdtover2+slip_dip_tplusdtover2*slip_dip_tplusdtover2);

    ///State Variable (t+dt)
    _statevar[_qp] = statevar_tplusdt;

    _interface_traction[_qp] = 0.0;
    _dinterface_traction_djump[_qp] = 0;   

}