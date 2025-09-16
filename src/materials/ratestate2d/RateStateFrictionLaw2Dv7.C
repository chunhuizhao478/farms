//Implementation of Rate-and-State Friction Law in 2D

#include "RateStateFrictionLaw2Dv7.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"

registerMooseObject("farmsApp", RateStateFrictionLaw2Dv7);

InputParameters
RateStateFrictionLaw2Dv7::validParams()
{ 
  //Tn_o,Ts_o,Vini,statevarini are defined in "CZMComputeLocalTractionTotalBaseRSF2D"
  InputParameters params = CZMComputeLocalTractionTotalBaseRSF2D::validParams();
  params.addClassDescription("Rate-and-State Frictional Law.");
  params.addParam<Real>("len",100,"element length");
  params.addParam<Real>("f_o",0.6,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_a",0.008,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_b",0.012,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_L",0.02,"rate-and-state friction coefficients");
  params.addParam<Real>("delta_o",1e-6,"slip rate parameter");
  params.addRequiredCoupledVar("reaction_rsf_x","reaction in x dir at time t");
  params.addRequiredCoupledVar("reaction_rsf_y","reaction in y dir at time t");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir at time t");
  return params;
}

RateStateFrictionLaw2Dv7::RateStateFrictionLaw2Dv7(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseRSF2D(parameters),
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
    _reaction_rsf_neighbor_x(coupledNeighborValue("reaction_rsf_x")),
    _reaction_rsf_neighbor_y(coupledNeighborValue("reaction_rsf_y")),
    _Ts_perturb(coupledValue("Ts_perturb")),
    _Ts_perturb_old(coupledValueOld("Ts_perturb")),
    _alongfaultvel_strike_plus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_strike_plus")),
    _alongfaultvel_strike_minus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_strike_minus")),
    _alongfaultvel_normal_plus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_normal_plus")),
    _alongfaultvel_normal_minus_old(getMaterialPropertyOldByName<Real>("alongfaultvel_normal_minus")),
    _alongfaultdisp_strike_plus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_strike_plus")),
    _alongfaultdisp_strike_minus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_plus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_normal_minus_old(getMaterialPropertyOldByName<Real>("alongfaultdisp_normal_minus")),
    _alongfaultdisp_strike_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_strike_plus")),
    _alongfaultdisp_strike_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_normal_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_minus")),
    _sliprate_strike_old(getMaterialPropertyOldByName<Real>("sliprate_strike")),
    _sliprate_normal_old(getMaterialPropertyOldByName<Real>("sliprate_normal")),
    _sliprate_mag_old(getMaterialPropertyOldByName<Real>("sliprate_mag")),
    _slip_strike_old(getMaterialPropertyOldByName<Real>("slip_strike")),
    _slip_normal_old(getMaterialPropertyOldByName<Real>("slip_normal")),
    _slip_mag_old(getMaterialPropertyOldByName<Real>("slip_mag")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _traction_strike_old(getMaterialPropertyOldByName<Real>("traction_strike")),
    _traction_normal_old(getMaterialPropertyOldByName<Real>("traction_normal"))
{
}

void
RateStateFrictionLaw2Dv7::computeInterfaceTractionAndDerivatives()
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
    
    //--------------------------------------------------------------------------------------------------//

    //*Restoration Force*
    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global(-_reaction_rsf_x[_qp],-_reaction_rsf_y[_qp],0);
    RealVectorValue R_minus_global(-_reaction_rsf_neighbor_x[_qp],-_reaction_rsf_neighbor_y[_qp],0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
    RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

    ///Get Components
    //current time step
    Real R_plus_local_strike  = R_plus_local(1);
    Real R_plus_local_normal  = R_plus_local(0);
    Real R_minus_local_strike = R_minus_local(1);
    Real R_minus_local_normal = R_minus_local(0);

    //--------------------------------------------------------------------------------------------------//

    //*Nodal Mass*
    ///QUAD4 Element
    Real M = _density[_qp] * len * len * 0.5;

    //Initialize Final Traction
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal

    //Initialize Shear Perturbation
    Real Ts_perturb = _Ts_perturb[_qp]; 

    //*Get slip rate and slip from time t-dt/2 and t
    ///strike direction
    Real sliprate_strike_tminusdtover2 = _sliprate_strike_old[_qp];
    Real slip_strike_t = _slip_strike_old[_qp];

    ///normal direction
    Real sliprate_normal_tminusdtover2 = _sliprate_normal_old[_qp];
    Real slip_normal_t = _slip_normal_old[_qp];

    ///mag
    Real sliprate_mag_tminusdtover2 = _sliprate_mag_old[_qp];

    //Get State Variable at Current Time Step
    Real statevar_t = _statevar_old[_qp];

    //*Compute Trial Normal Traction and Normal Traction*
    Real Tn_trial = ( -1.0 * (1.0/_dt) * M * M * ( sliprate_normal_tminusdtover2 + (1.0/_dt) * slip_normal_t) ) / ( len * (M + M) ) + ( M * R_minus_local_normal - M * R_plus_local_normal ) / ( len * (M + M) ) - Tn_o;
    if (Tn_trial<0)
    {
        Tn = Tn_trial;
    }else{
        Tn = 0;
    }

    ///Make Tn positive
    Tn = abs(Tn);

    //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
    Real Ts_trial = ( M * M * sliprate_strike_tminusdtover2 )/( len * _dt * (M + M) ) + (M * R_plus_local_strike - M * R_minus_local_strike) / ( len * ( M + M ) ) + Ts_o + Ts_perturb;
    Real Tmag_trial = sqrt(Ts_trial*Ts_trial);

    //const
    Real c = len * _dt * ( M + M ) / (M * M);
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
            mooseError("NOT CONVERGED!"); //strong convergence check
        }

        //update iterr
        iterr = iterr + 1;
    
    }

    Real sliprate_mag_tplusdtover2 = abs(solution); 

    //update state variable
    Real coeffD = exp(-sliprate_mag_tplusdtover2*_dt/rsf_L);
    Real statevar_tplusdt = statevar_t * coeffD + (rsf_L/sliprate_mag_tplusdtover2) * (1-coeffD);

    //*Compute shear traction at time t*
    Real T_mag = Tn * rsf_a * asinh( 0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2) * Z );

    ///Get Components
    Ts = T_mag * ( Ts_trial / Tmag_trial );

    //Update Slip Rate and Slip at t + dt/2 or t + dt
    ///Slip Rate
    Real sliprate_strike_tplusdtover2 = sliprate_mag_tplusdtover2;
    Real sliprate_normal_tplusdtover2 = 0.0;
    Real sliprate_mag_tplusdtover2_updated = sqrt(sliprate_strike_tplusdtover2 * sliprate_strike_tplusdtover2);

    ///Slip
    Real slip_strike_tplusdtover2 = slip_strike_t + sliprate_strike_tplusdtover2 * _dt;
    Real slip_normal_tplusdtover2 = 0.0;
    Real slip_mag_tplusdtover2    = sqrt(slip_strike_tplusdtover2 * slip_strike_tplusdtover2);

    //Save Results
    ///Traction (t)
    _traction_strike[_qp] = Ts;
    _traction_normal[_qp] = Tn;

    ///Slip Rate (t+dt/2)
    _sliprate_strike[_qp] = sliprate_strike_tplusdtover2;
    _sliprate_normal[_qp] = sliprate_normal_tplusdtover2;
    _sliprate_mag[_qp]    = sliprate_mag_tplusdtover2_updated;

    ///Slip (t+dt)
    _slip_strike[_qp] = slip_strike_tplusdtover2;
    _slip_normal[_qp] = slip_normal_tplusdtover2;
    _slip_mag[_qp]    = slip_mag_tplusdtover2;

    ///State Variable (t+dt)
    _statevar[_qp] = statevar_tplusdt;

    //Feed traction into the system
    //Assign back traction in CZM
    RealVectorValue traction;

    traction(0) = 0.0; 
    traction(1) = Ts - Ts_o - Ts_perturb; 
    traction(2) = 0.0;

    _interface_traction[_qp] = traction;
    _dinterface_traction_djump[_qp] = 0;  

}