//Implementation of Rate-and-State Friction Law in 2D

#include "RateStateFrictionLaw2DAsBCdev.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", RateStateFrictionLaw2DAsBCdev);

InputParameters
RateStateFrictionLaw2DAsBCdev::validParams()
{ 
  //Tn_o,Ts_o,Vini,statevarini are defined in "CZMComputeLocalTractionTotalBaseRSF2D"
  InputParameters params = CZMComputeLocalTractionTotalBaseRSF2D::validParams();
  params.addClassDescription("Rate-and-State Frictional Law.");
  params.addRequiredParam<Real>("f_o","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_a","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_b","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_L","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("delta_o","slip rate parameter");
  params.addRequiredParam<int>("RSFlaw","rate-and-state law options 1 (RS-A) : againg law + modified form; 2 (RS-S) : slip law + variant form");
  params.addParam<Real>("f_w",0.0,"weakened state friction coefficient");
  params.addParam<Real>("Vw",0.0,"characteristic weakening velocity");
//   params.addParam<Real>("elem_size", 1.0, "Value of element size");
//   params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("reaction_rsf_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_rsf_y","reaction in y dir");
  params.addRequiredCoupledVar("jacob_x","reaction in x dir");
  params.addRequiredCoupledVar("jacob_y","reaction in y dir");
//   params.addRequiredCoupledVar("reaction_rsf_pressure_x","reaction in x dir");
//   params.addRequiredCoupledVar("reaction_rsf_pressure_y","reaction in y dir");
//   params.addRequiredCoupledVar("jacob_pressure_x","reaction in x dir");
//   params.addRequiredCoupledVar("jacob_pressure_y","reaction in y dir");
//   params.addRequiredCoupledVar("interface_pressure","Pressure at sides of the fault");
//   params.addRequiredParam<std::string>("permeability_type", 
//     "Type of permeability condition (permeable/impermeable/semi_permeable)");
  params.addRequiredCoupledVar("reaction_damp_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_damp_y","reaction in y dir");
  params.addRequiredCoupledVar("jacob_damp_x","reaction in x dir");
  params.addRequiredCoupledVar("jacob_damp_y","reaction in y dir");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir");
  return params;
}

RateStateFrictionLaw2DAsBCdev::RateStateFrictionLaw2DAsBCdev(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseRSF2D(parameters),
    _len(coupledValue("nodal_area")),
    // _len2(getParam<Real>("elem_size")),
    _f_o(getParam<Real>("f_o")),
    _rsf_a(getParam<Real>("rsf_a")),
    _rsf_b(getParam<Real>("rsf_b")),
    _rsf_L(getParam<Real>("rsf_L")),
    _delta_o(getParam<Real>("delta_o")),
    _RSFlaw(getParam<int>("RSFlaw")),
    _f_w(getParam<Real>("f_w")),
    _Vw(getParam<Real>("Vw")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rhof(getMaterialPropertyByName<Real>(_base_name + "rhof")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _reaction_rsf_x(coupledValue("reaction_rsf_x")),
    _reaction_rsf_y(coupledValue("reaction_rsf_y")),
    _reaction_rsf_neighbor_x(coupledNeighborValue("reaction_rsf_x")),
    _reaction_rsf_neighbor_y(coupledNeighborValue("reaction_rsf_y")),
    _jacobian_rsf_x(coupledValue("jacob_x")),
    _jacobian_rsf_y(coupledValue("jacob_y")),
    _jacobian_rsf_neighbor_x(coupledNeighborValue("jacob_x")),
    _jacobian_rsf_neighbor_y(coupledNeighborValue("jacob_y")),
    // _reaction_rsf_pressure_x(coupledValue("reaction_rsf_pressure_x")),
    // _reaction_rsf_pressure_y(coupledValue("reaction_rsf_pressure_y")),
    // _reaction_rsf_neighbor_pressure_x(coupledNeighborValue("reaction_rsf_pressure_x")),
    // _reaction_rsf_neighbor_pressure_y(coupledNeighborValue("reaction_rsf_pressure_y")),
    // _jacobian_rsf_pressure_x(coupledValue("jacob_pressure_x")),
    // _jacobian_rsf_pressure_y(coupledValue("jacob_pressure_y")),
    // _jacobian_rsf_neighbor_pressure_x(coupledNeighborValue("jacob_pressure_x")),
    // _jacobian_rsf_neighbor_pressure_y(coupledNeighborValue("jacob_pressure_y")),
    // _permeability_type(getParam<std::string>("permeability_type")),
    // _interface_pressure_plus(coupledNeighborValue("interface_pressure")),
    // _interface_pressure_minus(coupledValue("interface_pressure")),
    _reaction_damp_x(coupledValue("reaction_damp_x")),
    _reaction_damp_y(coupledValue("reaction_damp_y")),
    _reaction_damp_neighbor_x(coupledNeighborValue("reaction_damp_x")),
    _reaction_damp_neighbor_y(coupledNeighborValue("reaction_damp_y")),
    _jacobian_damp_x(coupledValue("jacob_damp_x")),
    _jacobian_damp_y(coupledValue("jacob_damp_y")),
    _jacobian_damp_neighbor_x(coupledNeighborValue("jacob_damp_x")),
    _jacobian_damp_neighbor_y(coupledNeighborValue("jacob_damp_y")),
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
    _alongfaultdisp_normal_plus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_strike_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_minus_older(getMaterialPropertyOlderByName<Real>("alongfaultdisp_normal_minus")),
    _sliprate_strike_old(getMaterialPropertyOldByName<Real>("sliprate_strike")),
    _sliprate_normal_old(getMaterialPropertyOldByName<Real>("sliprate_normal")),
    _sliprate_mag_old(getMaterialPropertyOldByName<Real>("sliprate_mag")),
    _sliprate_predict(getMaterialPropertyOldByName<Real>("sliprate_predict")),
    _slip_strike_old(getMaterialPropertyOldByName<Real>("slip_strike")),
    _slip_normal_old(getMaterialPropertyOldByName<Real>("slip_normal")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _statevar_older(getMaterialPropertyOlderByName<Real>("statevar")),
    _traction_strike_old(getMaterialPropertyOldByName<Real>("traction_strike")),
    _traction_normal_old(getMaterialPropertyOldByName<Real>("traction_normal"))

{
}

void
RateStateFrictionLaw2DAsBCdev::computeInterfaceTractionAndDerivatives()
{   
    //Check at the beginning
    if ( _t == 0.0 ) {
        
        //Check if fw and Vw is not provided in _RSFlaw == 2
        if ( _RSFlaw == 2 and ( _f_w == 0.0 or _Vw == 0.0 ) ){
            mooseError(" Must provide valid fw, Vw for RSFlaw == 2 !");
        }
        if ( _RSFlaw == 1 and ( _f_w != 0.0 or _Vw != 0.0 ) ){
            mooseError(" Wrong RSFlaw for nonzero fw, Vw !");
        }
    }
    
    //Define Parameters
    Real len = _len[_qp];
    Real f_o = _f_o;
    Real rsf_a = _rsf_a;
    Real rsf_b = _rsf_b;
    Real rsf_L = _rsf_L;
    Real delta_o = _delta_o;
    Real Tn_o = _Tn_o;
    Real Ts_o = _Ts_o;
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal

    // //pressure state 
    // Real alpha = _biot_coefficient[_qp];
    // Real p_plus = _interface_pressure_plus[_qp];
    // Real p_minus = _interface_pressure_minus[_qp];
    // Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
    // Real p_old = std::max(_interface_pressure_plus_older[_qp], _interface_pressure_minus_older[_qp]);
    
    // // Define pressure based on fault permeability condition
    // if (_permeability_type == "permeable")
    // {
    //     Real p = 0.0;  // For permeable condition
    //     Real p_old = 0.0;
    // }

    //*Restoration Force*

    //Effective stress = Stress Divergence - pressure divergence Components (label as stsdivcomp)

    // ///Define in global coordinate
    // //current time step 
    // RealVectorValue R_plus_global_prsdivcomp(-_reaction_rsf_pressure_x[_qp],-_reaction_rsf_pressure_y[_qp], 0.0);
    // RealVectorValue R_minus_global_prsdivcomp(-_reaction_rsf_neighbor_pressure_x[_qp],-_reaction_rsf_neighbor_pressure_y[_qp], 0.0);

    // ///Rotate in local coordinate
    // //current time step
    // RealVectorValue R_plus_local_prsdivcomp = _rot[_qp].transpose() * R_plus_global_prsdivcomp;
    // RealVectorValue R_minus_local_prsdivcomp = _rot[_qp].transpose() * R_minus_global_prsdivcomp;

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_stsdivcomp(-_reaction_rsf_x[_qp],-_reaction_rsf_y[_qp], 0.0);
    RealVectorValue R_minus_global_stsdivcomp(-_reaction_rsf_neighbor_x[_qp],-_reaction_rsf_neighbor_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_stsdivcomp = _rot[_qp].transpose() * R_plus_global_stsdivcomp;
    RealVectorValue R_minus_local_stsdivcomp = _rot[_qp].transpose() * R_minus_global_stsdivcomp;

    ///Get Components
    //current time step
    Real R_plus_local_normal_stsdivcomp  = R_plus_local_stsdivcomp(0); //+ R_plus_local_prsdivcomp(0);
    Real R_plus_local_strike_stsdivcomp  = R_plus_local_stsdivcomp(1); //+ R_plus_local_prsdivcomp(1);
    
    Real R_minus_local_normal_stsdivcomp = R_minus_local_stsdivcomp(0); //+ R_minus_local_prsdivcomp(0);
    Real R_minus_local_strike_stsdivcomp = R_minus_local_stsdivcomp(1); //+ R_minus_local_prsdivcomp(1);

    //--------------------------------------------------------------------------------------------------//

    //Damping Components Contribution (label as dampingcomp)

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_dampingcomp(-_reaction_damp_x[_qp],-_reaction_damp_y[_qp], 0.0);
    RealVectorValue R_minus_global_dampingcomp(-_reaction_damp_neighbor_x[_qp],-_reaction_damp_neighbor_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_dampingcomp = _rot[_qp].transpose() * R_plus_global_dampingcomp;
    RealVectorValue R_minus_local_dampingcomp = _rot[_qp].transpose() * R_minus_global_dampingcomp;

    ///Get Components
    //current time step
    Real R_plus_local_normal_dampingcomp  = R_plus_local_dampingcomp(0);
    Real R_plus_local_strike_dampingcomp  = R_plus_local_dampingcomp(1);
    
    Real R_minus_local_normal_dampingcomp = R_minus_local_dampingcomp(0);
    Real R_minus_local_strike_dampingcomp = R_minus_local_dampingcomp(1);

    //--------------------------------------------------------------------------------------------------//

    //Add restoration forces from two contributions
    Real R_plus_local_normal  = R_plus_local_normal_stsdivcomp  + R_plus_local_normal_dampingcomp; 
    Real R_plus_local_strike  = R_plus_local_strike_stsdivcomp  + R_plus_local_strike_dampingcomp;
    Real R_minus_local_normal = R_minus_local_normal_stsdivcomp + R_minus_local_normal_dampingcomp; 
    Real R_minus_local_strike = R_minus_local_strike_stsdivcomp + R_minus_local_strike_dampingcomp; 

    //--------------------------------------------------------------------------------------------------//

    //*jacobian Restoration Force*

    //Effective stress = Stress Divergence - pressure divergence Components (label as stsdivcomp)

    // ///Define in global coordinate
    // //current time step 
    // RealVectorValue J_plus_global_prsdivcomp(-_jacobian_rsf_pressure_x[_qp],-_jacobian_rsf_pressure_y[_qp], 0.0);
    // RealVectorValue J_minus_global_prsdivcomp(-_jacobian_rsf_neighbor_pressure_x[_qp],-_jacobian_rsf_neighbor_pressure_y[_qp], 0.0);

    // ///Rotate in local coordinate
    // //current time step
    // RealVectorValue J_plus_local_prsdivcomp = _rot[_qp].transpose() * J_plus_global_prsdivcomp;
    // RealVectorValue J_minus_local_prsdivcomp = _rot[_qp].transpose() * J_minus_global_prsdivcomp;

    ///Define in global coordinate
    //current time step 
    RealVectorValue J_plus_global_stsdivcomp(-_jacobian_rsf_x[_qp],-_jacobian_rsf_y[_qp], 0.0);
    RealVectorValue J_minus_global_stsdivcomp(-_jacobian_rsf_neighbor_x[_qp],-_jacobian_rsf_neighbor_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue J_plus_local_stsdivcomp = _rot[_qp].transpose() * J_plus_global_stsdivcomp;
    RealVectorValue J_minus_local_stsdivcomp = _rot[_qp].transpose() * J_minus_global_stsdivcomp;

    ///Get Components
    //current time step
    Real J_plus_local_normal_stsdivcomp  = J_plus_local_stsdivcomp(0) ; //+ J_plus_local_prsdivcomp(0)
    Real J_plus_local_strike_stsdivcomp  = J_plus_local_stsdivcomp(1) ; // J_plus_local_prsdivcomp(1);
    
    Real J_minus_local_normal_stsdivcomp = J_minus_local_stsdivcomp(0); // J_minus_local_prsdivcomp(0);
    Real J_minus_local_strike_stsdivcomp = J_minus_local_stsdivcomp(1); // J_minus_local_prsdivcomp(1);

    //--------------------------------------------------------------------------------------------------//

    //Damping Components Contribution (label as dampingcomp)

    ///Define in global coordinate
    //current time step 
    RealVectorValue J_plus_global_dampingcomp(-_jacobian_damp_x[_qp],-_jacobian_damp_y[_qp], 0.0);
    RealVectorValue J_minus_global_dampingcomp(-_jacobian_damp_neighbor_x[_qp],-_jacobian_damp_neighbor_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue J_plus_local_dampingcomp = _rot[_qp].transpose() * J_plus_global_dampingcomp;
    RealVectorValue J_minus_local_dampingcomp = _rot[_qp].transpose() * J_minus_global_dampingcomp;

    ///Get Components
    //current time step
    Real J_plus_local_normal_dampingcomp  = J_plus_local_dampingcomp(0);
    Real J_plus_local_strike_dampingcomp  = J_plus_local_dampingcomp(1);
    
    Real J_minus_local_normal_dampingcomp = J_minus_local_dampingcomp(0);
    Real J_minus_local_strike_dampingcomp = J_minus_local_dampingcomp(1);

    //--------------------------------------------------------------------------------------------------//

    //Add restoration forces from two contributions
    Real J_plus_local_normal  = J_plus_local_normal_stsdivcomp  + J_plus_local_normal_dampingcomp; 
    Real J_plus_local_strike  = J_plus_local_strike_stsdivcomp  + J_plus_local_strike_dampingcomp; 
    Real J_minus_local_normal = J_minus_local_normal_stsdivcomp + J_minus_local_normal_dampingcomp; 
    Real J_minus_local_strike = J_minus_local_strike_stsdivcomp + J_minus_local_strike_dampingcomp; 

    //--------------------------------------------------------------------------------------------------//

    //*Nodal Mass*
    ///QUAD4 Element
    Real M = _density[_qp] * len * len * 0.5;  
    Real Mf = _rhof[_qp] * len * len * 0.5;

    //Initialize shear stress perturbation
    Real Ts_perturb = _Ts_perturb[_qp];

    ///Disp Plus Side
    Real alongfaultdisp_strike_plus_t = _alongfaultdisp_strike_plus_old[_qp];
    Real alongfaultdisp_normal_plus_t = _alongfaultdisp_normal_plus_old[_qp];
    ///Disp Minus Side
    Real alongfaultdisp_strike_minus_t = _alongfaultdisp_strike_minus_old[_qp];
    Real alongfaultdisp_normal_minus_t = _alongfaultdisp_normal_minus_old[_qp];

    //Older Value
    //Plus Side
    Real alongfaultdisp_strike_plus_tminust = _alongfaultdisp_strike_plus_older[_qp];
    Real alongfaultdisp_normal_plus_tminust = _alongfaultdisp_normal_plus_older[_qp];
    //Minus Side
    Real alongfaultdisp_strike_minus_tminust = _alongfaultdisp_strike_minus_older[_qp];
    Real alongfaultdisp_normal_minus_tminust = _alongfaultdisp_normal_minus_older[_qp];

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
    Tn = abs(Tn); // - p;

    //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
    Real Ts_trial = ( M * M * sliprate_strike_tminusdtover2 )/( len * _dt * (M + M) ) + (M * R_plus_local_strike - M * R_minus_local_strike) / ( len * ( M + M ) ) + Ts_o + Ts_perturb;
    Real Tmag_trial = sqrt(Ts_trial*Ts_trial);

    Real Z = 0.0;
    switch ( _RSFlaw )
    {
        case 1:

            Z = 0.5 / delta_o * exp((f_o + rsf_b * log(delta_o * statevar_t/rsf_L))/rsf_a); 
            break;

        case 2:

            Z = 0.5 / delta_o * exp( statevar_t / rsf_a );
            break; 

        default:

            mooseError("Must specify a valid RSFlaw Parameter!");

    }

    //const
    Real c = len * _dt * ( M + M ) / (M * M);
    //Real Z = 0.5 / delta_o * exp((f_o + rsf_b * log(delta_o * statevar_t/rsf_L))/rsf_a);    
    
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
    // theta_ref = ( theta_pre + _dt ) / ( 1.0 + _dt * dv_pre / rsf_L ); //febealg2dot1
    //theta_ref = ( _statevar_older[_qp] + 2 * _dt ) / ( 1.0 + 2 * _dt * dv_pre / rsf_L );

    Real statevar_ss = 0.0;
    Real f_LV = 0.0; //low-velocity steady state friction coefficient
    Real f_ss = 0.0; //steady-state friction coefficient
    switch (_RSFlaw)
    {
        case 1:

            statevar_ss = (rsf_L/sliprate_mag_tplusdtover2);
            break;

        case 2:
            
            f_LV = f_o - ( rsf_b - rsf_a ) * log( sliprate_mag_tplusdtover2 / delta_o );
            f_ss = _f_w + ( f_LV - _f_w ) / pow( 1 + pow( sliprate_mag_tplusdtover2 / _Vw , 8 ) , 0.125 );
            statevar_ss = rsf_a * log( 2 * delta_o / sliprate_mag_tplusdtover2 * sinh( f_ss / rsf_a ) );
            break;

        default:

            mooseError("Must specify a valid RSFlaw Parameter!");
        
    }

    Real coeffD = exp(-sliprate_mag_tplusdtover2*_dt/rsf_L);
    Real statevar_tplusdt = statevar_t * coeffD + statevar_ss * (1-coeffD);

    //*Compute shear traction at time t*
    Real T_mag = Tn * rsf_a * asinh( 0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2) * Z );

    ///Get Components
    Ts = T_mag * ( Ts_trial / Tmag_trial );

    ////DISP
    Real du_strike_plus_t  =  alongfaultdisp_strike_plus_t  - alongfaultdisp_strike_plus_tminust  + _dt * _dt / M * (R_plus_local_strike  - len * (Ts - Ts_o - Ts_perturb));
    Real du_normal_plus_t  =  alongfaultdisp_normal_plus_t  - alongfaultdisp_normal_plus_tminust  + _dt * _dt / M * (R_plus_local_normal  - len * (Tn - Tn_o));
    
    Real du_strike_minus_t =  alongfaultdisp_strike_minus_t  - alongfaultdisp_strike_minus_tminust + _dt * _dt / M * (R_minus_local_strike  + len * (Ts - Ts_o - Ts_perturb));
    Real du_normal_minus_t =  alongfaultdisp_normal_minus_t  - alongfaultdisp_normal_minus_tminust + _dt * _dt / M * (R_minus_local_normal  + len * (Tn - Tn_o));

    Real alongfaultdisp_strike_plus_tplusdt  = alongfaultdisp_strike_plus_t  + du_strike_plus_t;
    Real alongfaultdisp_normal_plus_tplusdt  = alongfaultdisp_normal_plus_t  + du_normal_plus_t;
    
    Real alongfaultdisp_strike_minus_tplusdt = alongfaultdisp_strike_minus_t + du_strike_minus_t;
    Real alongfaultdisp_normal_minus_tplusdt = alongfaultdisp_normal_minus_t + du_normal_minus_t;

    //Update Slip Rate and Slip at t + dt/2 or t + dt
    ///Slip Rate
    Real sliprate_strike_tplusdtover2 = sliprate_mag_tplusdtover2 * ( Ts_trial / Tmag_trial );
    Real sliprate_normal_tplusdtover2 = ( alongfaultdisp_normal_plus_tplusdt - alongfaultdisp_normal_plus_t ) / _dt - ( alongfaultdisp_normal_minus_tplusdt - alongfaultdisp_normal_minus_t ) / _dt;

    ///Slip
    Real slip_strike_tplusdtover2 = slip_strike_t + sliprate_strike_tplusdtover2 * _dt;
    Real slip_normal_tplusdtover2 = alongfaultdisp_normal_plus_tplusdt - alongfaultdisp_normal_minus_tplusdt;

    //Save Results
    ///Traction (t)
    _traction_strike[_qp] = Ts;
    _traction_normal[_qp] = Tn;

    ///Disp (t+dt)
    _alongfaultdisp_strike_plus[_qp] = alongfaultdisp_strike_plus_tplusdt;
    _alongfaultdisp_normal_plus[_qp] = alongfaultdisp_normal_plus_tplusdt;

    _alongfaultdisp_strike_minus[_qp] = alongfaultdisp_strike_minus_tplusdt;
    _alongfaultdisp_normal_minus[_qp] = alongfaultdisp_normal_minus_tplusdt;

    ///Slip Rate (t+dt/2)
    _sliprate_strike[_qp] = sliprate_strike_tplusdtover2;
    _sliprate_normal[_qp] = sliprate_normal_tplusdtover2;
    _sliprate_mag[_qp]    = sliprate_mag_tplusdtover2;

    ///Slip (t+dt)
    _slip_strike[_qp] = slip_strike_tplusdtover2;
    _slip_normal[_qp] = slip_normal_tplusdtover2;
    _slip_mag[_qp]    = sqrt(slip_strike_tplusdtover2*slip_strike_tplusdtover2);

    ///State Variable (t+dt)
    _statevar[_qp] = statevar_tplusdt;

    _interface_traction[_qp] = 0.0;
    _dinterface_traction_djump[_qp] = 0;   

    // //--------------------------------------------------------------------------------------------------//

    //Jacobians of the displacements
   
    // Diagonal and off diagonal jacobians with respect to displacements on both sides 
    // (This should consider the stress divergence, pressure divergence and damping jacobians)

    Real dR_plus_local_strike_du_strike_plus = J_plus_local_strike;
    Real dR_minus_local_strike_du_strike_plus = 0;
    Real dR_plus_local_normal_du_strike_plus  = 0;
    Real dR_minus_local_normal_du_strike_plus = 0;

    Real dR_plus_local_strike_du_strike_minus = 0;
    Real dR_minus_local_strike_du_strike_minus = J_minus_local_strike;
    Real dR_plus_local_normal_du_strike_minus  = 0;
    Real dR_minus_local_normal_du_strike_minus = 0;

    Real dR_plus_local_strike_du_normal_plus = 0;
    Real dR_minus_local_strike_du_normal_plus = 0;
    Real dR_plus_local_normal_du_normal_plus  = J_plus_local_normal;
    Real dR_minus_local_normal_du_normal_plus = 0;

    Real dR_plus_local_strike_du_normal_minus = 0;
    Real dR_minus_local_strike_du_normal_minus = 0;
    Real dR_plus_local_normal_du_normal_minus  = 0;
    Real dR_minus_local_normal_du_normal_minus = J_minus_local_normal;


    Real dTn_du_strike_plus =  0;
    Real dTn_du_strike_minus =  0;
    Real dTn_du_normal_plus  =  0;
    Real dTn_du_normal_minus =  0;

    // jacobians of the normal traction   
    if (Tn_trial<0)
    {
     // if we consider the old slip and slip rates

        Real dTn_du_strike_plus =  ( M * dR_minus_local_normal_du_strike_plus - M * dR_plus_local_normal_du_strike_plus  ) / ( len * (M + M) ) ;
        Real dTn_du_strike_minus =  ( M * dR_minus_local_normal_du_strike_minus - M * dR_plus_local_normal_du_strike_minus  ) / ( len * (M + M) ) ;

        Real dTn_du_normal_plus  =  (-1.0 * (1.0/_dt) * M * M * ( (1.0/_dt) * 1 + (1.0/_dt) * 1) ) / ( len * (M + M) ) + ( M * dR_minus_local_normal_du_normal_plus - M * dR_plus_local_normal_du_normal_plus ) / ( len * (M + M) );
        Real dTn_du_normal_minus =  (-1.0 * (1.0/_dt) * M * M * ( (1.0/_dt) * -1 + (1.0/_dt) * -1) ) / ( len * (M + M) ) + ( M * dR_minus_local_normal_du_normal_minus - M * dR_plus_local_normal_du_normal_minus ) / ( len * (M + M) );

     // if we don't consider the old slip and slip rates

        // dTn_du_strike_plus =  ( M * dR_minus_local_normal_du_strike_plus - M * dR_plus_local_normal_du_strike_plus  ) / ( len * (M + M) ) ;
        // dTn_du_strike_minus =  ( M * dR_minus_local_normal_du_strike_minus - M * dR_plus_local_normal_du_strike_minus  ) / ( len * (M + M) ) ;

        // dTn_du_normal_plus  =  ( M * dR_minus_local_normal_du_normal_plus - M * dR_plus_local_normal_du_normal_plus ) / ( len * (M + M) );
        // dTn_du_normal_minus =   ( M * dR_minus_local_normal_du_normal_minus - M * dR_plus_local_normal_du_normal_minus ) / ( len * (M + M) );

    }else{
        
    }


    Real dstatevar_ss_du_strike_plus = 0;
    Real dstatevar_ss_du_strike_minus = 0;
    Real dstatevar_ss_du_normal_plus = 0;
    Real dstatevar_ss_du_normal_minus = 0;

    Real dZ_du_strike_plus = 0;
    Real dZ_du_strike_minus = 0;
    Real dZ_du_normal_plus = 0;
    Real dZ_du_normal_minus = 0;

    // Then handle assignments in switch
    switch (_RSFlaw)
    {
        case 1:
            dstatevar_ss_du_strike_plus = -(rsf_L / pow(sliprate_mag_tplusdtover2, 2)) * (1 / _dt);
            dstatevar_ss_du_strike_minus = (rsf_L / pow(sliprate_mag_tplusdtover2, 2)) * (1 / _dt);
            // dstatevar_ss_du_normal_plus and dstatevar_ss_du_normal_minus stay 0

            dZ_du_strike_plus = (Z / statevar_t) * (rsf_b / rsf_a) * dstatevar_ss_du_strike_plus;
            dZ_du_strike_minus = (Z / statevar_t) * (rsf_b / rsf_a) * dstatevar_ss_du_strike_minus;
            dZ_du_normal_plus = (Z / statevar_t) * (rsf_b / rsf_a) * dstatevar_ss_du_normal_plus;
            dZ_du_normal_minus = (Z / statevar_t) * (rsf_b / rsf_a) * dstatevar_ss_du_normal_minus;
            break;

        case 2:
            // All variables already initialized to 0
            break;

        default:
            mooseError("Must specify a valid RSFlaw Parameter!");
    }  

 

       Real C =  0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2);

       Real dT_mag_du_strike_plus = dTn_du_strike_plus * rsf_a * asinh(C*Z)
                                  + Tn * rsf_a * dZ_du_strike_plus * C/std::sqrt(1+std::pow(C*Z,2))
                                  + Tn * rsf_a * (1.0/_dt) * Z * C/std::sqrt(1+std::pow(C*Z,2)); 
       Real dT_mag_du_strike_minus = dTn_du_strike_minus * rsf_a * asinh(C*Z)
                                   + Tn * rsf_a * dZ_du_strike_minus * C/std::sqrt(1+std::pow(C*Z,2))
                                   + Tn * rsf_a * (-1.0/_dt) * Z * C/std::sqrt(1+std::pow(C*Z,2));    
       Real dT_mag_du_normal_plus = dTn_du_normal_plus * rsf_a * asinh(C*Z)
                                  + Tn * rsf_a * dZ_du_normal_plus * C/std::sqrt(1+std::pow(C*Z,2));                              
       Real dT_mag_du_normal_minus = dTn_du_normal_minus * rsf_a * asinh(C*Z)
                                   + Tn * rsf_a * dZ_du_normal_minus * C/std::sqrt(1+std::pow(C*Z,2));
                                                   

       Real dTs_du_strike_plus = dT_mag_du_strike_plus * ( Ts_trial / Tmag_trial );
       Real dTs_du_strike_minus = dT_mag_du_strike_minus * ( Ts_trial / Tmag_trial );
       Real dTs_du_normal_plus = dT_mag_du_normal_plus * ( Ts_trial / Tmag_trial );
       Real dTs_du_normal_minus = dT_mag_du_normal_minus * ( Ts_trial / Tmag_trial );

        Real dalongfaultdisp_strike_plus_tplusdt_du_strike_plus = 1 + 1 + _dt * _dt / M * (dR_plus_local_strike_du_strike_plus - len * dTs_du_strike_plus);
        Real dalongfaultdisp_strike_plus_tplusdt_du_strike_minus = _dt * _dt / M * (dR_plus_local_strike_du_strike_minus - len * dTs_du_strike_minus);
        Real dalongfaultdisp_strike_plus_tplusdt_du_normal_plus = _dt * _dt / M * (dR_plus_local_strike_du_normal_plus - len * dTs_du_normal_plus);
        Real dalongfaultdisp_strike_plus_tplusdt_du_normal_minus = _dt * _dt / M * (dR_plus_local_strike_du_normal_minus - len * dTs_du_normal_minus);
        
        Real dalongfaultdisp_strike_minus_tplusdt_du_strike_plus = _dt * _dt / M * (dR_minus_local_strike_du_strike_plus + len * dTs_du_strike_plus);
        Real dalongfaultdisp_strike_minus_tplusdt_du_strike_minus = 1 + 1 + _dt * _dt / M * (dR_minus_local_strike_du_strike_minus + len * dTs_du_strike_minus);
        Real dalongfaultdisp_strike_minus_tplusdt_du_normal_plus = _dt * _dt / M * (dR_minus_local_strike_du_normal_plus + len * dTs_du_normal_plus);
        Real dalongfaultdisp_strike_minus_tplusdt_du_normal_minus = _dt * _dt / M * (dR_minus_local_strike_du_normal_minus + len * dTs_du_normal_minus);

        Real dalongfaultdisp_normal_plus_tplusdt_du_strike_plus = _dt * _dt / M * (dR_plus_local_normal_du_strike_plus - len * dTn_du_strike_plus);
        Real dalongfaultdisp_normal_plus_tplusdt_du_strike_minus = _dt * _dt / M * (dR_plus_local_normal_du_strike_minus - len * dTn_du_strike_minus);
        Real dalongfaultdisp_normal_plus_tplusdt_du_normal_plus = 1 + 1 + _dt * _dt / M * (dR_plus_local_normal_du_normal_plus - len * dTn_du_normal_plus);
        Real dalongfaultdisp_normal_plus_tplusdt_du_normal_minus = _dt * _dt / M * (dR_plus_local_normal_du_normal_minus - len * dTn_du_normal_minus);

        Real dalongfaultdisp_normal_minus_tplusdt_du_strike_plus = _dt * _dt / M * (dR_minus_local_normal_du_strike_plus + len * dTn_du_strike_plus);
        Real dalongfaultdisp_normal_minus_tplusdt_du_strike_minus = _dt * _dt / M * (dR_minus_local_normal_du_strike_minus + len * dTn_du_strike_minus);
        Real dalongfaultdisp_normal_minus_tplusdt_du_normal_plus = _dt * _dt / M * (dR_minus_local_normal_du_normal_plus + len * dTn_du_normal_plus);
        Real dalongfaultdisp_normal_minus_tplusdt_du_normal_minus = 1 + 1 + _dt * _dt / M * (dR_minus_local_normal_du_normal_minus + len * dTn_du_normal_minus);
 
        // _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus[_qp] = dalongfaultdisp_strike_plus_tplusdt_du_strike_plus;
        // _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus[_qp] = dalongfaultdisp_strike_plus_tplusdt_du_strike_minus;
        // _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus[_qp] = dalongfaultdisp_strike_plus_tplusdt_du_normal_plus;
        // _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus[_qp] = dalongfaultdisp_strike_plus_tplusdt_du_normal_minus;
        // _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus[_qp] = dalongfaultdisp_strike_minus_tplusdt_du_strike_plus;
        // _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus[_qp] = dalongfaultdisp_strike_minus_tplusdt_du_strike_minus;
        // _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus[_qp] = dalongfaultdisp_strike_minus_tplusdt_du_normal_plus;
        // _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus[_qp] = dalongfaultdisp_strike_minus_tplusdt_du_normal_minus;
        // _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus[_qp] = dalongfaultdisp_normal_plus_tplusdt_du_strike_plus;
        // _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus[_qp] = dalongfaultdisp_normal_plus_tplusdt_du_strike_minus;
        // _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus[_qp] = dalongfaultdisp_normal_plus_tplusdt_du_normal_plus;
        // _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus[_qp] = dalongfaultdisp_normal_plus_tplusdt_du_normal_minus;
        // _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus[_qp] = dalongfaultdisp_normal_minus_tplusdt_du_strike_plus;
        // _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus[_qp] = dalongfaultdisp_normal_minus_tplusdt_du_strike_minus;
        // _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus[_qp] = dalongfaultdisp_normal_minus_tplusdt_du_normal_plus;
        // _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus[_qp] = dalongfaultdisp_normal_minus_tplusdt_du_normal_minus;

        _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus[_qp] = 1;
        _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus[_qp] = 0;
        _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus[_qp] = 0;
        _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus[_qp] = 0;
        _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus[_qp] = 0;
        _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus[_qp] = 1;
        _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus[_qp] = 0;
        _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus[_qp] = 0;
        _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus[_qp] = 0;
        _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus[_qp] = 0;
        _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus[_qp] = 1;
        _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus[_qp] = 0;
        _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus[_qp] = 0;
        _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus[_qp] = 0;
        _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus[_qp] = 0;
        _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus[_qp] = 1;
 
 }

