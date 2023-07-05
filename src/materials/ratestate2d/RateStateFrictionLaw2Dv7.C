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

//Define Frictional Law (function of sliprate and statevar)
double mu_friction_law_2Dv7(Real sliprate, Real statevar, Real rsf_a, Real rsf_b, Real rsf_L, Real delta_o, Real f_o)
{
  double mu = 0;
  mu = rsf_a * asinh( sliprate/(2*delta_o) * exp((f_o + rsf_b * log(delta_o * statevar/rsf_L))/rsf_a) );
  //mu = f_o + rsf_a * log(sliprate/delta_o) + rsf_b * log(delta_o*statevar/rsf_L);
  return mu;
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

    //*Nodal Mass*
    ///QUAD4 Element
    //mass_martix assembly of 
    //_test[_i][_qp] * _density[_qp] * (*_du_dotdot_du)[_qp] * _phi[this->_j][_qp]
    //where _du_dotdot_du is the spatial derivative of acceleration u_dotdot
    //this turns out to be M / dt^2 
    //the following formulation is corrected and used inside moose

    Real M = _density[_qp] * len * len * 0.5;

    //Initialize Final Traction
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal

    //Initialize Shear Perturbation
    Real Ts_perturb = _Ts_perturb[_qp]; 
    //Real Ts_perturb = 0.0; 

    //Real x_coord =_q_point[_qp](0);
    //std::cout<<x_coord<<std::endl;

    //Specify Condition
    //_fe_problem.getCurrentExecuteOnFlag()==Moose::LINEAR : Before System Solve
    //_fe_problem.getCurrentExecuteOnFlag()==Moose::NONE : After System Solve
    //if (_fe_problem.getCurrentExecuteOnFlag()=="LINEAR"){
        
        /*
        Before System Solve, compute needed traction
        */ 

        // ///Disp Plus Side
        // Real alongfaultdisp_strike_plus_t = _alongfaultdisp_strike_plus_old[_qp];
        // Real alongfaultdisp_normal_plus_t = _alongfaultdisp_normal_plus_old[_qp];
        // ///Disp Minus Side
        // Real alongfaultdisp_strike_minus_t = _alongfaultdisp_strike_minus_old[_qp];
        // Real alongfaultdisp_normal_minus_t = _alongfaultdisp_normal_minus_old[_qp];

        // // // // ///Older Value
        // Real alongfaultdisp_strike_plus_tminust = _alongfaultdisp_strike_plus_older[_qp];
        // Real alongfaultdisp_normal_plus_tminust = _alongfaultdisp_normal_plus_older[_qp];
        // Real alongfaultdisp_strike_minus_tminust = _alongfaultdisp_strike_minus_older[_qp];
        // Real alongfaultdisp_normal_minus_tminust = _alongfaultdisp_normal_minus_older[_qp];

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

        //*Solve Nonlinear Equation for sliprate at new step t+dt/2*
        //guess - unknown sliprate of next time step
        
        ///Compute known values
        Real c = len * _dt * ( M + M ) / ( M * M );
        
        //Setup outer while loop
        Real err = 1;
        Real iter = 1;
        Real v_h_pre = sliprate_mag_tminusdtover2; //sliprate value at previous time step
        Real theta_pre = statevar_t; //state variable at previous time step
        Real dv_pre = abs(sliprate_mag_tminusdtover2); //sliprate predictor initialization 
        Real dv_corr;  //sliprate corrector initialization 
        //
        Real theta_ref; 
        Real solution; 
        while ( err > 1e-5 && iter < 50){
            
            //
            dv_pre = 0.5 * ( abs(sliprate_mag_tminusdtover2) + dv_pre ); //average slip rate predictor
            //
            theta_ref = theta_pre * exp(-dv_pre*_dt/rsf_L) + (rsf_L/dv_pre) * (1-exp(-dv_pre*_dt/rsf_L)); //update state variable
            //
            Real iterr = 1;
            Real max_iter = 10000;
            Real er = 1;
            //Setup inner while loop
            Real guess_i = sliprate_mag_tminusdtover2;
            Real residual;
            Real jacobian;
            Real guess_j;
            while ( er > 1e-8 && iterr < max_iter ){  
                
                //Compute Residual
                residual = guess_i + c * Tn * mu_friction_law_2Dv7(0.5*(guess_i+sliprate_mag_tminusdtover2), theta_ref, rsf_a, rsf_b, rsf_L, delta_o, f_o) - c * Ts_trial;

                //Compute Jacobian
                Real ratioup = rsf_a * exp((f_o+rsf_b*log((delta_o*theta_ref)/(rsf_L)))/(rsf_a));
                Real ratiodown = 2 * delta_o * sqrt((exp((2*f_o+2*rsf_b*log((delta_o*theta_ref)/(rsf_L)))/(rsf_a))*(0.5*(guess_i+sliprate_strike_tminusdtover2))*(0.5*(guess_i+sliprate_strike_tminusdtover2)))/(4*delta_o*delta_o)+1);
                jacobian = 1.0 + c * Tn * ratioup / ratiodown;

                //Compute New guess
                guess_j = guess_i - residual / jacobian;

                //save
                solution = guess_j;

                //Compute err
                er = abs(guess_j - guess_i)/abs(guess_j);

                //Update Old guess
                guess_i = guess_j;

                //printout
                if (iterr == max_iter){
                    std::cout<<"NOT CONVERGED!"<<std::endl;
                }

                //
                iterr = iterr + 1;
            
            }
            //
            dv_corr = 0.5 * ( abs(sliprate_mag_tminusdtover2) + abs(guess_j)); //update slip rate corrector
            // 
            err = abs(guess_j - v_h_pre)/abs(guess_j);
            //
            v_h_pre = guess_j; //pass previous value <- current value
            //
            dv_pre = dv_corr; //pass predictor <- corrector
            //
            iter = iter + 1;

            if (iterr == max_iter){
                std::cout<<"NOT CONVERGED!"<<std::endl;
            }

        }

        Real sliprate_mag_tplusdtover2 = solution; //2D

        Real statevar_tplusdt = theta_ref;

        //*Compute shear traction at time t*
        ///trapezoidal method
        Real mu_predict = mu_friction_law_2Dv7(0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2), statevar_tplusdt, rsf_a, rsf_b, rsf_L, delta_o, f_o);
        Real T_mag = Tn * mu_predict;

        //Get traction component for t
        Ts = T_mag;

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

    //}
    //else{

        /*
        After System Solve, do nothing"
        */

        //use the latest value
        // Ts = _traction_strike[_qp];
        // Tn = _traction_normal[_qp];

    //}

    //Feed traction into the system
    //Assign back traction in CZM
    RealVectorValue traction;

    traction(0) = 0.0; 
    traction(1) = -Ts+Ts_o+Ts_perturb; 
    traction(2) = 0.0;

    _interface_traction[_qp] = traction;
    _dinterface_traction_djump[_qp] = 0;  

}