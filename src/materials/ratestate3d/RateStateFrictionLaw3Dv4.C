//Implementation of Rate-and-State Friction Law in 3D

#include "RateStateFrictionLaw3Dv4.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", RateStateFrictionLaw3Dv4);

InputParameters
RateStateFrictionLaw3Dv4::validParams()
{ 
  //Tn_o,Ts_o,Vini,statevarini are defined in "CZMComputeLocalTractionTotalBaseRSF3D"
  InputParameters params = CZMComputeLocalTractionTotalBaseRSF3D::validParams();
  params.addClassDescription("Rate-and-State Frictional Law.");
  params.addParam<Real>("len",200,"element length");
  params.addParam<Real>("f_o",0.6,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_a",0.008,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_b",0.012,"rate-and-state friction coefficients");
  params.addParam<Real>("rsf_L",0.02,"rate-and-state friction coefficients");
  params.addParam<Real>("delta_o",1e-6,"slip rate parameter");
  params.addRequiredCoupledVar("reaction_rsf_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_rsf_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_rsf_z","reaction in z dir");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir");
  return params;
}

RateStateFrictionLaw3Dv4::RateStateFrictionLaw3Dv4(const InputParameters & parameters)
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
    _traction_strike_old(getMaterialPropertyOldByName<Real>("traction_strike")),
    _traction_normal_old(getMaterialPropertyOldByName<Real>("traction_normal")),
    _traction_dip_old(getMaterialPropertyOldByName<Real>("traction_dip"))
{
}

//Define Frictional Law (function of sliprate and statevar)
double mu_friction_law_3Dv4(Real sliprate, Real statevar, Real rsf_a, Real rsf_b, Real rsf_L, Real delta_o, Real f_o)
{
  double mu = 0;
  mu = rsf_a * asinh( sliprate/(2*delta_o) * exp((f_o + rsf_b * log(delta_o * statevar/rsf_L))/rsf_a) );
  return mu;
}

void
RateStateFrictionLaw3Dv4::computeInterfaceTractionAndDerivatives()
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
    
    //*Restoration Force*
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
    Real R_plus_local_strike  = R_plus_local(1);
    Real R_plus_local_normal  = R_plus_local(0);
    Real R_plus_local_dip     = R_plus_local(2);

    Real R_minus_local_strike = R_minus_local(1);
    Real R_minus_local_normal = R_minus_local(0);
    Real R_minus_local_dip    = R_minus_local(2);

    //*Nodal Mass*
    ///HEX8 Element
    //see "ActuallyExplicitEuler::solve()"
    //mass-matrix (len=200) = 1.068e+10 = 2670 * 200 * 200 * 100
    //the following equation is correct!
    Real M = _density[_qp] * len * len * len * 0.5;

    //Initialize Final Traction
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal
    Real Td = 0.0; //dip

    //Initialize shear stress perturbation
    Real Ts_perturb = 0.0;
    
    //By testing, the area used in CZM Residual is (half of element edge length)^2,
    //considering 4 quarture point in a quad4 x-z surface element,
    //it turns out to be the same length

    Real len_len = ( 1.0 * len ) * ( 1.0 * len );

    //Specify Condition
    //_fe_problem.getCurrentExecuteOnFlag()==Moose::NONE : After System Solve
    //_fe_problem.getCurrentExecuteOnFlag()==Moose::LINEAR : Before System Solve
    // if (_fe_problem.getCurrentExecuteOnFlag()=="LINEAR"){
        
        /*
        Before System Solve, quantities are allowed to be updated"
        */
        
        //Assign shear perturbation at t, t-dt
        Ts_perturb = _Ts_perturb[_qp];
    
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
        //
        Real c = len_len * _dt * ( M + M ) / (M * M);
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
                residual = guess_i + c * Tn * mu_friction_law_3Dv4(0.5*(guess_i+sliprate_mag_tminusdtover2), theta_ref, rsf_a, rsf_b, rsf_L, delta_o, f_o) - c * Tmag_trial;

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
                // if (x_coord > 0 && x_coord < 30){
                //     std::cout<<"inner er: "<< er <<std::endl;
                // }
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

            // if (x_coord > 0 && x_coord < 30){
            //     std::cout<<"outer err: "<< er <<std::endl;
            // }

        }

        Real sliprate_mag_tplusdtover2 = solution; 

        Real statevar_tplusdt = theta_ref;

        //*Compute shear traction at time t*
        ///trapezoidal method
        Real mu_predict = mu_friction_law_3Dv4(0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2), theta_ref, rsf_a, rsf_b, rsf_L, delta_o, f_o);
        Real T_mag = Tn * mu_predict;

        ///Get Components
        Ts = T_mag * ( Ts_trial / Tmag_trial );
        Td = T_mag * ( Td_trial / Tmag_trial );

        //*Compute quantities at time t + dt/2 or t + dt
        // _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual(type); : _JxW[_qp] = 10000 (len = 200)
        // From CentralDifference.C: _sys.hasVector(_u_dotdot_factor_tag)
        // expand u_dotdot := (u - 2 * u_old + u_old_old)/dt^2 = (u - u_old)/dt^2 + (u_old_old - u_old)/dt^2
        //                                                     = (du)/dt^2        + _u_dotdot_factor_tag
        // the tag goes into "InertialForceTempl<true>::computeResidualAdditional()" for "explicit" and "lumped" only
        //     
        //     this->_local_re(i) *= (*_u_dotdot_factor_dof)[i] + _eta[0] * (*_u_dot_factor_dof)[i]; 
        //
        // where _local_re(i) is density[_qp] * _test[_qp]
        // this gives M * (u_old_old - u_old)/dt^2
        //
        // by assembling "mass_matrix(M/dt^2)" (see ActuallyExplicitEuler::solve()), the lhs is M * (du)/dt^2, where du is displacement increment
        //
        // the "computeResidualAdditional()" function is evaluted in "TimeKernel::computeResidual()":

        // see ActuallyExplicitEuler::solve() -> 
        // //Perform the linear solve
        // bool converged = performExplicitSolve(mass_matrix);

        ////DISP
        // Real du_strike_plus_t  =  alongfaultdisp_strike_plus_t  - alongfaultdisp_strike_plus_tminust  + _dt * _dt / M * (R_plus_local_strike  - len_len * (Ts - Ts_o - Ts_perturb));
        // Real du_normal_plus_t  =  alongfaultdisp_normal_plus_t  - alongfaultdisp_normal_plus_tminust  + _dt * _dt / M * (R_plus_local_normal  - len_len * (Tn - Tn_o));
        // Real du_dip_plus_t     =  alongfaultdisp_dip_plus_t     - alongfaultdisp_dip_plus_tminust     + _dt * _dt / M * (R_plus_local_dip     - len_len * (Td - Td_o));
        
        // Real du_strike_minus_t =  alongfaultdisp_strike_minus_t  - alongfaultdisp_strike_minus_tminust + _dt * _dt / M * (R_minus_local_strike + len_len * (Ts - Ts_o - Ts_perturb));
        // Real du_normal_minus_t =  alongfaultdisp_normal_minus_t  - alongfaultdisp_normal_minus_tminust + _dt * _dt / M * (R_minus_local_normal + len_len * (Tn - Tn_o));
        // Real du_dip_minus_t    =  alongfaultdisp_dip_minus_t     - alongfaultdisp_dip_minus_tminust    + _dt * _dt / M * (R_minus_local_dip    + len_len * (Td - Td_o));

        // std::cout<<"material object: "<<(Ts - Ts_o - Ts_perturb)<<std::endl;

        //see ActuallyExplicitEuler::solve()
        //// Update the solution
        //*_nonlinear_implicit_system->solution = _nl.solutionOld();
        //*_nonlinear_implicit_system->solution += _solution_update;

        // Real alongfaultdisp_strike_plus_tplusdt  = alongfaultdisp_strike_plus_t  + du_strike_plus_t;
        // Real alongfaultdisp_normal_plus_tplusdt  = alongfaultdisp_normal_plus_t  + du_normal_plus_t;
        // Real alongfaultdisp_dip_plus_tplusdt     = alongfaultdisp_dip_plus_t     + du_dip_plus_t;
        
        // Real alongfaultdisp_strike_minus_tplusdt = alongfaultdisp_strike_minus_t + du_strike_minus_t;
        // Real alongfaultdisp_normal_minus_tplusdt = alongfaultdisp_normal_minus_t + du_normal_minus_t;
        // Real alongfaultdisp_dip_minus_tplusdt    = alongfaultdisp_dip_minus_t    + du_dip_minus_t;

        ////VEL
        // Real alongfaultvel_strike_plus_tplusdtover2  = (alongfaultdisp_strike_plus_tplusdt  - alongfaultdisp_strike_plus_t ) / ( 1 * _dt );
        // Real alongfaultvel_normal_plus_tplusdtover2  = (alongfaultdisp_normal_plus_tplusdt  - alongfaultdisp_normal_plus_t ) / ( 1 * _dt );
        // Real alongfaultvel_dip_plus_tplusdtover2     = (alongfaultdisp_dip_plus_tplusdt     - alongfaultdisp_dip_plus_t    ) / ( 1 * _dt );
        
        // Real alongfaultvel_strike_minus_tplusdtover2 = (alongfaultdisp_strike_minus_tplusdt - alongfaultdisp_strike_minus_t) / ( 1 * _dt );
        // Real alongfaultvel_normal_minus_tplusdtover2 = (alongfaultdisp_normal_minus_tplusdt - alongfaultdisp_normal_minus_t) / ( 1 * _dt );
        // Real alongfaultvel_dip_minus_tplusdtover2    = (alongfaultdisp_dip_minus_tplusdt    - alongfaultdisp_dip_minus_t   ) / ( 1 * _dt );

        //Update Slip Rate and Slip at t + dt/2 or t + dt
        ///Slip Rate
        Real sliprate_strike_tplusdtover2 = sliprate_mag_tplusdtover2 * ( Ts_trial / Tmag_trial );
        Real sliprate_normal_tplusdtover2 = 0.0;
        Real sliprate_dip_tplusdtover2    = sliprate_mag_tplusdtover2 * ( Td_trial / Tmag_trial );

        ///Slip
        Real slip_strike_tplusdtover2 = slip_strike_t + sliprate_strike_tplusdtover2 * _dt;
        Real slip_normal_tplusdtover2 = 0.0;
        Real slip_dip_tplusdtover2    = slip_dip_t + sliprate_dip_tplusdtover2 * _dt;

        //Save Results
        ///Traction (t)
        _traction_strike[_qp] = Ts;
        _traction_normal[_qp] = Tn;
        _traction_dip[_qp] = Td;

        // ///Vel (t+dt/2)
        // _alongfaultvel_strike_plus[_qp]  = alongfaultvel_strike_plus_tplusdtover2;
        // _alongfaultvel_normal_plus[_qp]  = alongfaultvel_normal_plus_tplusdtover2;
        // _alongfaultvel_dip_plus[_qp]     = alongfaultvel_dip_plus_tplusdtover2;

        // _alongfaultvel_strike_minus[_qp] = alongfaultvel_strike_minus_tplusdtover2;
        // _alongfaultvel_normal_minus[_qp] = alongfaultvel_normal_minus_tplusdtover2;
        // _alongfaultvel_dip_minus[_qp]    = alongfaultvel_dip_minus_tplusdtover2;

        ///Disp (t+dt)
        // _alongfaultdisp_strike_plus[_qp] = alongfaultdisp_strike_plus_tplusdt;
        // _alongfaultdisp_normal_plus[_qp] = alongfaultdisp_normal_plus_tplusdt;
        // _alongfaultdisp_dip_plus[_qp]    = alongfaultdisp_dip_plus_tplusdt;

        // _alongfaultdisp_strike_minus[_qp] = alongfaultdisp_strike_minus_tplusdt;
        // _alongfaultdisp_normal_minus[_qp] = alongfaultdisp_normal_minus_tplusdt;
        // _alongfaultdisp_dip_minus[_qp]    = alongfaultdisp_dip_minus_tplusdt;

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
        
        //Real statevar_ss = rsf_L / _sliprate_mag[_qp];
        //_statevar[_qp] = statevar_ss + (_statevar_old[_qp] - statevar_ss) * exp( (-1.0 * _sliprate_mag[_qp] * _dt) / rsf_L);

        ///State Variable (t+dt)
        _statevar[_qp] = statevar_tplusdt;

    // }
    // else{

        /*
        After System Solve, keep the traction same as last step
        */

        //get latest value
        // Ts = _traction_strike[_qp];
        // Tn = _traction_normal[_qp];
        // Td = _traction_dip[_qp];

    // }

    //Feed traction into the system
    //Assign back traction in CZM
    RealVectorValue traction;

    //see CZMInterfaceKernelBase::computeQpResidual(Moose::DGResidualType type)
    // [test_secondary-test_primary]*T where T represents the traction.
    //see ActuallyExplicitEuler::solve()
    // // Move the residual to the RHS
    // _explicit_residual *= -1.0;
    // This term becomes:
    // [test_primary-test_secondary]*T where T represents the traction.
    // then T is -(Ts - Ts_o - Ts_perturb) = -Ts+Ts_o+Ts_perturb

    traction(0) = 0.0; 
    traction(1) = -Ts+Ts_o+Ts_perturb; 
    traction(2) = -Td+Td_o;

    // traction(0) = 0.0; 
    // traction(1) = 0.0; 
    // traction(2) = 0.0;

    _interface_traction[_qp] = traction;
    _dinterface_traction_djump[_qp] = 0;   

}