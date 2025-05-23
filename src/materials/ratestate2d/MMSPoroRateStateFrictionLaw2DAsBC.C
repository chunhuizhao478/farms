//Implementation of Rate-and-State Friction Law in 2D

#include "MMSPoroRateStateFrictionLaw2DAsBC.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", MMSPoroRateStateFrictionLaw2DAsBC);

InputParameters
MMSPoroRateStateFrictionLaw2DAsBC::validParams()
{ 
  //Tn_o,Ts_o,Vini,statevarini are defined in "CZMComputeLocalTractionTotalBaseRSF2D"
  InputParameters params = PoroCZMComputeLocalTractionTotalBaseRSF2D::validParams();
  params.addClassDescription("Rate-and-State Frictional Law.");
  params.addRequiredParam<Real>("len","element length");
  params.addRequiredParam<Real>("f_o","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_a","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_b","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_L","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("delta_o","slip rate parameter");
  params.addParam<Real>("elem_size", 1.0, "Value of element size");
  params.addRequiredCoupledVar("nodal_area","nodal area");
  params.addRequiredCoupledVar("reaction_rsf_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_rsf_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_rsf_pressure_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_rsf_pressure_y","reaction in y dir");
  params.addRequiredCoupledVar("interface_pressure","Pressure at sides of the fault");
  params.addRequiredCoupledVar("reaction_damp_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_damp_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_pressdamp_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_pressdamp_y","reaction in y dir");
  params.addRequiredCoupledVar("fluid_disp_x","fluid displacement in x dir");
  params.addRequiredCoupledVar("fluid_disp_y","fluid displacement in y dir");
  params.addRequiredCoupledVar("fluid_vel_x","fluid velocity in x dir");
  params.addRequiredCoupledVar("fluid_vel_y","fluid velocity in y dir");
  params.addRequiredCoupledVar("mms_sliprate","fluid velocity in x dir");
  params.addRequiredCoupledVar("mms_shear","fluid velocity in y dir");
  params.addRequiredCoupledVar("mms_normal","fluid velocity in x dir");
  params.addRequiredCoupledVar("mms_pressure","fluid velocity in y dir");
  return params;
}

MMSPoroRateStateFrictionLaw2DAsBC::MMSPoroRateStateFrictionLaw2DAsBC(const InputParameters & parameters)
  : PoroCZMComputeLocalTractionTotalBaseRSF2D(parameters),
    _len(getParam<Real>("len")),
    _f_o(getParam<Real>("f_o")),
    _rsf_a(getParam<Real>("rsf_a")),
    _rsf_b(getParam<Real>("rsf_b")),
    _rsf_L(getParam<Real>("rsf_L")),
    _delta_o(getParam<Real>("delta_o")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rhof(getMaterialPropertyByName<Real>(_base_name + "rhof")),
    _len2(coupledValue("nodal_area")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _reaction_rsf_x(coupledValue("reaction_rsf_x")),
    _reaction_rsf_y(coupledValue("reaction_rsf_y")),
    _reaction_rsf_neighbor_x(coupledNeighborValue("reaction_rsf_x")),
    _reaction_rsf_neighbor_y(coupledNeighborValue("reaction_rsf_y")),
    _reaction_rsf_pressure_x(coupledValue("reaction_rsf_pressure_x")),
    _reaction_rsf_pressure_y(coupledValue("reaction_rsf_pressure_y")),
    _reaction_rsf_neighbor_pressure_x(coupledNeighborValue("reaction_rsf_pressure_x")),
    _reaction_rsf_neighbor_pressure_y(coupledNeighborValue("reaction_rsf_pressure_y")),
    _reaction_damp_x(coupledValue("reaction_damp_x")),
    _reaction_damp_y(coupledValue("reaction_damp_y")),
    _reaction_damp_neighbor_x(coupledNeighborValue("reaction_damp_x")),
    _reaction_damp_neighbor_y(coupledNeighborValue("reaction_damp_y")),
    _reaction_pressdamp_x(coupledValue("reaction_pressdamp_x")),
    _reaction_pressdamp_y(coupledValue("reaction_pressdamp_y")),
    _reaction_pressdamp_neighbor_x(coupledNeighborValue("reaction_pressdamp_x")),
    _reaction_pressdamp_neighbor_y(coupledNeighborValue("reaction_pressdamp_y")),
    _interface_pressure_plus(coupledNeighborValue("interface_pressure")),
    _interface_pressure_minus(coupledValue("interface_pressure")),
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
    _traction_normal_old(getMaterialPropertyOldByName<Real>("traction_normal")),
    _mms_sliprate(coupledValue("mms_sliprate")),
    _mms_shear(coupledValue("mms_shear")),
    _mms_normal(coupledValue("mms_normal")),
    _mms_pressure(coupledValue("mms_pressure")),
    _fluid_disp_x(coupledValue("fluid_disp_x")),
    _fluid_disp_neighbor_x(coupledNeighborValue("fluid_disp_x")),
    _fluid_disp_y(coupledValue("fluid_disp_y")),
    _fluid_disp_neighbor_y(coupledNeighborValue("fluid_disp_y")),
    _fluid_vel_x(coupledValue("fluid_vel_x")),
    _fluid_vel_neighbor_x(coupledNeighborValue("fluid_vel_x")),
    _fluid_vel_y(coupledValue("fluid_vel_y")),
    _fluid_vel_neighbor_y(coupledNeighborValue("fluid_vel_y")),
    _traction_strike_TSN(declarePropertyByName<Real>("traction_strike_TSN")),
    _traction_normal_TSN(declarePropertyByName<Real>("traction_normal_TSN"))
{
}

void
MMSPoroRateStateFrictionLaw2DAsBC::computeInterfaceTractionAndDerivatives()
{   
    //Define Parameters
    Real len = _len2[_qp];
    Real f_o = _f_o;
    Real rsf_a = _rsf_a;
    Real rsf_b = _rsf_b;
    Real rsf_L = _rsf_L;
    Real delta_o = _delta_o;
    Real Tn_o = _Tn_o;
    Real Ts_o = _Ts_o;
    Real Ts = 0.0; //strike
    Real Tn = 0.0; //normal

    Real len_2;
    if (_current_elem->default_order() == SECOND)  // QUAD9
        len_2 = _len/6.0;
    else  // QUAD4
        len_2 = _len/2.0;

    //pressure state 
    Real p_plus = _interface_pressure_plus[_qp];
    Real p_minus = _interface_pressure_minus[_qp];
    Real p = std::max(_interface_pressure_plus[_qp], _interface_pressure_minus[_qp]);
    
    //*Restoration Force*

    //Effective stress = Stress Divergence - pressure divergence Components (label as stsdivcomp)

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_prsdivcomp(-_reaction_rsf_pressure_x[_qp],-_reaction_rsf_pressure_y[_qp], 0.0);
    RealVectorValue R_minus_global_prsdivcomp(-_reaction_rsf_neighbor_pressure_x[_qp],-_reaction_rsf_neighbor_pressure_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_prsdivcomp = _rot[_qp].transpose() * R_plus_global_prsdivcomp;
    RealVectorValue R_minus_local_prsdivcomp = _rot[_qp].transpose() * R_minus_global_prsdivcomp;
    
    //*Restoration Force*

    //Stress Divergence Components (label as stsdivcomp)

    //--------------------------------------------------------------------------------------------------//

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
    Real R_plus_local_normal_stsdivcomp  = R_plus_local_stsdivcomp(0) + R_plus_local_prsdivcomp(0);
    Real R_plus_local_strike_stsdivcomp  = R_plus_local_stsdivcomp(1) + R_plus_local_prsdivcomp(1) ;
    
    Real R_minus_local_normal_stsdivcomp = R_minus_local_stsdivcomp(0) + R_minus_local_prsdivcomp(0) ;
    Real R_minus_local_strike_stsdivcomp = R_minus_local_stsdivcomp(1) + R_minus_local_prsdivcomp(1) ;

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

    ///Define in global coordinate
    //current time step 
    RealVectorValue R_plus_global_pressdampingcomp(-_reaction_pressdamp_x[_qp],-_reaction_pressdamp_y[_qp], 0.0);
    RealVectorValue R_minus_global_pressdampingcomp(-_reaction_pressdamp_neighbor_x[_qp],-_reaction_pressdamp_neighbor_y[_qp], 0.0);

    ///Rotate in local coordinate
    //current time step
    RealVectorValue R_plus_local_pressdampingcomp = _rot[_qp].transpose() * R_plus_global_pressdampingcomp;
    RealVectorValue R_minus_local_pressdampingcomp = _rot[_qp].transpose() * R_minus_global_pressdampingcomp;

    ///Get Components
    //current time step
    Real R_plus_local_normal_dampingcomp  = R_plus_local_dampingcomp(0) ;
    Real R_plus_local_strike_dampingcomp  = R_plus_local_dampingcomp(1) ;
    
    Real R_minus_local_normal_dampingcomp = R_minus_local_dampingcomp(0) ;
    Real R_minus_local_strike_dampingcomp = R_minus_local_dampingcomp(1) ;

    //--------------------------------------------------------------------------------------------------//

    //Add restoration forces from two contributions
    Real R_plus_local_normal  = R_plus_local_normal_stsdivcomp  + R_plus_local_normal_dampingcomp;
    Real R_plus_local_strike  = R_plus_local_strike_stsdivcomp  + R_plus_local_strike_dampingcomp;
    Real R_minus_local_normal = R_minus_local_normal_stsdivcomp + R_minus_local_normal_dampingcomp;
    Real R_minus_local_strike = R_minus_local_strike_stsdivcomp + R_minus_local_strike_dampingcomp;  

    //*Nodal Mass*
    ///QUAD4 Element
    Real M = _density[_qp] * len * len_2;
    Real Mf = _rhof[_qp] * len * len_2;

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

    //Global Fluid Displacement Jump 
    RealVectorValue fluid_disp_jump_global = (_fluid_disp_x[_qp]-_fluid_disp_neighbor_x[_qp],_fluid_disp_y[_qp]-_fluid_disp_neighbor_y[_qp]);
  
    //Global Fluid velocity Jump 
    RealVectorValue fluid_vel_jump_global = (_fluid_vel_x[_qp]-_fluid_vel_neighbor_x[_qp],_fluid_vel_y[_qp]-_fluid_vel_neighbor_y[_qp]);;

    //Local Fluid Displacement Jump 
    RealVectorValue fluid_disp_jump = _rot[_qp].transpose() * fluid_disp_jump_global;
    RealVectorValue fluid_vel_jump = _rot[_qp].transpose() * fluid_vel_jump_global;
  
    //*Get slip rate and slip from time t-dt/2 and t
    ///strike direction
    Real sliprate_strike_tminusdtover2 = _sliprate_strike_old[_qp];
    Real slip_strike_t = _slip_strike_old[_qp];

    ///normal direction
    Real sliprate_normal_tminusdtover2 = _sliprate_normal_old[_qp];
    Real slip_normal_t = _slip_normal_old[_qp];

    ///mag
    Real sliprate_mag_tminusdtover2 = _sliprate_mag_old[_qp];

    // Compute outputs that are meaningful at any y value
    Real shear_stress_mms = _mms_shear[_qp];
    Real normal_stress_mms = _mms_normal[_qp];
    Real source_term = _mms_pressure[_qp];
    Real slip_rate_mms = _mms_sliprate[_qp];

    //*Compute Trial Normal Traction and Normal Traction*
    Real Tn_trial = ( -1.0 * (1.0/_dt) * M * M * ( sliprate_normal_tminusdtover2 + (1.0/_dt) * slip_normal_t) ) / ( len * (M + M) ) 
                  + ( -1.0 * (1.0/_dt) * Mf * Mf * ( fluid_vel_jump(0) + (1.0/_dt) * fluid_disp_jump(0)) ) / ( len * (Mf + Mf) ) 
                  + ( M * R_minus_local_normal - M * R_plus_local_normal ) / ( len * (M + M) ) - Tn_o;

    if (Tn_trial<0)
    {
        Tn = Tn_trial;
    }else{
        Tn = 0;
    }

    // ///Make Tn positive
    // Tn = abs(Tn) - p;

    //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
    Real Ts_trial = ( M * M * sliprate_strike_tminusdtover2 )/( len * _dt * (M + M) ) 
                  + ( Mf * Mf * fluid_vel_jump(1) )/( len * _dt * (Mf + Mf) ) 
                  + (M * R_plus_local_strike - M * R_minus_local_strike) / ( len * ( M + M ) ) + Ts_o;

    // Moose::out << "fluid_vel_jump: " <<  fluid_vel_jump(0)  << std::endl;   
    // Moose::out << "fluid_disp_jump: " <<  fluid_disp_jump(0)<< std::endl;   

    // Moose::out << "shear_stress_mms: " << shear_stress_mms << std::endl; 
              
    Real Tmag_trial = sqrt(Ts_trial*Ts_trial);

    //const
    Real c = len * _dt * ( M + M ) / (M * M);
    //Get State Variable at Current Time Step
    Real statevar_t = _statevar_old[_qp];


    // // Calculate state variable based on rate-and-state friction law
    // // This should be inserted into your computeInterfaceValues function

    // // Add a small epsilon value to prevent division by zero
    // Real f_mms = (shear_stress_mms+Ts_o) / (normal_stress_mms + Tn_o - max_pressure_mms);

    // // Use the rate-and-state parameters from your class
    // Real V_mms = slip_rate_mms + 1e-12;          // Slip rate from MMS solution
    // Real V0 = _delta_o;                  // Reference velocity 
    // Real DRS = _rsf_L;                   // Characteristic slip distance
    // Real f0 = _f_o;                      // Reference friction coefficient
    // Real a_rsf = _rsf_a;                 // Direct effect parameter
    // Real b_rsf = _rsf_b;                 // Evolution effect parameter

    // // Calculate the intermediate term in the equation
    // Real sinh_term = std::sinh((f_mms - f0) / a_rsf);
    // Real inside_ln = sinh_term * (2.0 * V0 / V_mms);

    // // Calculate the exponent term
    // Real exponent = (a_rsf * std::log(inside_ln) - f0) / b_rsf;

    // // Calculate the state variable
    // Real statevar_mms = (DRS / V0) * std::exp(exponent);

    // statevar_t = source_term;
  
    Real Z = 0.5 / delta_o * exp((f_o + rsf_b * log(delta_o * (statevar_t)/rsf_L))/rsf_a);  

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
    Real statevar_tplusdt = statevar_t * coeffD + (1+source_term) * (rsf_L/sliprate_mag_tplusdtover2) * (1-coeffD);

    Moose::out << "statevar_tplusdt: " << source_term << std::endl; 


    // Real statevar_tplusdt = source_term;

    //*Compute shear traction at time t*
    Real T_mag = Tn * rsf_a * asinh( 0.5*(sliprate_mag_tminusdtover2+sliprate_mag_tplusdtover2) * Z ); 

    ///Get Components
    Ts = T_mag * ( Ts_trial / Tmag_trial );

    ////DISP
    Real du_strike_plus_t  =  alongfaultdisp_strike_plus_t  - alongfaultdisp_strike_plus_tminust  + _dt * _dt / M * (R_plus_local_strike  - len * (Ts - Ts_o));
    Real du_normal_plus_t  =  alongfaultdisp_normal_plus_t  - alongfaultdisp_normal_plus_tminust  + _dt * _dt / M * (R_plus_local_normal  - len * (Tn - Tn_o));
    
    Real du_strike_minus_t =  alongfaultdisp_strike_minus_t  - alongfaultdisp_strike_minus_tminust + _dt * _dt / M * (R_minus_local_strike  + len * (Ts - Ts_o));
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
    // _traction_strike_TSN[_qp] = Ts_trial - Ts_o ;
    // _traction_normal_TSN[_qp] =  ( -1.0 * (1.0/_dt) * M * M * (  slip_rate_mms * 0.1 / 5.2 + (1.0/_dt) * slip_rate_mms / 3 * _t * 0.1 / 5.2 ) ) / ( len * (M + M) ) 
    //               + ( -1.0 * (1.0/_dt) * Mf * Mf * ( shear_stress_mms + (1.0/_dt) * source_term) ) / ( len * (Mf + Mf) ) 
    //               + ( M * R_minus_local_normal - M * R_plus_local_normal ) / ( len * (M + M) )   ;
   
    // _traction_strike_TSN[_qp] = fluid_vel_jump(1);
    // _traction_normal_TSN[_qp] = fluid_disp_jump(1);


    //  Moose::out << "_traction_strike: " <<  _traction_strike[_qp] << std::endl;   
    //  Moose::out << "_traction_normal: " <<  _traction_normal_TSN[_qp]  << std::endl; 

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

}
  