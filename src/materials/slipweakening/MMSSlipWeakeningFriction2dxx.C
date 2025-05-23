#include "MMSSlipWeakeningFriction2dxx.h"
#include "InterfaceKernel.h"

// User-Specific App name
registerMooseObject("farmsApp", MMSSlipWeakeningFriction2dxx);

InputParameters
MMSSlipWeakeningFriction2dxx::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("MMS Slip Weakening Traction Separation Law with exponential time function.");
  params.addParam<Real>("T2_o", 50e6, "Background normal traction");
  params.addRequiredCoupledVar("nodal_area", "nodal area");
  params.addParam<Real>("mu_d", 0.7, "Value of dynamic friction parameter");
  params.addParam<Real>("Dc", 0.75, "Value of characteristic length");
  params.addParam<Real>("delta", 1.0, "Total slip amplitude (m)");
  params.addParam<Real>("R_m", 1000.0, "Characteristic length (m)");
  params.addParam<Real>("tbar", 0.25, "Reference time (s)");
  params.addParam<Real>("tw", 0.1, "Time scale parameter (s)");
  params.addParam<Real>("Vmin", 1e-9, "Minimum slip velocity (m/s)");
  params.addParam<Real>("tau_excess", 20e6, "Excess stress for nucleation (Pa)");
  params.addParam<Real>("R_nuc", 500.0, "Nucleation radius (m)");
  params.addParam<Real>("tau_background", 40e6, "Background stress (Pa)");
  params.addParam<bool>("exact_traction", false, "Whether to use exact analytical traction");
  params.addParam<Real>("lambda", 20e9, "First Lamé parameter");
  params.addParam<Real>("shear_modulus", 30e9, "Shear modulus (second Lamé parameter)");
  return params;
}

MMSSlipWeakeningFriction2dxx::MMSSlipWeakeningFriction2dxx(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
<<<<<<< HEAD
    _nodal_area(coupledValue("nodal_area")),
    _T2_o(getParam<Real>("T2_o")), 
=======
    _T2_o(getParam<Real>("T2_o")), 
    _nodal_area(coupledValue("nodal_area")),
>>>>>>> 6188d19945ce13b4debb058d2e709731e3f73bb8
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _delta(getParam<Real>("delta")),
    _R_m(getParam<Real>("R_m")),
    _tbar(getParam<Real>("tbar")),
    _tw(getParam<Real>("tw")),
    _Vmin(getParam<Real>("Vmin")),
    _tau_excess(getParam<Real>("tau_excess")),
    _R_nuc(getParam<Real>("R_nuc")),
    _tau_background(getParam<Real>("tau_background")),
    _exact_traction(getParam<bool>("exact_traction")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _interface_displacement_jump_old(getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_older(getMaterialPropertyOlder<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _lambda(getParam<Real>("lambda")),
    _mu(getParam<Real>("shear_modulus")),
    _width(_R_m / 10.0) // Based on the value 100.0 from the ParsedFunction
{
}

// Spatial distribution function phi(x,y) for MMS
Real
MMSSlipWeakeningFriction2dxx::spatialFunction(Real x, Real y) const
{
  // Use exponential function for spatial distribution
  return std::exp(-(x*x + y*y)/(2.0*_width*_width));
}

// Temporal function K(t) for MMS
Real
MMSSlipWeakeningFriction2dxx::temporalFunction(Real t) const
{
  // Based on the provided ParsedFunction: exp((t - tbar)/tw)
  return std::exp((t - _tbar)/_tw);
}

// Temporal derivative of K(t) for MMS
Real
MMSSlipWeakeningFriction2dxx::temporalDerivative(Real t) const
{
  // Derivative of exp((t - tbar)/tw) is (1/tw)*exp((t - tbar)/tw)
  return (1.0/_tw) * temporalFunction(t);
}

// Initial stress distribution for MMS
Real
MMSSlipWeakeningFriction2dxx::initialStress(Real x) const
{
  return _tau_background;
}

// Exact traction based on MMS (analytical formula)
Real
MMSSlipWeakeningFriction2dxx::exactTraction(Real x, Real y, Real t) const
{
  Real analytical_slip = _delta * temporalFunction(t) * spatialFunction(x, y);
  
  // Apply slip-weakening law
  Real mu_s = 1.0; // tau_s / T2_o = 50e6 / 50e6 = 1.0
  Real mu_d = _mu_d;
  Real Dc = _Dc;
  Real T2 = -_T2_o;
  
  Real tau_f;
  if (analytical_slip < Dc) {
    tau_f = (mu_s - (mu_s - mu_d) * analytical_slip / Dc) * (-T2);
  } else {
    tau_f = mu_d * (-T2);
  }
  
  return tau_f;
}


void
MMSSlipWeakeningFriction2dxx::computeInterfaceTractionAndDerivatives()
{   
  // Obtain coordinates of current qp and time
  Real x_coord = _q_point[_qp](0);
  Real y_coord = _q_point[_qp](1);
  Real t = _t;
  
  // Calculate analytical slip and slip rate from MMS solution
  Real analytical_slip = _delta * temporalFunction(t) * spatialFunction(x_coord, y_coord);
  Real analytical_slip_rate = _delta * temporalDerivative(t) * spatialFunction(x_coord, y_coord);
  
  // Determine background shear stress
  Real T1_o = 0.0;
  if (t < 1e-6) {
    T1_o = initialStress(x_coord);
  } else {
    T1_o = _tau_background; // Use background stress after nucleation
  }
  
  // Get the current stress state from the material
  RealVectorValue stress_y = _stress[_qp].row(1);  
  RealVectorValue dstress_y = _dstress[_qp].row(1);  
  Real shear = stress_y(0); 
  Real dshear = dstress_y(0); 
  Real normal = stress_y(1); // Normal stress component
  
  
  // Calculate the derivatives of the spatial function
  Real dphi_dx = -x_coord/(_width*_width) * spatialFunction(x_coord, y_coord);
  Real dphi_dy = -y_coord/(_width*_width) * spatialFunction(x_coord, y_coord);
    
  // Contribution to normal stress from the analytical solution
  Real sigma_xy = _mu * _delta/2 * temporalFunction(t) * (-y_coord/(_width*_width) * spatialFunction(x_coord, y_coord));

  // Calculate the trial traction using analytical slip rate
  Real area = _nodal_area[_qp];
  Real M = _density[_qp] * area * area/2;
  
  // Use analytical slip rate instead of numerical one
  Real T_trial = shear + T1_o -(1/_dt)*M*analytical_slip_rate/(2*area);
  
  // Normal stress (always negative in compression)
  Real T2 = -_T2_o;
  
  // Slip-weakening friction parameters
  Real mu_s = 1.0; // Static friction coefficient (tau_s / T2_o)
  Real mu_d = _mu_d; // Dynamic friction coefficient
  Real Dc = _Dc; // Critical slip distance
  
  // Calculate friction coefficient based on analytical slip
  Real mu;
  if (analytical_slip < Dc) {
    mu = mu_s - (mu_s - mu_d) * analytical_slip / Dc;
  } else {
    mu = mu_d;
  }
  
  // Calculate strength (maximum allowed friction)
  Real tau_strength = mu * (-T2);
  
  // Compute the actual traction based on trial value and strength
  Real T1;
  Real dT1d_disp;
  
  if (std::abs(T_trial) <= tau_strength) {
    // Stick condition - use trial traction
    T1 = T_trial;
    dT1d_disp = -(1/_dt)*M*(1/_dt)/(2*area) + dshear;
  } else {
    // Slip condition - traction limited by strength
    T1 = tau_strength * std::copysign(1.0, T_trial);
    
    // Calculate derivative for slip condition
    Real dtau_f = 0.0;
    if (analytical_slip < Dc) {
      dtau_f = -(mu_s - mu_d) / Dc * (-T2);
    }
    
    dT1d_disp = dtau_f;
  }

  // Calculate the source term needed to satisfy the slip-weakening friction law
  Real source = sigma_xy - (T1-T1_o);
  
  // Assign traction as BC in CZM
  RealVectorValue traction;
  traction(0) = 0; // No opening component
  traction(1) = 0; // Apply calculated traction plus source term
  traction(2) = 0; // No out-of-plane component
  
  // Set traction and derivative for solver
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = -dT1d_disp;
}