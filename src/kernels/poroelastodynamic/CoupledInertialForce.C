

#include "CoupledInertialForce.h"
#include "MooseVariable.h"

registerMooseObject("farmsApp", CoupledInertialForce);

InputParameters
CoupledInertialForce::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("Time derivative Kernel that acts on a coupled variable. Weak form: "
                             "$(\\psi_i, \\frac{\\partial v_h}{\\partial t})$.");
  params.addRequiredCoupledVar("v", "Coupled variable");
  params.addParam<Real>("eta",
                                        0.0,
                                        "Name of material property or a constant real "
                                        "number defining the eta parameter for the "
                                        "Rayleigh damping.");
  params.addParam<unsigned int>(
      "fluid_component", 0, "The index corresponding to the component for this kernel"); 
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

CoupledInertialForce::CoupledInertialForce(const InputParameters & parameters)
  : TimeKernel(parameters), 
    _v_dot_dot(coupledDotDot("v")), 
    _v_dot(coupledDot("v")), 
    _dv_dot_dot(coupledDotDotDu("v")),   
    _dv_dot(coupledDotDu("v")),  
    _v_var(coupled("v")),
    _eta(getParam<Real>("eta")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_density_qp(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp"))
{
}

Real
CoupledInertialForce::computeQpResidual()
{

  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }
  return _test[_i][_qp] * dens * (_v_dot_dot[_qp] + _eta  *  _v_dot[_qp]);
  
}

Real
CoupledInertialForce::computeQpJacobian()
{
  return 0.0;
}

Real
CoupledInertialForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real dens = 0.0;

  if (jvar == _v_var)

   for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }
    return _test[_i][_qp]  * dens * _dv_dot_dot[_qp] * _phi[_j][_qp] ; 
            + _eta  * _test[_i][_qp] * dens * _dv_dot[_qp] * _phi[_j][_qp];
  

  return 0.0;
}
