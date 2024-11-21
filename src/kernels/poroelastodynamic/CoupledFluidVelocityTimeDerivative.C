

#include "CoupledFluidVelocityTimeDerivative.h"
#include "MooseVariable.h"

registerMooseObject("farmsApp", CoupledFluidVelocityTimeDerivative);

InputParameters
CoupledFluidVelocityTimeDerivative::validParams()
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

CoupledFluidVelocityTimeDerivative::CoupledFluidVelocityTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters), 
    _v_dot(coupledDot("v")), 
    _dv_dot(coupledDotDu("v")),  
    _v_var(coupled("v")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_density_qp(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp"))
{
}

Real
CoupledFluidVelocityTimeDerivative::computeQpResidual()
{

  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }
  return _test[_i][_qp] * dens * _v_dot[_qp];
  
}

Real
CoupledFluidVelocityTimeDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
CoupledFluidVelocityTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real dens = 0.0;
   for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }

  if (jvar == _v_var)

    return _test[_i][_qp]  * dens * _dv_dot[_qp] * _phi[_j][_qp];
  
  return 0.0;
}
 