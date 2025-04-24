#include "CheckAlpha.h"
#include "SystemBase.h"

registerMooseObject("farmsApp", CheckAlpha);

InputParameters
CheckAlpha::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Copies the specified variable into an auxiliary variable.");
  // We use the key "source" for the variable to be copied.
  params.addRequiredCoupledVar("source", "Variable to read from.");
  params.addRequiredCoupledVar("initial_damage_aux", "Initial damage variable");
  // Optionally you can also require another coupled var if needed,
  // but make sure the naming is consistent.
  MooseEnum stateEnum("CURRENT=0 OLD=1 OLDER=2", "CURRENT");
  params.addParam<MooseEnum>("state", stateEnum,
      "The state to copy (CURRENT, OLD, or OLDER).");
  return params;
}

CheckAlpha::CheckAlpha(const InputParameters & parameters)
  : AuxKernel(parameters),
    _state(getParam<MooseEnum>("state"))
{
  // Store pointers from the coupled variable interface
  _v = &coupledValue("source");

  _initial_damage_aux=&coupledValue("initial_damage_aux");
  
  // getVar returns a pointer, so we store it as such.
  _source_variable = getVar("source", 0);
  
  // Check that the finite element families and orders match.
  if (_var.feType().family != _source_variable->feType().family)
    paramError("source",
               "Source (" + Moose::stringify(_source_variable->feType().family) + 
               ") and target (" + Moose::stringify(_var.feType().family) +
               ") variable have different families. You may use a ProjectionAux instead of CheckAlpha");
  if (_var.feType().order != _source_variable->feType().order)
    paramError("source",
               "Source (" + Moose::stringify(_source_variable->feType().order) +
               ") and target (" + Moose::stringify(_var.feType().order) +
               ") variable are of different orders. You may use a ProjectionAux instead of CheckAlpha");
}

Real
CheckAlpha::computeValue()
{

  if ( (*_v)[_qp] < (*_initial_damage_aux)[_qp] )
  {
    // If the value is less than the initial damage, set it to the initial damage.
    return (*_initial_damage_aux)[_qp];
  }
  else
  {
    // Otherwise, set it to the current value.
    return (*_v)[_qp];
  }
}