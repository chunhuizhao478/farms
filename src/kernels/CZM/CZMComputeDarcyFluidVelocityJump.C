
#include "CZMComputeDarcyFluidVelocityJump.h"

registerMooseObject("farmsApp", CZMComputeDarcyFluidVelocityJump);

InputParameters
CZMComputeDarcyFluidVelocityJump::validParams()
{
  InputParameters params = CZMComputeDarcyFluidVelocityJumpBase::validParams();
  params.addClassDescription("Compute the total fluid velocity jump across a czm interface in local "
                             "coordinates ");

  return params;
}

CZMComputeDarcyFluidVelocityJump::CZMComputeDarcyFluidVelocityJump(
    const InputParameters & parameters)
  : CZMComputeDarcyFluidVelocityJumpBase(parameters)
{
}

void
CZMComputeDarcyFluidVelocityJump::computeLocalFluidVelocityJump()
{
  _interface_fluid_vel_jump[_qp] =
      _czm_total_rotation[_qp].transpose() * _fluid_vel_jump_global[_qp];
}

