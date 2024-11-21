

#pragma once

#include "CZMComputeDarcyFluidVelocityJumpBase.h"

/**
 * Compute the interface displacement jump across a cohesive zone under the small strain
 * assumption
 */

class CZMComputeDarcyFluidVelocityJump : public CZMComputeDarcyFluidVelocityJumpBase
{
public:
  static InputParameters validParams();
  CZMComputeDarcyFluidVelocityJump(const InputParameters & parameters);

protected:
  /// compute the total displacement jump in interface coordinates
  void computeLocalFluidVelocityJump() override;

  usingCZMComputeDarcyFluidVelocityJumpBaseMembers;
};

