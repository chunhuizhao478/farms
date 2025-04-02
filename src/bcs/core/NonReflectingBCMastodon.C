/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     MASTODON                  */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/
#include "Function.h"
#include "MooseError.h"
#include "MooseMesh.h"
#include "NonReflectingBCMastodon.h"

registerMooseObject("farmsApp", NonReflectingBCMastodon);

InputParameters
NonReflectingBCMastodon::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Applies Lysmer damper in the normal and "
                             "tangential directions to soil boundary.");
  params += NonReflectingBCMastodon::commonParameters();
  params.addRequiredParam<unsigned int>("component",
                                        "The direction in which the Lysmer damper is applied.");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

InputParameters
NonReflectingBCMastodon::commonParameters()
{
  InputParameters params = emptyInputParameters();
  params.addCoupledVar("displacements",
                       "The vector of displacement variables. "
                       "The size of this vector must be same "
                       "as the number of dimensions.");
  params.addRequiredRangeCheckedParam<Real>(
      "p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_wave_speed", "shear_wave_speed>0.0", "shear wave speed of the material.");
  return params;
}

NonReflectingBCMastodon::NonReflectingBCMastodon(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _disp_var(3),
    _disp_older(3),
    _density(getMaterialPropertyByName<Real>("density")),
    _p_wave_speed(getParam<Real>("p_wave_speed")),
    _shear_wave_speed(getParam<Real>("shear_wave_speed"))
{

  // Error checking on variable vectors
  if (_ndisp != _mesh.dimension())
    mooseError("The number of variables listed in the 'displacements' parameter in \"",
               name(),
               "\" block must match the mesh dimension.");

  if (_component >= _mesh.dimension())
    mooseError(
        "The 'component' parameter in \"", name(), "\" block should be less than mesh dimension.");

  // Populate coupled variable information
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _disp_var[i] = coupled("displacements", i);
    _disp_older[i] = &coupledValueOlder("displacements", i);
  }
}

Real
NonReflectingBCMastodon::computeQpResidual()
{
  std::vector<Real> vel(3, 0.0);

  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    // accel = 1. / _beta * (((*_disp[i])[_qp] - (*_disp_old[i])[_qp]) / (_dt * _dt) -
    //                       (*_vel_old[i])[_qp] / _dt - (*_accel_old[i])[_qp] * (0.5 - _beta));
    // vel[i] =
    //     (*_vel_old[i])[_qp] + (_dt * (1. - _gamma)) * (*_accel_old[i])[_qp] + _gamma * _dt * accel;
    // vel[i] = (1. + _alpha) * vel[i] - _alpha * (*_vel_old[i])[_qp]; // HHT time integration
    // normal_vel += vel[i] * _normals[_qp](i);

    // central difference
    vel[i] = ((*_disp[i])[_qp] - (*_disp_older[i])[_qp])/(2*_dt);
    normal_vel += vel[i] * _normals[_qp](i);

  }
  // residual is test[i][_qp] *( density* V_p * normal component of velocity +
  // density * V_s* tangential component of velocity)
  return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed * normal_vel * _normals[_qp](_component) +
          _shear_wave_speed * (vel[_component] - normal_vel * _normals[_qp](_component)));
}

Real
NonReflectingBCMastodon::computeQpJacobian()
{
  return 0.0;
}

Real
NonReflectingBCMastodon::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}