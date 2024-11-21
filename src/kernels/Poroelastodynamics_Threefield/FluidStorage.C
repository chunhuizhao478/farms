//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidStorage.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "TimeIntegrator.h"

#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", FluidStorage);

InputParameters
FluidStorage::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");

  return params;
}

FluidStorage::FluidStorage(const InputParameters & parameters)
  : TimeKernel(parameters),
   _lumping(getParam<bool>("lumping")),
   _Biot_modulus(getMaterialProperty<Real>("biot_modulus")),
   _time_integrator(*_sys.getTimeIntegrator())
{
    _du_dot_du = &(_var.duDotDu());

    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
    

    if (_time_integrator.isLumped())
    {
    _u_dot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotFactorTag());
    }
}

Real
FluidStorage::computeQpResidual()
{
  if (_dt == 0)
  return 0.0;

  //return _test[_i][_qp] * _u_dot[_qp]* 1/_Biot_modulus[_qp];
  return  _test[_i][_qp] * 1/_Biot_modulus[_qp]*(*_u_dot_factor)[_qp];
 // return _test[_i][_qp] * 1/_Biot_modulus[_qp]*(_u[_qp]-(*_u_old)[_qp])/2/_dt;
    if (_time_integrator.isLumped())
  {
    return _test[_i][_qp] * 1/_Biot_modulus[_qp];
    for (unsigned int i = 0; i < _test.size(); ++i)
    this->_local_re(i) *= (*_u_dot_factor_dof)[_i] ;
   };
    
  

}
Real
FluidStorage::computeQpJacobian()
{
  if (_dt == 0)
  return 0.0;
  //return _test[_i][_qp] * _phi[_j][_qp] * 1/2/_dt * 1/_Biot_modulus[_qp];
  return _test[_i][_qp] * _phi[_j][_qp] * (*_du_dot_du)[_qp] * 1/_Biot_modulus[_qp];
}
void
FluidStorage::computeJacobian()
{
  if (_lumping)
  {
    prepareMatrixTag(_assembly, _var.number(), _var.number());

    precalculateJacobian();
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
          _local_ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

    accumulateTaggedLocalMatrix();
  }
  else
    TimeKernel::computeJacobian();
}