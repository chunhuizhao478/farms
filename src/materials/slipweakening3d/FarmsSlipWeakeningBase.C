//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Assembly.h"
#include "FarmsSlipWeakeningBase.h"
InputParameters
FarmsSlipWeakeningBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Base class for rigid cohesive zone mateirla models");
  return params;
}

FarmsSlipWeakeningBase::FarmsSlipWeakeningBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _normals(_assembly.normals()),
    _traction_on_interface(declarePropertyByName<RealVectorValue>("traction_on_interface")),
    _traction_total_global(declarePropertyByName<RealVectorValue>("traction_total_global")),
    _traction_total_local(declarePropertyByName<RealVectorValue>("traction_total_local")),
    _material_tangent_modulus_on_interface(
        declareProperty<RealTensorValue>("material_tangent_modulus_on_interface")),
    _rotation_matrix(declareProperty<RealTensorValue>("rotation_matrix")),
    _displacement_jump_global(declareProperty<RealVectorValue>("displacement_jump_global")),
    _displacement_jump_rate_global(declareProperty<RealVectorValue>("displacement_jump_rate_global")),
    _displacements_plus_global(declareProperty<RealVectorValue>("displacements_plus_global")),
    _displacements_minus_global(declareProperty<RealVectorValue>("displacements_minus_global")),
    _displacements_plus_local(declareProperty<RealVectorValue>("displacements_plus_local")),
    _displacements_minus_local(declareProperty<RealVectorValue>("displacements_minus_local")),
    _velocities_plus_local(declareProperty<RealVectorValue>("velocities_plus_local")),
    _velocities_minus_local(declareProperty<RealVectorValue>("velocities_minus_local")),
    _absolute_slip(declareProperty<Real>("absolute_slip")),
    _slip_total(declareProperty<Real>("slip_total")),
    _below_strength_marker(declareProperty<RealVectorValue>("below_strength_marker")),
    _R_plus_local_vec(declareProperty<RealVectorValue>("R_plus_local_vec")),
    _R_minus_local_vec(declareProperty<RealVectorValue>("R_minus_local_vec"))
{
}

void
FarmsSlipWeakeningBase::computeQpProperties()
{

  /**
   * Compute rotation matrix
   */
  _rotation_matrix[_qp] = computeRotationMatrix(_normals[_qp]);

  /**
   * Compute traction on the interface and its derivatives
   * Compute displacements/velocities on the interface 
   */
  Real dummy = computeTractionAndDisplacements();
  _material_tangent_modulus_on_interface[_qp] = computeTractionDerivatives();

}

RealTensorValue
FarmsSlipWeakeningBase::computeRotationMatrix(const RealVectorValue n_vector)
{
  RealVectorValue n_vector_x(0.0,0.0,0.0);
  RealVectorValue n_vector_y(0.0,0.0,0.0);
  RealTensorValue Q_mat(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  
  Real sum_power = n_vector(0) * n_vector(0) + n_vector(1) * n_vector(1);

  if (sum_power < 1.e-8) {
      n_vector_x(0) = 1.0;
      n_vector_x(1) = 0.0;
      n_vector_x(2) = 0.0;
      n_vector_y(0) = 0.0;
      n_vector_y(1) = 1.0; 
      n_vector_y(2) = 0.0;
  } else {
      RealVectorValue global_z(0.0,0.0,0.0);
      global_z(0) = 0.0;
      global_z(1) = 0.0;
      global_z(2) = 1.0;
      n_vector_x = crossProduct(global_z, n_vector);
      n_vector_x = normalize(n_vector_x);
      n_vector_y = crossProduct(n_vector, n_vector_x);
      n_vector_y = normalize(n_vector_y);
  }

  for (unsigned int i = 0; i < 3; ++i) {
      Q_mat(i,0) = n_vector_x(i);
      Q_mat(i,1) = n_vector_y(i);
      Q_mat(i,2) = n_vector(i);
  }

  return Q_mat;
}

RealVectorValue
FarmsSlipWeakeningBase::normalize(const RealVectorValue v)
{ 
  RealVectorValue v_normalized(0.0,0.0,0.0);
  Real norm = std::sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2));
  v_normalized(0) = v(0) / norm;
  v_normalized(1) = v(1) / norm;
  v_normalized(2) = v(2) / norm;
  return v_normalized;
}

RealVectorValue
FarmsSlipWeakeningBase::crossProduct(const RealVectorValue a, const RealVectorValue b)
{
  RealVectorValue c_crossproduct(0.0,0.0,0.0);
  c_crossproduct(0) = a(1) * b(2) - a(2) * b(1);
  c_crossproduct(1) = a(2) * b(0) - a(0) * b(2);
  c_crossproduct(2) = a(0) * b(1) - a(1) * b(0);
  return c_crossproduct;
}

RealVectorValue
FarmsSlipWeakeningBase::GlobaltoLocalVector(const RealVectorValue globalVec, const RealTensorValue Q)
{
  RealTensorValue Q_T = MatrixTranspose(Q);
  return MatrixVectorMultiply(Q_T, globalVec);
}

RealVectorValue
FarmsSlipWeakeningBase::LocaltoGlobalVector(const RealVectorValue localVec, const RealTensorValue Q)
{
  return MatrixVectorMultiply(Q, localVec);
}

RealTensorValue
FarmsSlipWeakeningBase::GlobaltoLocalMatrix(const RealTensorValue globalMat, const RealTensorValue Q)
{
  RealTensorValue Q_T = MatrixTranspose(Q);
  return MatrixMatrixMultiply(MatrixMatrixMultiply(Q_T, globalMat), Q);
}

RealTensorValue
FarmsSlipWeakeningBase::LocaltoGlobalMatrix(const RealTensorValue localMat, const RealTensorValue Q)
{
  RealTensorValue Q_T = MatrixTranspose(Q);
  return MatrixMatrixMultiply(MatrixMatrixMultiply(Q, localMat), Q_T);
}

RealTensorValue
FarmsSlipWeakeningBase::MatrixTranspose(const RealTensorValue Q)
{
  RealTensorValue Q_T(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
          Q_T(j,i) = Q(i,j);
      }
  }
  return Q_T;
}

RealVectorValue
FarmsSlipWeakeningBase::MatrixVectorMultiply(const RealTensorValue matrix, const RealVectorValue vec)
{
  RealVectorValue result_vec(0.0,0.0,0.0);
  for (unsigned int i = 0; i < 3; ++i) {
      result_vec(i) = matrix(i,0) * vec(0) + matrix(i,1) * vec(1) + matrix(i,2) * vec(2);
  }
  return result_vec;
}

RealTensorValue
FarmsSlipWeakeningBase::MatrixMatrixMultiply(const RealTensorValue matA, const RealTensorValue matB)
{
  RealTensorValue result_mat(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
          for (unsigned int k = 0; k < 3; ++k) {
              result_mat(i,j) += matA(i,k) * matB(k,j);
          }
      }
  }
  return result_mat;
}

RealTensorValue
FarmsSlipWeakeningBase::RankTwoTensor2RealTensorValue(const RankTwoTensor mat)
{
  RealTensorValue result_mat(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      result_mat(i, j) = mat(i, j);
    }
  }
  return result_mat;
}