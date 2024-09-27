//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"
//c
//c This is the base Material class for implementing a rigid traction separation material model.
//c
class FarmsSlipWeakeningBase: public InterfaceMaterial
{
public:
  static InputParameters validParams();
  FarmsSlipWeakeningBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /* normal to the interface */
  const MooseArray<Point> & _normals;

  /**
   * method returning the traction and displacements in the interface coordinate system.
   */ 
  virtual Real computeTractionAndDisplacements() = 0;

  /**
   * method returning the traction derivitaves wrt local displacement jump.
   */ 
  virtual RealTensorValue computeTractionDerivatives() = 0;

  /**
   * method to compute rotation matrix from global to local coordinates
   */
  virtual RealTensorValue computeRotationMatrix(const RealVectorValue n_vector);

  /* helper functions */
  //normalize vector
  virtual RealVectorValue normalize(const RealVectorValue v);

  //perform cross product of two vectors
  virtual RealVectorValue crossProduct(const RealVectorValue a, const RealVectorValue b);

  /**
   * method to perform global to local transformation of a vector
   */
  virtual RealVectorValue GlobaltoLocalVector(const RealVectorValue globalVec, const RealTensorValue Q);

  /**
   * method to perform local to global transformation of a vector
   */
  virtual RealVectorValue LocaltoGlobalVector(const RealVectorValue localVec, const RealTensorValue Q);

  /**
   * method to perform global to local transformation of a matrix
   */
  virtual RealTensorValue GlobaltoLocalMatrix(const RealTensorValue globalMat, const RealTensorValue Q);

  /**
   * method to perform local to global transformation of a matrix
   */
  virtual RealTensorValue LocaltoGlobalMatrix(const RealTensorValue localMat, const RealTensorValue Q);  

  /* helper functions */
  //transpose matrix
  virtual RealTensorValue MatrixTranspose(const RealTensorValue Q);

  //matrix vector multiplication
  virtual RealVectorValue MatrixVectorMultiply(const RealTensorValue matrix, const RealVectorValue vec);

  //matrix matrix multiplication
  virtual RealTensorValue MatrixMatrixMultiply(const RealTensorValue matA, const RealTensorValue matB);

  //RankTwoTensor to RealTensorValue
  virtual RealTensorValue RankTwoTensor2RealTensorValue(const RankTwoTensor mat);

  /* the value of the traction in global (difference) */ 
  MaterialProperty<RealVectorValue> & _traction_on_interface;

  /* the value of the traction in global (total) */
  MaterialProperty<RealVectorValue> & _traction_total_global;

  /* the value of the traction in local (total) */
  MaterialProperty<RealVectorValue> & _traction_total_local;

  /* the tangent modulus on interface */
  MaterialProperty<RealTensorValue> & _material_tangent_modulus_on_interface;

  /* the rotation matrix */
  MaterialProperty<RealTensorValue> & _rotation_matrix;

  /* the displacement jump in global */
  MaterialProperty<RealVectorValue> & _displacement_jump_global;

  /* the displacement jump rate in global */
  MaterialProperty<RealVectorValue> & _displacement_jump_rate_global;  

  /* the displacements on the plus side (global) */
  MaterialProperty<RealVectorValue> & _displacements_plus_global;

  /* the displacements on the minus side (global) */
  MaterialProperty<RealVectorValue> & _displacements_minus_global;  

  /* the displacements on the plus side (local) */
  MaterialProperty<RealVectorValue> & _displacements_plus_local;

  /* the displacements on the minus side (local) */
  MaterialProperty<RealVectorValue> & _displacements_minus_local;

  /* the velocity on the plus side (local) */
  MaterialProperty<RealVectorValue> & _velocities_plus_local;

  /* the velocity on the minus side (local) */
  MaterialProperty<RealVectorValue> & _velocities_minus_local;  

  /* the absolute slip (local) */
  MaterialProperty<Real> & _absolute_slip;  

  /* the instaneous slip (local) */
  MaterialProperty<Real> & _slip_total; 

  /* marker for all elements below strength */
  MaterialProperty<RealVectorValue> & _below_strength_marker;

  /* the restoration force on primary surface */
  MaterialProperty<RealVectorValue> & _R_plus_local_vec;

  /* the restoration force on secondary surface */
  MaterialProperty<RealVectorValue> & _R_minus_local_vec;

};
