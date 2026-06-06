//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "NDSmallDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class NDSmallDeformationIsotropicElasticity : public NDSmallDeformationElasticityModel,
                                            public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  NDSmallDeformationIsotropicElasticity(const InputParameters & parameters);

  virtual RankTwoTensor computeStress(const RankTwoTensor & strain) override;

  // Compute Jacobian of the stress w/r/t strain
  virtual RankFourTensor computeJacobian(const RankTwoTensor & strain) override;

protected:
private:
  // @{ Decomposition methods
  virtual RankTwoTensor computeStressNoDecomposition(const RankTwoTensor & strain);
  virtual RankTwoTensor computeStressSpectralDecomposition(const RankTwoTensor & strain);
  virtual RankTwoTensor computeStressVolDevDecomposition(const RankTwoTensor & strain);
  // @}

  //add jacobain of the stress w/r/t strain with decomposition methods
  virtual RankFourTensor computeJacobianNoDecomposition(const RankTwoTensor & strain);
  virtual RankFourTensor computeJacobianSpectralDecomposition(const RankTwoTensor & strain);
  virtual RankFourTensor computeJacobianVolDevDecomposition(const RankTwoTensor & strain);

  // @{ add additional functions for porous flow coupling
  virtual void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);
  virtual void updatePermeabilityForCracking();
  // @}

  // @{ Helper functions
  Real Macaulay(const Real x, const bool deriv = false);
  std::vector<Real> Macaulay(const std::vector<Real> & v, const bool deriv = false);
  RankTwoTensor spectralDecomposition(const RankTwoTensor & r2t);
  // Return the characteristic length h_c used in the normal-strain permeability
  // model. Source depends on `characteristic_length_type`.
  Real getCharacteristicLength() const;
  // Liu et al. 2024 CMAME eqs. (29)-(30): return the unit eigenvector e_1 of the
  // MAXIMUM principal strain of `strain` (the strain-based crack normal n_F).
  // The maximum principal strain value eps_1 is returned via `eps1`; callers gate
  // the fracture-permeability enhancement on eps_1 > 0 (tensile opening) so that a
  // closed/compressed crack -- where the plane-strain max eigenvalue is the
  // out-of-plane eps_zz = 0 and e_1 = e_z -- is not given a spurious aperture.
  RealVectorValue maxPrincipalStrainDirection(const RankTwoTensor & strain,
                                              Real & eps1) const;
  // Bulk modulus extracted from the degraded SPECTRAL elastic tangent C(d, eps)
  // via the volumetric contraction K = (1/9) I:C:I. Builds C with the TRUE
  // I (x) I (outerProduct), NOT the diagonal-only RankFourTensor(initIdentity).
  Real computeSpectralBulkModulus(const RankTwoTensor & strain);
  // @}

  // Compute g and its derivatives
  void computeGDerivatives();

  /// The bulk modulus
  const MaterialProperty<Real> & _K;

  /// The shear modulus
  const MaterialProperty<Real> & _G;

  // Model type
  const std::string _model_type;

  /// Name of parameters for PF_CZM model
  //------------------------------------//
  // Store pointers to material properties that are only needed for PF_CZM
  const MaterialProperty<Real> * _a1_prop;
  const MaterialProperty<Real> * _a2_prop;
  const MaterialProperty<Real> * _a3_prop;
  const MaterialProperty<Real> * _p_prop;

  // Store the property names
  const MaterialPropertyName _a1_name;
  const MaterialPropertyName _a2_name;
  const MaterialPropertyName _a3_name;
  const MaterialPropertyName _p_name;
  //------------------------------------//

  /// Name of the phase-field variable
  const VariableValue & _d;

  // @{ Strain energy density and its derivative w/r/t damage
  MaterialProperty<Real> & _psie;
  MaterialProperty<Real> & _psie_active;
  MaterialProperty<Real> & _psie_inactive;
  MaterialProperty<Real> & _dpsie_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  MaterialProperty<Real> & _g;
  MaterialProperty<Real> & _dg_dd;
  MaterialProperty<Real> & _d2g_dd2;
  // @}

  // Constants
  const Real _eta;

  /// Decomposittion types
  const enum class Decomposition { none, spectral, voldev } _decomposition;

  /// Add additional material properties for porous flow coupling
  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  MaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  /// @brief define the effective permeability
  MaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old;

  const bool _porous_flow_coupling; // flag to indicate if porous flow coupling is enabled
  const Real _intrinsic_permeability;

  // Exponential permeability model
  const bool _exponential_permeability_model; // flag to indicate if exponential permeability model is used
  const Real _coeff_b; // coefficient for the exponential function in the effective permeability
  // Darcy-Poiseuille permeability model
  const bool _darcy_poiseuille_permeability_model; // flag to indicate if Darcy-Poiseuille permeability model is used
  const Real _wc; // characteristic width for the Darcy-Poiseuille model
  const Real _perm_exponent; // exponent for the Darcy-Poiseuille model

  /// Permeability model selection (canonical storage for the enum); legacy
  /// boolean flags above are mapped into this enum in the constructor.
  enum class PermeabilityModel { none, exponential, darcy_poiseuille, normal_strain };
  const PermeabilityModel _permeability_model;

  /// Source of the unit crack normal n_d for the normal-strain model.
  enum class CrackNormalSource { damage_gradient, principal_strain };
  const CrackNormalSource _normal_source;

  /// Source of the characteristic length h_c in w_c = h_c * |1 + eps_nn|.
  enum class LcType { element_size, regularization_length, constant };
  const LcType _lc_type;

  // Damage gradient (used by the normal-strain model when
  // _normal_source == damage_gradient). Bound to coupledGradient("phase_field")
  // in the constructor; zero-cost when not used.
  const VariableGradient & _grad_d;

  // Regularization length material property (only valid when
  // _lc_type == regularization_length; nullptr otherwise).
  const MaterialProperty<Real> * _l_mat_prop;

  // Local element size coupled variable (only valid when
  // _lc_type == element_size; nullptr otherwise).
  const VariableValue * _h_elem;

  // Constant characteristic length (only valid when _lc_type == constant).
  const Real _lc_const;

  // Anisotropy, threshold, roughness, tolerance for the normal-strain model.
  const bool _perm_anisotropic;
  const Real _d_perm_threshold;
  const Real _fc;
  const Real _grad_d_tol;

  // Regularized crack-normal option (damage_gradient source only). When true,
  // n_d = grad(d) / (|grad(d)| + _crack_normal_reg_eps) replaces the hard
  // |grad(d)| > _grad_d_tol cutoff. As |grad(d)| -> 0 (crack core, d -> 1),
  // n_d -> 0 so the tangential projector (I - n_d (x) n_d) -> I, giving
  // isotropic fracture permeability at the core instead of matrix-perm
  // fallback. _crack_normal_reg_eps > 0 (range-checked) avoids div-by-zero.
  const bool _regularize_crack_normal;
  const Real _crack_normal_reg_eps;

  // Residual aperture w_r (Heider 2021 eq. 46 closed-crack branch). When the
  // open-crack aperture w_c shrinks below w_r under crack closure, w_h floors
  // at f_c * w_r * chi_d so the fracture conductivity does not collapse to
  // zero. Default 0.0 reproduces the open-only formulation (legacy behavior).
  const Real _w_res;

  // Heider eq. (47) uses the total linearized strain ε^S = ½(∇u + ∇^T u),
  // i.e. the "mechanical_strain" property declared by ComputeSmallStrain. We
  // bind it explicitly here (equal to _elastic_strain in pure elasticity, but
  // the right thing if a plasticity model is later attached and `_elastic_strain`
  // diverges from the kinematic strain).
  const MaterialProperty<RankTwoTensor> * _total_strain;

  // Damaged solid bulk compliance C_s(d) = 1 / (g(d) * K)
  MaterialProperty<Real> & _solid_bulk_compliance_damaged;

  /// Bulk modulus K = (1/9) I:C:I of the degraded SPECTRAL elastic tangent.
  MaterialProperty<Real> & _bulk_modulus_degraded;
};
