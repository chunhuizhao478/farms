#pragma once

#include "Material.h"

/**
 * Class to evaluate the phase field psi function and other related properties.
 * This class is used for non-AD implementation of the phase field model.
 * Created by Chunhui Zhao on 2025-5-11.
 */
class NDPhaseFieldPsiDerivatives : public Material
{
public:
  static InputParameters validParams();

  NDPhaseFieldPsiDerivatives(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

private:
  Real computeNormalizationConstant();
  Real normalizationIntegrand(const Real & d);

  void computeAlphaDerivatives(); // alpha(d) and dalpha/dd
  void computeGDerivatives(); // g(d) and dg/dd
  void computePsiDerivatives(); // psi(d) and dpsi/dd, d2psi/dd2

  const VariableValue & _d; // phase field damage variable
  const VariableValue & _psie_active; // active strain energy density

  const std::string _model_type; // model type: AT1, AT2, PF_CZM

  MaterialProperty<Real> & _c0; //normalization constant
  MaterialProperty<Real> & _alpha; //alpha(d)
  MaterialProperty<Real> & _dalpha_dd; //dalpha/dd derivative wrt d
  MaterialProperty<Real> & _d2alpha_dd2; //d2alpha/dd2 derivative wrt d
  MaterialProperty<Real> & _g; //g(d)
  MaterialProperty<Real> & _dg_dd; //dg/dd derivative wrt d
  MaterialProperty<Real> & _d2g_dd2; //d2g/dd2 derivative wrt d
  MaterialProperty<Real> & _psi; //free energy density
  MaterialProperty<Real> & _dpsi_dd; //dpsi/dd derivative wrt d
  MaterialProperty<Real> & _d2psi_dd2; //d2psi/dd2 derivative wrt d

  const Real _eta; // parameter in the degradation function
  const MaterialProperty<Real> & _Gc; // energy release rate
  const MaterialProperty<Real> & _l; // regularization length

  const Real _tolerance;
  const unsigned int _max_its;

  Real _c0_0;
};
