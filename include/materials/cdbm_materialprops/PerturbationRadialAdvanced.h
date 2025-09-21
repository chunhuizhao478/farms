//* Advanced radial perturbation material with smooth temporal and spatial tapers
#pragma once
#include "Material.h"

class PerturbationRadialAdvanced : public Material {
public:
  static InputParameters validParams();
  PerturbationRadialAdvanced(const InputParameters & params);
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
private:
  // material properties
  MaterialProperty<Real> & _damage_perturbation;
  MaterialProperty<Real> & _shear_stress_perturbation;
  MaterialProperty<Real> & _mean_stress_perturbation;
  // stateful old (requested to make stateful)
  const MaterialProperty<Real> & _damage_perturbation_old;
  const MaterialProperty<Real> & _shear_stress_perturbation_old;
  const MaterialProperty<Real> & _mean_stress_perturbation_old;
  MaterialProperty<std::vector<Real>> & _nucl_center_mat;
  MaterialProperty<Real> & _thickness_mat;
  MaterialProperty<Real> & _length_mat;

  // input parameters
  std::vector<Real> _nucl_center;
  Real _peak_value;
  Real _thickness;
  Real _length;
  Real _duration;
  std::string _perturbation_type; // damage | shear_stress | mean_stress
  Real _sigma_divisor;
  std::string _temporal_ramp;     // linear | cosine | smoothstep | quintic
  bool _hold_after_duration;      // keep peak after duration
  std::string _thickness_taper;   // none | cosine
  bool _zero_outside_band;        // if true, no carry over outside thickness

  // helpers
  Real rampFactor(Real t) const;
  Real thicknessFactor(Real y) const;
};
