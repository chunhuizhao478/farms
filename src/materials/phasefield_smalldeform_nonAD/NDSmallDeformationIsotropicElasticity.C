//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDSmallDeformationIsotropicElasticity.h"
#include "RaccoonUtils.h"
#include "SymmetricRankFourTensor.h"

registerMooseObject("farmsApp", NDSmallDeformationIsotropicElasticity);

InputParameters
NDSmallDeformationIsotropicElasticity::validParams()
{
  InputParameters params = NDSmallDeformationElasticityModel::validParams();
  params.addClassDescription("Isotropic elasticity under small strain asumptions.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  //material property names
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("strain_energy_density_active",
                                        "psie_active",
                                        "Name of the active strain energy density");
  params.addParam<MaterialPropertyName>("strain_energy_density_inactive",
                                        "psie_inactive",
                                        "Name of the inactive strain energy density");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density_derivative",
      "dpsie_dd",
      "Name of the strain energy density derivative w/r/t damage");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MaterialPropertyName>(
      "degradation_function_derivative",
      "dg_dd",
      "Name of the degradation function derivative w/r/t damage");
  params.addParam<MaterialPropertyName>(
      "degradation_function_second_derivative",
      "d2g_dd2",
      "Name of the degradation function second derivative w/r/t damage");
  params.addParam<MooseEnum>(
      "decomposition", MooseEnum("NONE SPECTRAL VOLDEV", "NONE"), "The decomposition method");

  params.addParam<std::string>(
      "model_type",
      "AT1",
      "The type of the model: AT1, AT2, PF_CZM");

  params.addRequiredParam<Real>("eta", "Parameter in the degradation function");
  params.addParam<bool>("porous_flow_coupling", false, "Enable porous flow coupling");
  params.addParam<Real>("intrinsic_permeability", 5e-19, "Intrinsic permeability in m^2");

  //Permeability models
  //Exponential permeability model
  params.addParam<bool>("exponential_permeability_model", false,
                        "Use an exponential function for the effective permeability");
  params.addParam<Real>("coeff_b", -1.0,
                        "Coefficient for the exponential function in the effective permeability");
  //Darcy-Poiseuille permeability model
  params.addParam<bool>("darcy_poiseuille_permeability_model",
                        false,
                        "Use Darcy-Poiseuille model for the effective permeability");
  params.addParam<Real>("wc", -1.0, "ultimate crack width for Darcy-Poiseuille model");
  params.addParam<Real>("perm_exponent", -1.0,
                        "Exponent for the Darcy-Poiseuille model for the effective permeability");

  // Canonical permeability-model enum. When set to anything other than "none"
  // it overrides the legacy booleans above.
  params.addParam<MooseEnum>(
      "permeability_model",
      MooseEnum("none exponential darcy_poiseuille normal_strain", "none"),
      "Permeability enhancement model to apply when porous_flow_coupling=true");

  // Normal-strain (Heider 2021 eqs. 46-48) parameters
  params.addParam<MooseEnum>(
      "crack_normal_source",
      MooseEnum("damage_gradient principal_strain", "damage_gradient"),
      "Source of the unit crack normal n_F. 'damage_gradient' uses "
      "grad(d)/|grad(d)| (Heider 2021 eq. 46), optionally regularized via "
      "regularize_crack_normal. 'principal_strain' uses the strain-based normal "
      "n_F = e_1, the eigenvector of the maximum principal strain of the model's "
      "mechanical strain (Liu et al. 2024 CMAME eqs. 29-30); it is unit-norm "
      "everywhere and needs no regularization.");
  params.addParam<Real>(
      "damage_gradient_tolerance", 1e-30,
      "When |grad(d)| is below this tolerance, fall back to K = K_poro "
      "(no enhancement). Avoids division by zero in undamaged regions. "
      "Ignored when regularize_crack_normal = true.");
  params.addParam<bool>(
      "regularize_crack_normal", false,
      "When true (and crack_normal_source = damage_gradient), compute the crack "
      "normal as n_d = grad(d) / (|grad(d)| + eps) instead of grad(d)/|grad(d)| "
      "with a hard |grad(d)| cutoff. As |grad(d)| -> 0 at the fully-damaged "
      "crack core (d -> 1), n_d -> 0, so the tangential projector "
      "(I - n_d (x) n_d) -> I and the fracture permeability becomes ISOTROPIC "
      "there instead of falling back to the matrix permeability k0*I. The "
      "epsilon is set by `crack_normal_regularization`. Default false preserves "
      "the legacy hard-cutoff behavior. NOTE: as n_d -> 0 the normal strain "
      "eps_nn = n_d.eps.n_d -> 0 too, so at the core the aperture w_c -> h_c "
      "loses its strain dependence and K_frac becomes the strain-independent "
      "isotropic value d^b*(h_c^2/12)*I.");
  params.addRangeCheckedParam<Real>(
      "crack_normal_regularization", 1e-8, "crack_normal_regularization > 0.0",
      "Regularization epsilon added to |grad(d)| in the denominator of the "
      "regularized crack normal n_d = grad(d)/(|grad(d)| + eps) (only used when "
      "regularize_crack_normal = true). Has units of 1/length (same as "
      "|grad(d)|). Choose it small relative to the typical |grad(d)| ~ 1/l so "
      "that n_d stays ~unit away from the crack core, yet large enough that "
      "n_d -> 0 (isotropic permeability) as the core (grad(d) -> 0) is "
      "approached. Must be > 0 to avoid division by zero.");
  params.addParam<MooseEnum>(
      "characteristic_length_type",
      MooseEnum("element_size regularization_length constant", "element_size"),
      "Source of h_c in w_c = h_c*(...). 'element_size' matches the paper's "
      "1-D line element (default). 'regularization_length' uses l. 'constant' "
      "uses a user-supplied value.");
  params.addParam<MaterialPropertyName>(
      "regularization_length_name", "l",
      "Name of the material property holding the regularization length "
      "(used when characteristic_length_type = regularization_length).");
  params.addCoupledVar(
      "element_size_variable",
      "Name of an aux variable holding the local element size h (used when "
      "characteristic_length_type = element_size).");
  params.addParam<Real>(
      "characteristic_length_value", -1.0,
      "Constant h_c (used when characteristic_length_type = constant).");
  params.addParam<bool>(
      "permeability_anisotropic", true,
      "If true, K_frac = (w^2/12)*(I - n_d (x) n_d) (Heider eq. 46). "
      "If false, K_frac = (w^2/12)*I (isotropic fallback).");
  params.addParam<Real>(
      "damage_threshold_for_permeability", 0.5,
      "Heaviside gate chi_d = H(d - threshold) (Heider eq. 46). "
      "K_frac is zero where d < threshold. Default 0.5 matches the paper.");
  params.addParam<Real>(
      "correction_factor_fc", 1.0,
      "Roughness correction factor f_c in w_h = f_c*w_c*chi_d (Heider eq. 46). "
      "Default 1.0 (smooth walls).");
  params.addRangeCheckedParam<Real>(
      "residual_aperture", 0.0, "residual_aperture >= 0.0",
      "Residual (closed-crack) aperture w_r in Heider eq. (46): "
      "w_h = max{(f_c*w_c)*chi_d, (f_c*w_r)*chi_d}. When the open-crack "
      "aperture collapses under closure (eps_nn -> -1), w_h floors at "
      "f_c*w_r*chi_d so the fracture conductivity decays to "
      "(f_c*w_r)^2/12 instead of zero. Default 0.0 reproduces the open-only "
      "formulation. Typical jointed-rock value: 1e-5 m.");

  //only used for PF_CZM model
  params.addParam<MaterialPropertyName>("a1", "", "a1 (only needed for PF_CZM)");
  params.addParam<MaterialPropertyName>("a2", "", "a2 (only needed for PF_CZM)");
  params.addParam<MaterialPropertyName>("a3", "", "a3 (only needed for PF_CZM)");
  params.addParam<MaterialPropertyName>("p",  "",  "p (only needed for PF_CZM)");

  return params;
}

NDSmallDeformationIsotropicElasticity::NDSmallDeformationIsotropicElasticity(
    const InputParameters & parameters)
  : NDSmallDeformationElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _K(getMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),

    // model type
    _model_type(getParam<std::string>("model_type")),

    // Only retrieve these properties if we're using PF_CZM model
    _a1_prop(_model_type == "PF_CZM" ?
             &getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("a1")) : nullptr),
    _a2_prop(_model_type == "PF_CZM" ?
             &getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("a2")) : nullptr),
    _a3_prop(_model_type == "PF_CZM" ?
             &getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("a3")) : nullptr),
    _p_prop(_model_type == "PF_CZM" ?
             &getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("p")) : nullptr),
    // Store the property names too (only used if model_type is PF_CZM)
    _a1_name(getParam<MaterialPropertyName>("a1")),
    _a2_name(getParam<MaterialPropertyName>("a2")),
    _a3_name(getParam<MaterialPropertyName>("a3")),
    _p_name(getParam<MaterialPropertyName>("p")),

    // The phase-field variable
    _d(coupledValue("phase_field")),

    // The strain energy density and its derivatives
    _psie(declareProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density"))),
    _psie_active(declareProperty<Real>(getParam<MaterialPropertyName>(
        "strain_energy_density_active"))),
    _psie_inactive(declareProperty<Real>(getParam<MaterialPropertyName>(
        "strain_energy_density_inactive"))),
    _dpsie_dd(declareProperty<Real>(getParam<MaterialPropertyName>(
        "strain_energy_density_derivative"))),

    // The degradation function and its derivatives
    _g(declareProperty<Real>(getParam<MaterialPropertyName>("degradation_function"))),
    _dg_dd(declareProperty<Real>(getParam<MaterialPropertyName>(
        "degradation_function_derivative"))),
    _d2g_dd2(declareProperty<Real>(getParam<MaterialPropertyName>(
        "degradation_function_second_derivative"))),

    // Constants
    _eta(getParam<Real>("eta")),

    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>()),

    //porous flow coupling
    _crack_rotation(declareProperty<RankTwoTensor>("crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOldByName<RankTwoTensor>("crack_rotation")),
    _effective_perm(declareProperty<RealTensorValue>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RealTensorValue>("effective_perm")),
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),

    // Exponential permeability model
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    // Darcy-Poiseuille permeability model
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent")),
    // Canonical enum, overrides legacy bools if non-"none"
    _permeability_model([&]() -> PermeabilityModel {
      const std::string choice = getParam<MooseEnum>("permeability_model");
      if (choice == "exponential") return PermeabilityModel::exponential;
      if (choice == "darcy_poiseuille") return PermeabilityModel::darcy_poiseuille;
      if (choice == "normal_strain") return PermeabilityModel::normal_strain;
      // "none": fall back to legacy boolean flags (backward compatibility).
      if (getParam<bool>("exponential_permeability_model"))
        return PermeabilityModel::exponential;
      if (getParam<bool>("darcy_poiseuille_permeability_model"))
        return PermeabilityModel::darcy_poiseuille;
      return PermeabilityModel::none;
    }()),
    _normal_source(getParam<MooseEnum>("crack_normal_source") == "damage_gradient"
                     ? CrackNormalSource::damage_gradient
                     : CrackNormalSource::principal_strain),
    _lc_type([&]() -> LcType {
      const std::string choice = getParam<MooseEnum>("characteristic_length_type");
      if (choice == "regularization_length") return LcType::regularization_length;
      if (choice == "constant") return LcType::constant;
      return LcType::element_size;
    }()),
    _grad_d(coupledGradient("phase_field")),
    _l_mat_prop(nullptr),
    _h_elem(nullptr),
    _lc_const(getParam<Real>("characteristic_length_value")),
    _perm_anisotropic(getParam<bool>("permeability_anisotropic")),
    _d_perm_threshold(getParam<Real>("damage_threshold_for_permeability")),
    _fc(getParam<Real>("correction_factor_fc")),
    _grad_d_tol(getParam<Real>("damage_gradient_tolerance")),
    _regularize_crack_normal(getParam<bool>("regularize_crack_normal")),
    _crack_normal_reg_eps(getParam<Real>("crack_normal_regularization")),
    _w_res(getParam<Real>("residual_aperture")),
    _total_strain(nullptr),
    _solid_bulk_compliance_damaged(declareProperty<Real>("solid_bulk_compliance_damaged")),
    _bulk_modulus_degraded(declareProperty<Real>("bulk_modulus_degraded"))
{
  // Warn only when the new enum and the legacy boolean flags *disagree*.
  // Same-choice redundancy (e.g. enum=darcy_poiseuille + legacy_dp=true) is a
  // common user mistake and should be silent; only emit when the enum picks
  // one branch while the other legacy bool is also true (or normal_strain is
  // picked but a legacy bool is also on).
  {
    const std::string enum_choice = getParam<MooseEnum>("permeability_model");
    const bool leg_exp = getParam<bool>("exponential_permeability_model");
    const bool leg_dp = getParam<bool>("darcy_poiseuille_permeability_model");
    const bool enum_is_exp = (enum_choice == "exponential");
    const bool enum_is_dp = (enum_choice == "darcy_poiseuille");
    const bool enum_is_ns = (enum_choice == "normal_strain");
    const bool inconsistent =
        (enum_is_exp && leg_dp) ||
        (enum_is_dp && leg_exp) ||
        (enum_is_ns && (leg_exp || leg_dp));
    if (inconsistent)
      mooseDoOnce(mooseWarning(
          "permeability_model = '", enum_choice,
          "' conflicts with the legacy boolean flags "
          "(exponential=", leg_exp, ", darcy_poiseuille=", leg_dp,
          "); the permeability_model enum wins."));
  }

  // Look up optional l_c sources based on the selected lc_type.
  // Only required when the normal-strain permeability model is active; legacy
  // branches (exponential, darcy_poiseuille) do not consume h_c.
  if (_permeability_model == PermeabilityModel::normal_strain)
  {
    if (_lc_type == LcType::regularization_length)
      _l_mat_prop = &getMaterialPropertyByName<Real>(
          getParam<MaterialPropertyName>("regularization_length_name"));
    if (_lc_type == LcType::element_size)
    {
      if (!isCoupled("element_size_variable"))
        paramError("characteristic_length_type",
                   "characteristic_length_type = element_size requires "
                   "`element_size_variable` to be coupled. Set "
                   "`element_size_variable = <aux var>` or pick a different "
                   "characteristic_length_type.");
      _h_elem = &coupledValue("element_size_variable");
    }
    // Heider eq. (47) explicitly uses the total linearized strain
    // ε^S = ½(∇u + ∇^T u). Bind it via the "mechanical_strain" property declared
    // by ComputeSmallStrain. In pure elasticity this equals _elastic_strain
    // exactly (see NDSmallDeformationElasticityModel.C:67). If a plasticity
    // model is later attached, both this _total_strain read AND the
    // _elastic_strain eigendecomposition in computeCrackStrainAndOrientation
    // (used for the principal-strain fallback for n_d) need to be revisited
    // for full Heider-eq.-(47) consistency — n_d would be derived from
    // (mechanical_strain - plastic_strain) while eps_nn would be derived
    // from mechanical_strain.
    if (!hasMaterialPropertyByName<RankTwoTensor>(prependBaseName("mechanical_strain")))
      paramError("permeability_model",
                 "permeability_model = normal_strain requires the kinematic "
                 "strain tensor 'mechanical_strain' to be declared by a "
                 "ComputeSmallStrain (or compatible) material in [Materials]. "
                 "Add `[strain] type = ComputeSmallStrain []` to your input.");
    _total_strain = &getMaterialPropertyByName<RankTwoTensor>(
        prependBaseName("mechanical_strain"));
  }

  // Generic validity checks (apply to any enabled permeability path).
  if (_porous_flow_coupling &&
      _permeability_model == PermeabilityModel::none)
    paramError("porous_flow_coupling",
               "Porous flow coupling is enabled, but no permeability model is "
               "selected. Set `permeability_model` to exponential, "
               "darcy_poiseuille, or normal_strain, or enable the corresponding "
               "legacy boolean flag.");

  if (_permeability_model == PermeabilityModel::darcy_poiseuille &&
      (_wc <= 0.0 || _perm_exponent <= 0.0))
    paramError("darcy_poiseuille_permeability_model",
               "Darcy-Poiseuille permeability model is enabled, but wc and perm_exponent "
               "must be positive values. Please check the input parameters.");

  if (_permeability_model == PermeabilityModel::exponential && _coeff_b <= 0.0)
    paramError("exponential_permeability_model",
               "Exponential permeability model is enabled, but coeff_b must be a positive value. "
               "Please check the input parameters.");

  // Normal-strain model validity checks (Heider 2021 eqs. 46-48).
  if (_permeability_model == PermeabilityModel::normal_strain)
  {
    if (!_porous_flow_coupling)
      paramError("permeability_model",
                 "permeability_model = normal_strain requires "
                 "porous_flow_coupling = true.");
    if (_perm_exponent <= 0.0)
      paramError("perm_exponent",
                 "permeability_model = normal_strain requires perm_exponent > 0.");
    if (_lc_type == LcType::constant && _lc_const <= 0.0)
      paramError("characteristic_length_value",
                 "characteristic_length_type = constant requires "
                 "characteristic_length_value > 0.");
    if (_fc <= 0.0 || _fc > 1.0)
      paramError("correction_factor_fc",
                 "correction_factor_fc must satisfy 0 < f_c <= 1.");
    if (_d_perm_threshold < 0.0 || _d_perm_threshold >= 1.0)
      paramError("damage_threshold_for_permeability",
                 "damage_threshold_for_permeability must satisfy "
                 "0 <= threshold < 1.");
    if (_regularize_crack_normal &&
        _normal_source != CrackNormalSource::damage_gradient)
      paramError("regularize_crack_normal",
                 "regularize_crack_normal = true only applies to "
                 "crack_normal_source = damage_gradient. The principal-strain "
                 "normal is already a well-defined unit eigenvector and needs "
                 "no regularization.");
    // The default crack_normal_regularization (1e-8 /length) is far below the
    // phase-field gradient scale |grad(d)| ~ 1/l, so leaving it at the default
    // makes the shrink factor |grad d|/(|grad d|+eps) ~ 1 at every quadrature
    // point and the regularization a silent no-op. Require an explicit value so
    // the user picks eps on the gradient scale they want to isotropize.
    if (_regularize_crack_normal &&
        !isParamSetByUser("crack_normal_regularization"))
      paramError("crack_normal_regularization",
                 "regularize_crack_normal = true requires "
                 "crack_normal_regularization to be set explicitly. The default "
                 "(1e-8) is far below the phase-field gradient scale "
                 "|grad(d)| ~ 1/l, so the regularization would have no effect at "
                 "the quadrature points. Set eps on the order of the near-core "
                 "|grad(d)| you want to isotropize (e.g. a fraction of 1/l).");
  }
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStress(const RankTwoTensor & strain)
{
  RankTwoTensor stress;

  // Evaluate g and derivatives
  computeGDerivatives();

  if (_decomposition == Decomposition::none)
    stress = computeStressNoDecomposition(strain);
  else if (_decomposition == Decomposition::spectral)
    stress = computeStressSpectralDecomposition(strain);
  else if (_decomposition == Decomposition::voldev)
    stress = computeStressVolDevDecomposition(strain);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  // Bulk modulus from the degraded elastic tangent. g(d) is already set by
  // computeGDerivatives() above. The SPECTRAL tangent is contracted exactly
  // (K_eff = (1/9) I:C:I); for other decompositions we report g*K (exact for
  // NONE) and warn once, rather than mislead with the spectral split value.
  if (_decomposition == Decomposition::spectral)
    _bulk_modulus_degraded[_qp] = computeSpectralBulkModulus(strain);
  else
  {
    mooseDoOnce(mooseWarning(
        "bulk_modulus_degraded uses g*K for decomposition != SPECTRAL; the "
        "exact contraction is only computed for the spectral tangent."));
    _bulk_modulus_degraded[_qp] = _g[_qp] * _K[_qp];
  }

  // Reciprocal damaged solid bulk compliance C_s(d) = 1 / K_eff (populates the
  // previously never-assigned property).
  _solid_bulk_compliance_damaged[_qp] = 1.0 / std::max(_bulk_modulus_degraded[_qp], 1e-30);

  return stress;
}

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobian(const RankTwoTensor & strain)
{
  RankFourTensor Jacobian;

  // Evaluate g and derivatives
  computeGDerivatives();

  if (_decomposition == Decomposition::none)
    Jacobian = computeJacobianNoDecomposition(strain);
  else if (_decomposition == Decomposition::spectral)
    Jacobian = computeJacobianSpectralDecomposition(strain);
  else if (_decomposition == Decomposition::voldev)
    Jacobian = computeJacobianVolDevDecomposition(strain);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return Jacobian;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressNoDecomposition(const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress = _g[_qp] * stress_intact;

  _psie_active[_qp] = 0.5 * stress_intact.doubleContraction(strain);
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressSpectralDecomposition(
    const RankTwoTensor & strain)
{
  const Real lambda = _K[_qp] - 2 * _G[_qp] / LIBMESH_DIM;
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real strain_tr = strain.trace();
  Real strain_tr_pos = NDSmallDeformationIsotropicElasticity::Macaulay(strain_tr);

  // Spectral decomposition
  RankTwoTensor strain_pos = NDSmallDeformationIsotropicElasticity::spectralDecomposition(strain);

  // Stress
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress_pos = lambda * strain_tr_pos * I2 + 2 * _G[_qp] * strain_pos;
  RankTwoTensor stress_neg = stress_intact - stress_pos;
  RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  Real psie_intact =
      0.5 * lambda * strain_tr * strain_tr + _G[_qp] * strain.doubleContraction(strain);
  _psie_active[_qp] = 0.5 * lambda * strain_tr_pos * strain_tr_pos +
                      _G[_qp] * strain_pos.doubleContraction(strain_pos);
  Real psie_inactive = psie_intact - _psie_active[_qp];
  _psie_inactive[_qp] = psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  //Porous flow coupling
  /* Populate _crack_rotation[_qp] via computeCrackStrainAndOrientation's
     side effect. The principal-strain out-param is not consumed by the
     permeability update (which reads _crack_rotation and _grad_d directly),
     so the local is intentionally unused. */
  RealVectorValue unused_principal_strains;
  computeCrackStrainAndOrientation(unused_principal_strains);

  // Compute effective permeability.
  updatePermeabilityForCracking();

  return stress;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressVolDevDecomposition(
    const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);

  // Volumetric-deviatoric decomposition
  Real strain_tr = strain.trace();
  Real strain_tr_pos = NDSmallDeformationIsotropicElasticity::Macaulay(strain_tr);
  Real strain_tr_neg = strain_tr - strain_tr_pos;
  RankTwoTensor strain_dev = strain.deviatoric();

  // Stress
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress_neg = _K[_qp] * strain_tr_neg * I2;
  RankTwoTensor stress_pos = stress_intact - stress_neg;
  RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  Real psie_intact =
      0.5 * _K[_qp] * strain_tr * strain_tr + _G[_qp] * strain_dev.doubleContraction(strain_dev);
  Real psie_inactive = 0.5 * _K[_qp] * strain_tr_neg * strain_tr_neg;
  _psie_active[_qp] = psie_intact - psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

//Jacobian of the stress w/r/t strain with decomposition methods
// no decomposition: σ = g·(K tr ε I + 2G ε_dev)
//
// ⇒ J_ijkl = g · [ K δ_ij δ_kl + 2G ( ½(δ_ik δ_jl + δ_il δ_jk) − 1/3 δ_ij δ_kl ) ]
RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianNoDecomposition(
    const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4 = RankFourTensor(RankFourTensor::initIdentity);
  RankFourTensor I4_sym = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

  RankFourTensor Jacobian_intact = _K[_qp] * I4 + 2 * _G[_qp] * (I4_sym - I4 / 3.0);
  RankFourTensor Jacobian = _g[_qp] * Jacobian_intact;

  return Jacobian;
}

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianSpectralDecomposition(
    const RankTwoTensor & strain)
{
  //--------------------------------------------------------------------
  // 1.  Some handy constants and fourth–order identity tensors
  //--------------------------------------------------------------------
  const Real  lambda = _K[_qp] - 2.0 * _G[_qp] / LIBMESH_DIM;

  const RankFourTensor I4 = RankFourTensor(RankFourTensor::initIdentity);
  const RankFourTensor I4_sym = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

  //--------------------------------------------------------------------
  // 2.  Intact (undegraded) isotropic‑elastic tangent
  //--------------------------------------------------------------------
  RankFourTensor C_intact =
      _K[_qp] * I4 + 2.0 * _G[_qp] * (I4_sym - I4 / 3.0);

  //--------------------------------------------------------------------
  // 3.  Positive‑part projector  P⁺  and volumetric Heaviside term
  //--------------------------------------------------------------------
  RankTwoTensor   eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  RankFourTensor P_pos =
      strain.positiveProjectionEigenDecomposition(eigvals, eigvecs); // P⁺₍ᵢⱼₖₗ₎

  const Real H_tr = (strain.trace() > 0.0) ? 1.0 : 0.0;              // H(tr ε)

  RankFourTensor C_pos =
      lambda * H_tr * I4           // volumetric part
    + 2.0 * _G[_qp] * P_pos;       // deviatoric spectral part

  //--------------------------------------------------------------------
  // 4.  Final tangent  C = C_intact + (g-1) C_pos
  //--------------------------------------------------------------------
  RankFourTensor Jacobian =
      C_intact + (_g[_qp] - 1.0) * C_pos;

  return Jacobian;
}

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianVolDevDecomposition(
    const RankTwoTensor & strain)
{
  RankFourTensor Jacobian;
  return Jacobian;
}

//helper function, grab from RaccoonUtils, make it non-AD
Real
NDSmallDeformationIsotropicElasticity::Macaulay(const Real x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}

std::vector<Real>
NDSmallDeformationIsotropicElasticity::Macaulay(const std::vector<Real> & v, const bool deriv)
{
  std::vector<Real> m = v;
  for (auto & x : m)
    x = Macaulay(x, deriv);
  return m;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::spectralDecomposition(const RankTwoTensor & r2t)
{
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  r2t.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  RankTwoTensor eigvals_pos;
  eigvals_pos.fillFromInputVector(Macaulay(eigvals));
  return eigvecs * eigvals_pos * eigvecs.transpose();
}

RealVectorValue
NDSmallDeformationIsotropicElasticity::maxPrincipalStrainDirection(
    const RankTwoTensor & strain, Real & eps1) const
{
  // Liu et al. 2024 CMAME eq. (29): eps = sum_i eps_i e_i, with the eigenvalues
  // returned in ASCENDING order by symmetricEigenvaluesEigenvectors.
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;
  strain.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  // Liu et al. 2024 CMAME eq. (30): n_F = e_1 = eigenvector of the largest
  // (most-tensile) principal strain. Ascending order => column(2). The
  // eigenvectors are unit-norm by construction, so n_F is a valid unit normal.
  eps1 = eigval[2]; // maximum principal strain eps_1 (caller gates on eps_1 > 0)
  return eigvec.column(2);
}

void
NDSmallDeformationIsotropicElasticity::computeGDerivatives()
{
  if (_model_type == "AT2"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "AT1"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "PF_CZM"){

    // Reference: Gupta et al. (2022) An adaptive mesh refinement algorithm for phase-field fracture  models: Application to brittle, cohesive, and dynamic fracture
    // Get the parameters
    // a1, a2, a3, p, eta
    // read in the real properties on‐the‐fly
    const Real a1 = (*_a1_prop)[_qp];
    const Real a2 = (*_a2_prop)[_qp];
    const Real a3 = (*_a3_prop)[_qp];
    const Real p = (*_p_prop)[_qp];
    const Real d = _d[_qp];
    const Real eta = _eta;

    // degradation function
    _g[_qp] = std::pow((1-d),p)/(std::pow(1-d,p)+a1*d*(1+a2*d+a2*a3*std::pow(d,2)))*(1-_eta)+_eta;

    // here we break down _g into two parts: D = U + V
    // and copmute the derivatives separately, then combine them
    // U = (1-d)^p
    // V = a1 * (d + a2*d^2 + a2*a3*d^3)
    // D = (1-d)^p + a1 * (d + a2*d^2 + a2*a3*d^3)
    // Derivative of the degradation function w/r/t damage
    Real U   = std::pow(1-d, p);
    Real Up  = -p * std::pow(1-d, p-1);
    Real Up2 =  p*(p-1) * std::pow(1-d, p-2);

    Real V   = a1 * (d + a2*d*d + a2*a3*d*d*d);
    Real Vp  = a1 * (1 + 2*a2*d     + 3*a2*a3*d*d);
    Real Vpp = a1 * (    2*a2       + 6*a2*a3*d     );

    Real D   = U + V;
    Real Dp  = Up + Vp;
    Real Dpp = Up2 + Vpp;

    // first derivative g'
    Real N1  = Up*D - U*Dp;                       // numerator for g0'
    Real g0p = N1/(D*D);                          // g0'
    Real dg  = g0p * (1-eta);                     // g'

    // second derivative g''
    Real N2  = Up2*D  - U*Dpp;                    // numerator for N'
    Real g0pp = (N2*D - 2*N1*Dp)/(D*D*D);         // g0''
    Real d2g  = g0pp * (1-eta);                   // g''

    // store
    _dg_dd[_qp]    = dg;
    _d2g_dd2[_qp]  = d2g;
  }
  else
    mooseError("Unknown model type: " + _model_type);
}

void
NDSmallDeformationIsotropicElasticity::computeCrackStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _elastic_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}

void
NDSmallDeformationIsotropicElasticity::updatePermeabilityForCracking()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix (used by legacy exponential/Darcy branches)
  const RankTwoTensor & R = _crack_rotation[_qp];

  //Compute the intrinsic permeability
  RankTwoTensor perm_intrinsic = _intrinsic_permeability * RankTwoTensor::Identity();

  // Dispatch on the canonical permeability model enum.
  if (_permeability_model == PermeabilityModel::exponential)
  {
    // Legacy exponential permeability model (unchanged)
    RankTwoTensor effective_perm_new = perm_intrinsic * std::exp(_d[_qp] * _coeff_b);
    effective_perm_new.rotate(R);
    _effective_perm[_qp] = effective_perm_new;
  }
  else if (_permeability_model == PermeabilityModel::darcy_poiseuille)
  {
    // Legacy Darcy-Poiseuille model (unchanged): w = d * wc
    Real w = _d[_qp] * _wc;
    RankTwoTensor kf = std::pow(w, 2) / (12.0) * RankTwoTensor::Identity();
    RankTwoTensor effective_perm_new =
        perm_intrinsic + std::pow(_d[_qp], _perm_exponent) * (kf - perm_intrinsic);
    effective_perm_new.rotate(R);
    _effective_perm[_qp] = effective_perm_new;
  }
  else if (_permeability_model == PermeabilityModel::normal_strain)
  {
    // Heider 2021 eqs. 46-48: normal-strain-driven aperture with tangential
    // projector. K = K_poro + (d^b) * K_frac, where
    //   K_frac = (w_h^2 / 12) * (I - n_d (x) n_d)     (anisotropic), or
    //   K_frac = (w_h^2 / 12) * I                      (isotropic fallback).
    // w_h = f_c * w_c * chi_d,
    // w_c = h_c * |1 + n_d . eps . n_d|,
    // chi_d = H(d - d_threshold)   (strict: H(0) = 0),
    // n_d is either grad(d)/|grad(d)| or the most-tensile principal eigenvector.

    // (a) Determine crack normal n_d. NOTE: with regularize_crack_normal = true
    //     this is grad(d)/(|grad(d)| + eps), whose magnitude is < 1 (it tends to
    //     0 at the crack core); it is a true unit vector only on the legacy
    //     hard-cutoff and principal-strain paths.
    RealVectorValue n_d;
    bool have_normal = false;

    if (_normal_source == CrackNormalSource::damage_gradient)
    {
      const Real gnorm = _grad_d[_qp].norm();
      if (_regularize_crack_normal)
      {
        // Regularized crack normal n_d = grad(d) / (|grad(d)| + eps). As
        // |grad(d)| -> 0 at the fully-damaged crack core (d -> 1), n_d -> 0,
        // so the tangential projector (I - n_d (x) n_d) below tends to I and
        // the fracture permeability becomes isotropic there, instead of
        // falling back to the matrix permeability k0*I. eps > 0 (enforced by
        // the range check on crack_normal_regularization) guarantees no
        // division by zero, so have_normal is always true on this path.
        //
        // Note both d -> 0 (undamaged) and d -> 1 (crack core) have
        // |grad(d)| -> 0 and hence n_d -> 0 here, but the two are SEPARATED
        // downstream by the chi_d / d>0 gate: at d -> 0 the Heaviside gate
        // chi_d = 0 routes the point to matrix perm k0*I (undamaged
        // formulation unchanged), while at d -> 1 chi_d = 1 keeps the
        // isotropic fracture perm. So this regularized normal only changes
        // behavior at the damaged core, never in the undamaged bulk.
        n_d = _grad_d[_qp] / (gnorm + _crack_normal_reg_eps);
        have_normal = true;
      }
      else if (gnorm > _grad_d_tol)
      {
        n_d = _grad_d[_qp] / gnorm;
        have_normal = true;
      }
    }
    else // principal_strain: strain-based crack normal, Liu 2024 eqs. (29)-(30)
    {
      // Eq. (29)-(30): n_F = e_1 = eigenvector of the maximum principal strain.
      // Derive it from the model's mechanical strain `_total_strain` (the same
      // strain used for eps_nn below), NOT from _crack_rotation/_elastic_strain,
      // so the normal and the normal-strain aperture are computed from one
      // self-consistent strain tensor (identical in pure elasticity; the correct
      // total strain if plasticity is later attached). Eigenvectors are
      // unit-norm, so n_d is a valid unit normal everywhere (no |grad d|
      // division, no regularization needed -- contrast
      // crack_normal_source = damage_gradient).
      //
      // This re-eigendecomposes the TOTAL strain on purpose (self-consistent
      // with eps_nn); _crack_rotation (from _elastic_strain) is intentionally
      // not reused here -- it is dead work on this path but still feeds the
      // legacy exponential/Darcy branches.
      Real eps1 = 0.0;
      n_d = maxPrincipalStrainDirection((*_total_strain)[_qp], eps1);
      // Only enhance permeability when there is a tensile opening (eps_1 > 0).
      // Under in-plane compression (plane strain eps_zz = 0 is the max eigenvalue
      // => e_1 = e_z) or zero strain, eps_1 <= 0; have_normal is then false and
      // the chi_d / d-gate below routes the point to matrix perm k0*I, avoiding a
      // spurious aperture for a closed crack and a non-deterministic e_1 for a
      // degenerate (zero) strain.
      have_normal = (eps1 > 0.0);
    }

    const Real k0 = _intrinsic_permeability;
    const RankTwoTensor I2 = RankTwoTensor::Identity();
    const Real d = _d[_qp];

    // Heaviside gate chi_d (Heider eq. 46). Plan convention: H(0) = 1, i.e.
    // "1 if d >= threshold, else 0". The fracture contributes permeability
    // once damage *reaches* the threshold — required for the Phase-1
    // acceptance criterion at d = 0.5 (equal to the default threshold).
    const Real chi_d = (d >= _d_perm_threshold) ? 1.0 : 0.0;

    // Fall back to matrix perm if the normal is ill-defined, the Heaviside
    // gate is closed, or the damage is zero. This is also the d=0 vs d=1
    // separator for the regularized-normal path: with regularize_crack_normal
    // = true, have_normal is always true, so the fallback is driven purely by
    // damage -- chi_d == 0 (d < threshold) or d <= 0 sends the UNDAMAGED bulk
    // to k0*I, while the fully-damaged core (d >= threshold, chi_d = 1, d > 0)
    // proceeds to the isotropic fracture perm below.
    if (!have_normal || chi_d == 0.0 || d <= 0.0)
    {
      _effective_perm[_qp] = k0 * I2;
      return;
    }

    // (b) Normal strain eps_nn = n_d . eps^S . n_d (explicit double
    // contraction). Heider eq. (47) defines eps^S as the total linearized
    // (kinematic) strain, eq. (63). We bind _total_strain to the
    // "mechanical_strain" property of ComputeSmallStrain so this remains the
    // kinematic strain even if a plasticity model is later attached.
    Real eps_nn = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        eps_nn += n_d(i) * (*_total_strain)[_qp](i, j) * n_d(j);

    // Small-strain regime check (plan Risk §3). Warn once per run when the
    // normal strain enters the non-small regime; the actual clamp below is
    // what guarantees w_c = 0 for eps_nn < -1.
    if (eps_nn < -0.5)
      mooseDoOnce(mooseWarning(
          "Normal strain eps_nn = ", eps_nn,
          " < -0.5 violates small-strain assumptions in the normal-strain "
          "permeability model; further occurrences suppressed."));

    // Plan §Edge Case 3: clamp 1 + eps_nn to 0 for eps_nn < -1, so the
    // aperture does not spuriously re-grow under strong compression (Risk §3).
    const Real one_plus = std::max(1.0 + eps_nn, 0.0);

    // (c) Aperture (Heider eq. 47 with compression clamp):
    //     w_c = h_c * max(1 + eps_nn, 0)
    const Real h_c = getCharacteristicLength();

    // Guard against h_c <= 0: can happen at INITIAL with an ElementLengthAux
    // that has not yet been populated for a pre-damaged configuration. Fall
    // back to matrix perm rather than silently producing K_frac = 0.
    if (h_c <= 0.0)
    {
      _effective_perm[_qp] = k0 * I2;
      return;
    }

    const Real w_c = h_c * one_plus; // one_plus already >= 0 from clamp

    // (d) Roughness-corrected aperture (Heider eq. 46):
    //     w_h = max{ (f_c * w_c) * chi_d,    (open branch)
    //                (f_c * w_r) * chi_d }   (closed branch).
    // The residual aperture w_r floors w_h once damage exceeds the threshold,
    // so K_frac decays to (f_c*w_r)^2/12 under closure rather than to zero.
    // With _w_res = 0 (default), the closed branch is inactive and the
    // formula reduces to the open-only legacy form.
    const Real w_h = std::max(_fc * w_c * chi_d, _fc * _w_res * chi_d);
    const Real k_w = w_h * w_h / 12.0;

    // (e) Fracture permeability tensor.
    RankTwoTensor K_frac;
    if (_perm_anisotropic)
    {
      // Tangential projector K_frac = k_w * (I - n_d (x) n_d).
      RankTwoTensor n_outer_n;
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
          n_outer_n(i, j) = n_d(i) * n_d(j);
      K_frac = k_w * (I2 - n_outer_n);
    }
    else
    {
      K_frac = k_w * I2;
    }

    // (f) Damage-weighted total perm (Heider eq. 48).
    const Real weight = std::pow(d, _perm_exponent);
    _effective_perm[_qp] = k0 * I2 + weight * K_frac;
  }
  else
  {
    mooseError("Unknown permeability model type.");
  }

}

Real
NDSmallDeformationIsotropicElasticity::getCharacteristicLength() const
{
  switch (_lc_type)
  {
    case LcType::regularization_length:
      return (*_l_mat_prop)[_qp];
    case LcType::element_size:
      return (*_h_elem)[_qp];
    case LcType::constant:
      return _lc_const;
  }
  mooseError("Unknown characteristic_length_type.");
  return 0.0; // unreachable; silences -Wreturn-type
}

Real
NDSmallDeformationIsotropicElasticity::computeSpectralBulkModulus(const RankTwoTensor & strain)
{
  // Bulk modulus extracted from the degraded SPECTRAL elastic tangent C(d, eps)
  // by the volumetric contraction K = (1/9) I:C:I = (1/9) sum_{i,k} C_iikk.
  //
  // C        = C_intact + (g - 1) * C_pos
  // C_intact = K (I (x) I) + 2G (I4_sym - (1/3) I (x) I)
  // C_pos    = lambda * H(tr eps) * (I (x) I) + 2G * P_pos,  lambda = K - 2G/LIBMESH_DIM
  //
  // IMPORTANT: build the TRUE second-order identity dyad I (x) I = delta_ij delta_kl
  // via outerProduct. Do NOT use RankFourTensor(initIdentity): in this MOOSE build
  // that constructor is diagonal-only ((i,i,i,i)=1), so I:(initIdentity):I = 3 != 9
  // and the extracted bulk modulus would be wrong. We therefore do not reuse
  // computeJacobianSpectralDecomposition (which builds C with initIdentity).
  const Real g = _g[_qp]; // set by computeGDerivatives() before this call
  const Real K = _K[_qp];
  const Real G = _G[_qp];
  const Real lambda = K - 2.0 * G / LIBMESH_DIM; // matches the spectral stress/Jacobian

  const RankTwoTensor I2 = RankTwoTensor::Identity();
  const RankFourTensor IxI = I2.outerProduct(I2); // delta_ij delta_kl
  const RankFourTensor I4_sym(RankFourTensor::initIdentitySymmetricFour);

  // Intact isotropic tangent.
  const RankFourTensor C_intact = K * IxI + 2.0 * G * (I4_sym - IxI / 3.0);

  // Positive-projection part (spectral split), consistent with
  // computeStressSpectralDecomposition / computeJacobianSpectralDecomposition.
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  const RankFourTensor P_pos = strain.positiveProjectionEigenDecomposition(eigvals, eigvecs);
  const Real H_tr = (strain.trace() > 0.0) ? 1.0 : 0.0;
  const RankFourTensor C_pos = lambda * H_tr * IxI + 2.0 * G * P_pos;

  // Degraded spectral elastic tangent.
  const RankFourTensor C = C_intact + (g - 1.0) * C_pos;

  // Volumetric contraction K_eff = (1/9) I:C:I = (1/9) sum_{i,k} C_iikk.
  Real ICI = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      ICI += C(i, i, k, k);
  return ICI / 9.0;
}
