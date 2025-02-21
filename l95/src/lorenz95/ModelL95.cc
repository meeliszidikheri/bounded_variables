/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/config/Configuration.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
static oops::interface::ModelMaker<L95Traits, ModelL95> makermodel_("L95");

// -----------------------------------------------------------------------------

ModelL95::ModelL95(const Resolution & resol, const eckit::Configuration & config)
  : resol_(resol), f_(config.getDouble("f")),
    tstep_(util::Duration(config.getString("tstep"))),
    dt_(tstep_.toSeconds()/432000.0), vars_({"x"})
{
  // set default parameter values for aerosol model so that model does not fail
  // when parameters not specified in config file; otherwise use values in config file
  include_aerosol_ = false;
  if (config.has("include_aerosol")) include_aerosol_ = config.getBool("include_aerosol");
  dust_model_ = false;
  if (config.has("dust_model")) dust_model_ = config.getBool("dust_model");
  alpha_ = 1.0;
  if (config.has("alpha")) alpha_ = config.getDouble("alpha");
  S_ = 1.0;
  if (config.has("S")) S_ = config.getDouble("S");
  grid_spacing_ = 1.0;
  if (config.has("grid_spacing")) grid_spacing_ = config.getDouble("grid_spacing");
  exponent_ = 1.0;
  if (config.has("exponent")) exponent_ = config.getDouble("exponent");
  oops::Log::info() << *this << std::endl;
  oops::Log::trace() << "ModelL95::ModelL95 created" << std::endl;
}

// -----------------------------------------------------------------------------
ModelL95::~ModelL95()
{
  oops::Log::trace() << "ModelL95::~ModelL95 destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelL95::initialize(StateL95 &) const {}
void ModelL95::finalize(StateL95 &) const {}
// -----------------------------------------------------------------------------
void ModelL95::step(StateL95 & xx, const ModelBias & bias) const {
  ModelTrajectory traj(false);
  this->stepRK(xx.getField(), bias, traj);
  xx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------

void ModelL95::stepRK(FieldL95 & xx, const ModelBias & bias,
                      ModelTrajectory & traj) const {
  FieldL95 dx(xx, false);
  FieldL95 zz(xx, false);
  FieldL95 dz(xx, false);

  zz = xx;
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx = dz;

  zz = xx;
  zz.axpy(0.5, dz);
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx.axpy(2.0, dz);

  zz = xx;
  zz.axpy(0.5, dz);
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx.axpy(2.0, dz);

  zz = xx;
  zz += dz;
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx += dz;

  const double zt = 1.0/6.0;
  xx.axpy(zt, dx);
}

// -----------------------------------------------------------------------------

#ifdef __INTEL_COMPILER
#pragma optimize("", off)
#endif
void ModelL95::tendencies(const FieldL95 & xx, const double & bias,
                          FieldL95 & dx) const {
  const int nn = resol_.npoints();
  // if including aerosol then divide domain into half u-field, half c-field;
  // (therefore, need to specify twice the resolution in config file to maintain
  // same configuration of the u-field when no aerosol present)
  int nh = nn;
  if (include_aerosol_) {
      nh = nn / 2;
  }
  // intel 19 is doing some agressive optimization of this loop that
  // is modifying the solution.
  for (int jj = 0; jj < nh; ++jj) {
    int jm2 = jj - 2;
    int jm1 = jj - 1;
    int jp1 = jj + 1;
    int jp2 = jj + 2;
    if (jm2 < 0) jm2 += nh;
    if (jm1 < 0) jm1 += nh;
    if (jp1 >= nh) jp1 -= nh;
    if (jp2 >= nh) jp2 -= nh;
    const double dxdt = -xx[jm2] * xx[jm1] + xx[jm1] * xx[jp1] - xx[jj] + f_ - bias;
    dx[jj] = dt_ * dxdt;
    // start aerosol tendency
    if (include_aerosol_) {
        double cc_adv = 0.0;
        if (xx[jj] > 0) {
             cc_adv = ( xx[jj+nh] - xx[jm1+nh] ) / grid_spacing_;
        } else {
             cc_adv = (-xx[jj+nh] + xx[jp1+nh] ) / grid_spacing_;
        }
        // source term, S: either two-point source (at fixed locations), or
        // "dust model" with S = S_ * u^exponent dependent on u-field
        double S;
        if (dust_model_) {
              S = S_ * pow(xx[jj], exponent_);
        } else {
              // source locations hardcoded initially at one-quarter and
              // three-quarter positions (will not be exact if number of
              // grid points not a multiple of 4).
              int loc1 = nh / 4;
              int loc2 = (3 * nh) / 4;
              S = ((jj == loc1) || (jj == loc2)) ? S_ : 0.0;
        }
        const double dcdt = -xx[jj]*cc_adv - alpha_*xx[jj+nh] + S;
        dx[jj+nh] = dt_ * dcdt;
    } // end aerosol tendency
  }
}
#ifdef __INTEL_COMPILER
#pragma optimize("", on)
#endif

// -----------------------------------------------------------------------------

void ModelL95::print(std::ostream & os) const {
  os << "ModelL95: resol = " << resol_ << ", f = " << f_ << ", tstep = " << tstep_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
