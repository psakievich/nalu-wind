/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "UnitTestAlgorithm.h"
#include "UnitTestKokkosUtils.h"
#include "UnitTestFieldUtils.h"
#include "UnitTestAlgorithmUtils.h"

#include "TurbViscKsgsAlgorithm.h"
#include "TurbViscSmagorinskyAlgorithm.h"
#include "TurbViscWaleAlgorithm.h"

TEST_F(TestTurbulenceAlgorithm, turbviscksgsalgorithm)
{
  sierra::nalu::Realm& realm = this->create_realm();

  fill_mesh_and_init_fields();

  // Execute
  sierra::nalu::TurbViscKsgsAlgorithm alg(realm, meshPart_);
  alg.execute();

  // Perform tests
  const double tol = 1e-14;
  double norm = field_norm(*tvisc_);
  const double gold_norm = 0.0285191520668428;
  EXPECT_NEAR(norm, gold_norm, tol);
}

TEST_F(TestTurbulenceAlgorithm, turbviscsmagorinskyalgorithm)
{
  sierra::nalu::Realm& realm = this->create_realm();

  fill_mesh_and_init_fields();

  // Execute
  sierra::nalu::TurbViscSmagorinskyAlgorithm alg(realm, meshPart_);
  alg.execute();

  // Perform tests
  const double tol = 1e-14;
  double norm = field_norm(*tvisc_);
  const double gold_norm = 0.0015635636790984;
  EXPECT_NEAR(norm, gold_norm, tol);
}

TEST_F(TestTurbulenceAlgorithm, turbviscwalealgorithm)
{
  sierra::nalu::Realm& realm = this->create_realm();

  fill_mesh_and_init_fields();

  // Execute
  sierra::nalu::TurbViscWaleAlgorithm alg(realm, meshPart_);
  alg.execute();

  // Perform tests
  const double tol = 1e-14;
  double norm = field_norm(*tvisc_);
  const double gold_norm = 0.0094154596233012953;
  EXPECT_NEAR(norm, gold_norm, tol);
}
