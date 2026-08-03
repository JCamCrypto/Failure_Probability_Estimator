/**************************************************************************************************/
/** \brief    Compute the failure probability in ID-ML-KEM_MNTRU
 *
 *  \author   Julien CAM
 *
 *  \date     2025/09/22
 *
 *  \file
 **************************************************************************************************/

/* ---------------------------------------------------------------------------------------------- */
/* IMPORTS                                                                                        */
/* ---------------------------------------------------------------------------------------------- */

#include "../distributions.hpp"

#include "id-ml-kem_mntru.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_PARAM_N            1024u
#define C_PARAM_Q         8380417u
#define C_PARAM_K               2u
#define C_PARAM_ETA1            3u
#define C_PARAM_ETA2            2u
#define C_PARAM_DU             19u
#define C_PARAM_DV              2u
#define C_PARAM_SIGMA         325.0

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL VARIABLES                                                                                */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - PROTOTYPE                                                                    */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC VARIABLES                                                                               */
/* ---------------------------------------------------------------------------------------------- */

size_t gModulus = C_PARAM_Q;

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC FUNCTIONS - IMPLEMENTATION                                                              */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \implements computeFailureProbabilityOfIDMLKEM
 *
 **************************************************************************************************/
double computeFailureProbabilityOfIDMLKEM
(
  void
)
{
  TPDistribution pCBD1  = allocateDistribution();
  TPDistribution pCBD2  = allocateDistribution();
  TPDistribution pCompU = allocateDistribution();
  TPDistribution pCompV = allocateDistribution();
  TPDistribution pNorm  = allocateDistribution();

  TPDistribution pTemp1 = allocateDistribution();
  TPDistribution pTemp2 = allocateDistribution();
  TPDistribution pTemp3 = allocateDistribution();

  initCenteredBinomialDistribution(C_PARAM_ETA1,  pCBD1);
  initCenteredBinomialDistribution(C_PARAM_ETA2,  pCBD2);
  initCenteredNormalDistribution(C_PARAM_SIGMA,   pNorm);
  initCompressionErrorDistribution(C_PARAM_DU,    pCompU);
  initCompressionErrorDistribution(C_PARAM_DV,    pCompV);

  // The final error polynomial is computed as:
  //    E = e*y - <s|e1+eu> + e2 + ev
  // where:
  //  * e, s are sampled from N(sigma)
  //  * y is sampled from CBD(eta1)
  //  * e1, e2 are sampled from CBD(eta2)
  //  * eu is sampled from CompU
  //  * ev is sampled from CompV
  // So, the distribution of each coefficient of E can be computed as:
  //    D = n * [N(sigma) * CBD(eta1) + k * N(sigma) * [CBD(eta2) + CompU]] + CBD(eta2) + CompV

  multiplyDistributions(pNorm, pCBD1, pTemp1);
  freeDistribution(pCBD1);

  // Temp1 = N(sigma) * CBD(eta1)

  addDistributions(pCBD2, pCompU, pTemp2);
  freeDistribution(pCompU);

  // Temp1 = N(sigma) * CBD(eta1)
  // Temp2 = CBD(eta2) + CompU

  multiplyDistributions(pNorm, pTemp2, pTemp3);
  freeDistribution(pNorm);

  // Temp1 = N(sigma) * CBD(eta1)
  // Temp3 = N(sigma) * [CBD(eta2) + CompU]

  applyScalarProduct(C_PARAM_K, pTemp3, pTemp2);

  // Temp1 = N(sigma) * CBD(eta1)
  // Temp2 = k * N(sigma) * [CBD(eta2) + CompU]

  addDistributions(pTemp1, pTemp2, pTemp3);

  // Temp3 = N(sigma) * CBD(eta1) + k * N(sigma) * [CBD(eta2) + CompU]

  applyScalarProduct(C_PARAM_N, pTemp3, pTemp2);

  // Temp2 = n * [N(sigma) * CBD(eta1) + k * N(sigma) * [CBD(eta2) + CompU]]

  addDistributions(pCBD2, pCompV, pTemp3);
  freeDistribution(pCBD2);
  freeDistribution(pCompV);

  // Temp2 = n * [N(sigma) * CBD(eta1) + k * N(sigma) * [CBD(eta2) + CompU]]
  // Temp3 = CBD(eta2) + CompV

  addDistributions(pTemp2, pTemp3, pTemp1);
  freeDistribution(pTemp2);
  freeDistribution(pTemp3);

  // Temp1 = D

  saveDistribution(pTemp1, "FinalErrorDistribution.save");
  
  // The probability that one given coefficient of the error polynomial is not rounded to 0
  double failureProbability = computeRoundingToOneProbability(pTemp1);
  freeDistribution(pTemp1);
  // The probability that all the coefficients of the error polynomial are not rounded to 0
  return failureProbability * ((double) C_PARAM_N);
}

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - IMPLEMENTATION                                                               */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */