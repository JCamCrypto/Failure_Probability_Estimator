/**************************************************************************************************/
/** \brief    Compute the failure probability in ML-KEM
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

#include "ml-kem.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#if (C_VAL_SECURITY == 128u)
# define C_PARAM_N            256u
# define C_PARAM_Q           3329u
# define C_PARAM_K              2u
# define C_PARAM_ETA1           3u
# define C_PARAM_ETA2           2u
# define C_PARAM_DU            10u
# define C_PARAM_DV             4u
#elif (C_VAL_SECURITY == 192u)
# define C_PARAM_N            256u
# define C_PARAM_Q           3329u
# define C_PARAM_K              3u
# define C_PARAM_ETA1           2u
# define C_PARAM_ETA2           2u
# define C_PARAM_DU            10u
# define C_PARAM_DV             4u
#elif (C_VAL_SECURITY == 256u)
# define C_PARAM_N            256u
# define C_PARAM_Q           3329u
# define C_PARAM_K              4u
# define C_PARAM_ETA1           2u
# define C_PARAM_ETA2           2u
# define C_PARAM_DU            11u
# define C_PARAM_DV             5u
#endif

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
/** \implements computeFailureProbabilityOfMLKEM
 *
 **************************************************************************************************/
double computeFailureProbabilityOfMLKEM
(
  void
)
{
  TPDistribution pCBD1  = allocateDistribution();
  TPDistribution pCBD2  = allocateDistribution();
  TPDistribution pCompU = allocateDistribution();
  TPDistribution pCompV = allocateDistribution();

  TPDistribution pTemp1 = allocateDistribution();
  TPDistribution pTemp2 = allocateDistribution();
  TPDistribution pTemp3 = allocateDistribution();

  initCenteredBinomialDistribution(C_PARAM_ETA1,  pCBD1);
  initCenteredBinomialDistribution(C_PARAM_ETA2,  pCBD2);
  initCompressionErrorDistribution(C_PARAM_DU,    pCompU);
  initCompressionErrorDistribution(C_PARAM_DV,    pCompV);

  // The final error polynomial is computed as:
  //    E = <e|y> - <s|e1+eu> + e2 + ev
  // where:
  //  * e, s and y are sampled from CBD(eta1)
  //  * e1, e2 are sampled from CBD(eta2)
  //  * eu is sampled from CompU
  //  * ev is sampled from CompV
  // So, the distribution of each coefficient of E can be computed as:
  //    D = k * n * [CBD(eta1) * CBD(eta1) + CBD(eta1) * [CBD(eta2) + CompU]] + CBD(eta2) + CompV

  multiplyDistributions(pCBD1, pCBD1, pTemp1);

  // Temp1 = CBD(eta1) * CBD(eta1)

  addDistributions(pCBD2, pCompU, pTemp2);
  freeDistribution(pCompU);

  // Temp1 = CBD(eta1) * CBD(eta1)
  // Temp2 = CBD(eta2) + CompU

  multiplyDistributions(pCBD1, pTemp2, pTemp3);
  freeDistribution(pCBD1);

  // Temp1 = CBD(eta1) * CBD(eta1)
  // Temp3 = CBD(eta1) * [CBD(eta2) + CompU]

  addDistributions(pTemp1, pTemp3, pTemp2);

  // Temp2 = CBD(eta1) * CBD(eta1) + CBD(eta1) * [CBD(eta2) + CompU]

  applyScalarProduct(C_PARAM_K * C_PARAM_N, pTemp2, pTemp3);

  // Temp3 = k * n * [CBD(eta1) * CBD(eta1) + CBD(eta1) * [CBD(eta2) + CompU]]

  addDistributions(pCBD2, pCompV, pTemp2);
  freeDistribution(pCBD2);
  freeDistribution(pCompV);

  // Temp2 = CBD(eta2) + CompV
  // Temp3 = k * n * [CBD(eta1) * CBD(eta1) + CBD(eta1) * [CBD(eta2) + CompU]]

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