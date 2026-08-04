/**************************************************************************************************/
/** \brief    Compute the failure probability in DLP14
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

#include "dlp14.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#if (C_VAL_SECURITY == 80u)
# define C_PARAM_N             512u
# define C_PARAM_Q         8399873u
# define C_PARAM_ETA             1u
# define C_PARAM_DU             24u
# define C_PARAM_DV              3u
# define C_PARAM_SIGMA        4605.0
#elif (C_VAL_SECURITY == 192u)
# define C_PARAM_N            1024u
# define C_PARAM_Q       134246401u
# define C_PARAM_ETA             1u
# define C_PARAM_DU             28u
# define C_PARAM_DV              3u
# define C_PARAM_SIGMA       26583.0
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
/** \implements computeFailureProbabilityOfDLP14
 *
 **************************************************************************************************/
double computeFailureProbabilityOfDLP14
(
  void
)
{
  TPDistribution pCompU = allocateDistribution();
  TPDistribution pCompV = allocateDistribution();
  TPDistribution pNorm  = allocateDistribution();
  TPDistribution pUnif  = allocateDistribution();

  TPDistribution pTemp1 = allocateDistribution();
  TPDistribution pTemp2 = allocateDistribution();
  TPDistribution pTemp3 = allocateDistribution();

  initCenteredNormalDistribution(C_PARAM_SIGMA, pNorm);
  initCenteredUniformDistribution(C_PARAM_ETA,  pUnif);
  initCompressionErrorDistribution(C_PARAM_DU,  pCompU);
  initCompressionErrorDistribution(C_PARAM_DV,  pCompV);

  // The final error polynomial is computed as:
  //    E = s1*r - s2*(e1 + eu) + e2 + ev
  // where:
  //  * s1, s2 are sampled from N(sigma)
  //  * r, e1, e2 are sampled from U(eta)
  //  * eu, ev are compression errors, resulting from dropping low-order bits
  // So, the distribution of each coefficient of E can be computed as:
  //    D = n * [N(sigma) * U(eta) + N(sigma) * [U(eta) + CompU]] + U(eta) + CompV

  multiplyDistributions(pNorm, pUnif, pTemp1);

  // Temp1 = N(sigma) * U(eta)

  addDistributions(pUnif, pCompU, pTemp2);
  freeDistribution(pCompU);

  // Temp1 = N(sigma) * U(eta)
  // Temp2 = U(eta) + CompU

  multiplyDistributions(pNorm, pTemp2, pTemp3);
  freeDistribution(pNorm);

  // Temp1 = N(sigma) * U(eta)
  // Temp3 = N(sigma) * [U(eta) + CompU]

  addDistributions(pTemp1, pTemp3, pTemp2);

  // Temp2 = N(sigma) * U(eta) + N(sigma) * [U(eta) + CompU]

  applyScalarProduct(C_PARAM_N, pTemp2, pTemp3);

  // Temp3 = n * [N(sigma) * U(eta) + N(sigma) * [U(eta) + CompU]]

  addDistributions(pUnif, pCompV, pTemp2);
  freeDistribution(pUnif);
  freeDistribution(pCompV);

  // Temp2 = U(eta) + CompV
  // Temp3 = n * [N(sigma) * U(eta) + N(sigma) * [U(eta) + CompU]]

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