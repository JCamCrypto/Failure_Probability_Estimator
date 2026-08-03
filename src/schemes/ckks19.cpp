/**************************************************************************************************/
/** \brief    Compute the failure probability in CKKS19
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

#include "ckks19.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_PARAM_N             512u
#define C_PARAM_Q          524288u
#define C_PARAM_D               3u
#define C_PARAM_ETA             1u
#define C_PARAM_COMP            3u
#define C_PARAM_SIGMA         128.0

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
/** \implements computeFailureProbabilityOfCKKS19
 *
 **************************************************************************************************/
double computeFailureProbabilityOfCKKS19
(
  void
)
{
  TPDistribution pComp  = allocateDistribution();
  TPDistribution pNorm  = allocateDistribution();
  TPDistribution pUnif  = allocateDistribution();

  TPDistribution pTemp1 = allocateDistribution();
  TPDistribution pTemp2 = allocateDistribution();
  TPDistribution pTemp3 = allocateDistribution();

  initCenteredNormalDistribution(C_PARAM_SIGMA,   pNorm);
  initCenteredUniformDistribution(C_PARAM_ETA,    pUnif);
  initCompressionErrorDistribution(C_PARAM_COMP,  pComp);

  // The final error polynomial is computed as:
  //    E = r*s[0] - e[1]*s[1] - ... - e[d-1]*s[d-1] + e[0] + c
  // where:
  //  * r, e[0], ..., e[d-1] are sampled from U(eta)
  //  * s[0], ..., s[d-1] are sampled from N(sigma)
  //  * c is a compression error, resulting from dropping low-order bits of c[0]
  // So, the distribution of each coefficient of E can be computed as:
  //    D = d * n * U(eta) * N(sigma) + U(eta) + Comp

  multiplyDistributions(pUnif, pNorm, pTemp1);
  freeDistribution(pNorm);

  // Temp1 = U(eta) * N(sigma)

  applyScalarProduct(C_PARAM_D * C_PARAM_N, pTemp1, pTemp2);

  // Temp2 = d * n * U(eta) * N(sigma)

  addDistributions(pUnif, pComp, pTemp3);
  freeDistribution(pUnif);
  freeDistribution(pComp);

  // Temp2 = d * n * U(eta) * N(sigma)
  // Temp3 = U(eta) + Comp

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