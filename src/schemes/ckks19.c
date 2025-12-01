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

#include "../distributions.h"

#include "ckks19.h"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_PARAM_N             512u
#define C_PARAM_Q          524288u
#define C_PARAM_K               2u
#define C_PARAM_ETA             1u
#define C_PARAM_DV              3u
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
  TPDistribution pCompV = allocateDistribution();
  TPDistribution pNorm  = allocateDistribution();
  TPDistribution pUnif  = allocateDistribution();

  TPDistribution pTemp1 = allocateDistribution();
  TPDistribution pTemp2 = allocateDistribution();
  TPDistribution pTemp3 = allocateDistribution();

  char aSaveFile[] = "saved/DistributionA.save";

  initCenteredNormalDistribution(C_PARAM_SIGMA,   pNorm);
  initCenteredUniformDistribution(C_PARAM_ETA,    pUnif);
  initCompressionErrorDistribution(C_PARAM_DV,    pCompV);

  // Compute the distribution of the error E = r*e - <sk|e1> + e2 + ev
  //  * e, sk are sampled from N(sigma)
  //  * r, e1, e2 are sampled from U(eta)
  //  * ev is sampled from CompV

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is empty

  multiplyDistributions(pUnif, pNorm, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of r*e
  // Temp2 is empty
  // Temp3 is empty

  multiplyDistributions(pNorm, pUnif, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of r*e
  // Temp2 is empty
  // k * n * Temp3 is the distribution of <sk|e1>

  // k = 2   =>   Temp2 = Temp3 + Temp3
  addDistributions(pTemp3, pTemp3, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of r*e
  // n * Temp2 is the distribution of <sk|e1>
  // Temp3 is empty

  // By symmetry of s, -Temp2 = Temp2. Therefore, Temp3 = Temp1 - Temp2 = Temp1 + Temp2
  addDistributions(pTemp1, pTemp2, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // Temp2 is empty
  // n * Temp3 is the distribution of r*e - <sk|e1>

  // Since n is a power of 4, we use log4(n) pairs of doublings
  for (size_t i = 1u; i < C_PARAM_N; i <<= 2u)
  {
    addDistributions(pTemp3, pTemp3, pTemp2);
    saveDistribution(pTemp2, aSaveFile);
    aSaveFile[18u]++;
    addDistributions(pTemp2, pTemp2, pTemp3);
    saveDistribution(pTemp3, aSaveFile);
    aSaveFile[18u]++;
  }

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is the distribution of r*e - <sk|e1>

  addDistributions(pUnif, pCompV, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u] = '0';

  // Temp1 is empty
  // Temp2 is the distribution of e2 + ev
  // Temp3 is the distribution of r*e - <sk|e1>

  addDistributions(pTemp2, pTemp3, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  
  // Temp1 is the distribution of the coefficients of E = r*e - <sk|e1> + e2 + ev
  // Temp2 is empty
  // Temp3 is empty

  // The probability that one given coefficient of the error polynomial is not rounded to 0
  double failureProbability = computeRoundingToOneProbability(pTemp1);
  // The probability that all the coefficients of the error polynomial are not rounded to 0
  failureProbability *= (double) C_PARAM_N;

  return failureProbability;
}

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - IMPLEMENTATION                                                               */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */