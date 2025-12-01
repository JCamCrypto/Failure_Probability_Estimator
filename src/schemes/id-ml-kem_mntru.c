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

#include "../distributions.h"

#include "id-ml-kem_mntru.h"

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

  char aSaveFile[] = "saved/DistributionA.save";

  initCenteredBinomialDistribution(C_PARAM_ETA1,  pCBD1);
  initCenteredBinomialDistribution(C_PARAM_ETA2,  pCBD2);
  initCenteredNormalDistribution(C_PARAM_SIGMA,   pNorm);
  initCompressionErrorDistribution(C_PARAM_DU,    pCompU);
  initCompressionErrorDistribution(C_PARAM_DV,    pCompV);

  // Compute the distribution of the error E = e*y - <s|e1+eu> + e2 + ev
  //  * e, s are sampled from N(sigma)
  //  * y is sampled from CBD(eta1)
  //  * e1, e2 are sampled from CBD(eta2)
  //  * eu is sampled from CompU
  //  * ev is sampled from CompV

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is empty

  multiplyDistributions(pNorm, pCBD1, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of e*y
  // Temp2 is empty
  // Temp3 is empty

  addDistributions(pCBD2, pCompU, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of e*y
  // Temp2 is the distribution of e1 + eu
  // Temp3 is empty

  multiplyDistributions(pNorm, pTemp2, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of e*y
  // Temp2 is empty
  // k * n * Temp3 is the distribution of <s|e1+eu>

  // k = 2   =>   Temp2 = Temp3 + Temp3
  addDistributions(pTemp3, pTemp3, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of e*y
  // n * Temp2 is the distribution of <s|e1+eu>
  // Temp3 is empty

  // By symmetry of s, -Temp2 = Temp2. Therefore, Temp3 = Temp1 - Temp2 = Temp1 + Temp2
  addDistributions(pTemp1, pTemp2, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // Temp2 is empty
  // n * Temp3 is the distribution of e*y - <s|e1+eu>

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
  // Temp3 is the distribution of e*y - <s|e1+eu>

  addDistributions(pCBD2, pCompV, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u] = '0';

  // Temp1 is empty
  // Temp2 is the distribution of e2 + ev
  // Temp3 is the distribution of e*y - <s|e1+eu>

  addDistributions(pTemp2, pTemp3, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  
  // Temp1 is the distribution of the coefficients of E = e*y - <s|e1+eu> + e2 + ev
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