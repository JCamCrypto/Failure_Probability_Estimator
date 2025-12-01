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

#include "../distributions.h"

#include "ml-kem.h"

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

  char aSaveFile[] = "saved/DistributionA.save";

  initCenteredBinomialDistribution(C_PARAM_ETA1,  pCBD1);
  initCenteredBinomialDistribution(C_PARAM_ETA2,  pCBD2);
  initCompressionErrorDistribution(C_PARAM_DU,    pCompU);
  initCompressionErrorDistribution(C_PARAM_DV,    pCompV);

  // Compute the distribution of the error E = <e|y> - <s|e1+eu> + e2 + ev
  //  * e, s and y are sampled from CBD(eta1)
  //  * e1, e2 are sampled from CBD(eta2)
  //  * eu is sampled from CompU
  //  * ev is sampled from CompV

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is empty

  multiplyDistributions(pCBD1, pCBD1, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  aSaveFile[18u]++;

  // k * n * Temp1 is the distribution of <e|y>
  // Temp2 is empty
  // Temp3 is empty

  addDistributions(pCBD2, pCompU, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // k * n * Temp1 is the distribution of <e|y>
  // Temp2 is the distribution of e1 + eu
  // Temp3 is empty

  multiplyDistributions(pCBD1, pTemp2, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // k * n * Temp1 is the distribution of <e|y>
  // Temp2 is empty
  // k * n * Temp3 is the distribution of <s|e1+eu>

  addDistributions(pTemp1, pTemp3, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // k * n * Temp2 is the distribution of <e|y> + <s|e1+eu> = <e|y> - <s|e1+eu>   (by symmetry of s)
  // Temp3 is empty

#if (C_PARAM_K == 2u)
  addDistributions(pTemp2, pTemp2, pTemp3);
#elif (C_PARAM_K == 3u)
  addDistributions(pTemp2, pTemp2, pTemp1);
  addDistributions(pTemp1, pTemp2, pTemp3);
#elif (C_PARAM_K == 4u)
  addDistributions(pTemp2, pTemp2, pTemp1);
  addDistributions(pTemp1, pTemp1, pTemp3);
#endif
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // Temp2 is empty
  // n * Temp3 is the distribution of <e|y> - <s|e1+eu>

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
  // Temp3 is the distribution of <e|y> - <s|e1+eu>

  addDistributions(pCBD2, pCompV, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u] = '0';

  // Temp1 is empty
  // Temp2 is the distribution of e2 + ev
  // Temp3 is the distribution of <e|y> - <s|e1+eu>

  addDistributions(pTemp2, pTemp3, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  
  // Temp1 is the distribution of the coefficients of E = <e|y> - <s|e1+eu> + e2 + ev
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