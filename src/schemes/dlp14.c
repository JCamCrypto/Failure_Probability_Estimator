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

#include "../distributions.h"

#include "dlp14.h"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#if (C_VAL_SECURITY == 80u)
# define C_PARAM_N             512u
# define C_PARAM_Q         8399873u
# define C_PARAM_ETA             1u
# define C_PARAM_DV              3u
# define C_PARAM_SIGMA        4605.0
#elif (C_VAL_SECURITY == 192u)
# define C_PARAM_N            1024u
# define C_PARAM_Q       134246401u
# define C_PARAM_ETA             1u
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

  // Compute the distribution of the error E = r*e - e1*sk + e2 + ev
  //  * e, s are sampled from N(sigma)
  //  * r, e1, e2 are sampled from U(eta)
  //  * ev is sampled from CompV

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is empty

  multiplyDistributions(pUnif, pNorm, pTemp1);
  saveDistribution(pTemp1, aSaveFile);
  aSaveFile[18u]++;

  // n * Temp1 is the distribution of r*e, and of e1*sk
  // Temp2 is empty
  // Temp3 is empty

  // By symmetry of s, -Temp1 = Temp1. Therefore, Temp3 = Temp1 - Temp1 = Temp1 + Temp1
  addDistributions(pTemp1, pTemp1, pTemp3);
  saveDistribution(pTemp3, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // Temp2 is empty
  // n * Temp3 is the distribution of r*e - e1*sk

  // Since n is a power of 2, we use log2(n) doublings
  size_t i = 1u;
  while (C_PARAM_N > i)
  {
    addDistributions(pTemp3, pTemp3, pTemp2);
    saveDistribution(pTemp2, aSaveFile);
    aSaveFile[18u]++;
    i <<= 1u;
    if (C_PARAM_N > i) {
      addDistributions(pTemp2, pTemp2, pTemp3);
      saveDistribution(pTemp3, aSaveFile);
      aSaveFile[18u]++;
      i <<= 1u;
    }
    else
    {
      // After the loop, we expect the result to be in Temp3
      copyDistribution(pTemp2, pTemp3);
    }
  }

  // Temp1 is empty
  // Temp2 is empty
  // Temp3 is the distribution of r*e - e1*sk

  addDistributions(pUnif, pCompV, pTemp2);
  saveDistribution(pTemp2, aSaveFile);
  aSaveFile[18u]++;

  // Temp1 is empty
  // Temp2 is the distribution of e2 + ev
  // Temp3 is the distribution of r*e - e1*sk

  addDistributions(pTemp2, pTemp3, pTemp1);
  aSaveFile[18u] = '0';
  saveDistribution(pTemp1, aSaveFile);
  
  // Temp1 is the distribution of the coefficients of E = r*e - e1*sk + e2 + ev
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