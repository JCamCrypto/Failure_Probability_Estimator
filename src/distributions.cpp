/**************************************************************************************************/
/** \brief    Computations on symmetric distributions
 *
 *  \author   Julien CAM
 *
 *  \date     2026/07/22
 *
 *  \file
 **************************************************************************************************/

/* ---------------------------------------------------------------------------------------------- */
/* IMPORTS                                                                                        */
/* ---------------------------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "distributions.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_LEN_ARRAYS  ((gModulus / 2u) + 2u)    //!< Length of the arrays storing distributions
#define C_OFS_BOUND   ((gModulus / 2u) + 1u)    //!< Offset to the bound on the distribution

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL VARIABLES                                                                                */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - PROTOTYPE                                                                    */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC VARIABLES                                                                               */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC FUNCTIONS - IMPLEMENTATION                                                              */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \implements allocateDistribution
 * 
 **************************************************************************************************/
TPDistribution allocateDistribution
(
  void
)
{
  return (TPDistribution) malloc(C_LEN_ARRAYS * sizeof(double));
}

/**************************************************************************************************/
/** \implements applyScalarProduct
 * 
 **************************************************************************************************/
void applyScalarProduct
(
  const size_t xScalar,
  const TPDistribution xpDistribution,
  TPDistribution xpResult
)
{
  copyDistribution(xpDistribution, xpResult);

  TPDistribution pTemp = allocateDistribution();

  for (size_t i = (size_t) log2(xScalar); 0u < i; i--)
  {
    addDistributions(xpResult, xpResult, pTemp);
    copyDistribution(pTemp, xpResult);

    if (1u == ((xScalar >> (i - 1u)) & 1u))
    {
      addDistributions(xpResult, xpDistribution, pTemp);
      copyDistribution(pTemp, xpResult);
    }
  }
}

/**************************************************************************************************/
/** \implements computeRoundingToOneProbability
 * 
 **************************************************************************************************/
double computeRoundingToOneProbability
(
  const TPDistribution xpDistribution
)
{
  double result = 0.0;

  // Sum first the lowest probabilities to increase the precision
  // Because gModulus is prime,
  // * q is not a multiple of 4, so: [q/4, 3q/4] = [floor(q/4)+1, q-floor(q/4)-1]
  // * q is not a multiple of 2, so: floor(q/2) =/= q-floor(q/2)
  for (size_t x = (size_t) xpDistribution[C_OFS_BOUND]; gModulus / 4u < x; x--)
  {
    result += 2.0 * xpDistribution[x];
  }

  return result;
}

/**************************************************************************************************/
/** \implements copyDistribution
 * 
 **************************************************************************************************/
void copyDistribution
(
  TPDistribution xpInput,
  TPDistribution xpOutput
)
{
  for (size_t k = 0u; k <= xpInput[C_OFS_BOUND]; k++)
  {
    xpOutput[k] = xpInput[k];
  }

  xpOutput[C_OFS_BOUND] = xpInput[C_OFS_BOUND];
}

/**************************************************************************************************/
/** \implements freeDistribution
 * 
 **************************************************************************************************/
void freeDistribution
(
  const TPDistribution xpDistribution
)
{
  free(xpDistribution);
}

/**************************************************************************************************/
/** \implements initCenteredBinomialDistribution
 * 
 **************************************************************************************************/
void initCenteredBinomialDistribution
(
  const size_t xEta,
  TPDistribution xpDistribution
)
{
  size_t precomputed1 = 2u * xEta;
  double precomputed2 = 1.0 / ((double) (1u << precomputed1));

  // For k between 0 and xEta, the probability of k is [k+xEta among 2*xEta]/2^(2*xEta)
  for (size_t k = 0u; xEta >= k; k++)
  {
    // Initialize the probability to 1/2^(2*xEta)
    xpDistribution[k] = precomputed2;

    // Multiply the probability by [k+xEta among 2*xEta]
    for (size_t i = 0u; xEta + k > i; i++)
    {
      xpDistribution[k] *= ((double) (precomputed1 - i)) / ((double) (i + 1u));
    }
  }

  // The distribution is only zero in [-xEta, xEta]
  xpDistribution[C_OFS_BOUND] = (double) xEta;
}

/**************************************************************************************************/
/** \implements initCenteredNormalDistribution
 * 
 **************************************************************************************************/
void initCenteredNormalDistribution
(
  const double xStdDev,
  TPDistribution xpDistribution
)
{
  size_t x;
  double weight;
  // The probability of 0 is proportional to: exp(-(0 / xStdDev)^2) / 2) = 1
  xpDistribution[0u] = 1.0;
  double totalWeight = 1.0;
  // If no value has a zero probability, the bound is q/2
  xpDistribution[C_OFS_BOUND] = (double) (gModulus / 2u);

  for (x = 1u; gModulus / 2u >= x; x++)
  {
    // The probability of x is proportional to: exp(-(x / xStdDev)^2) / 2)
    weight = ((double) x) / xStdDev;
    weight = exp(-(weight * weight) / 2.0);

    // If the probability of x is zero, then all remaining probabilities will also be 0
    if (0.0 == weight)
    {
      xpDistribution[C_OFS_BOUND] = (double) (x - 1u);
      break;
    }

    xpDistribution[x] = weight;
    totalWeight += 2.0 * weight;
  }

  // Normalize the distribution, and note the last value with non-zero probability
  for (x = 0u; x <= (size_t) xpDistribution[C_OFS_BOUND]; x++)
  {
    weight = xpDistribution[x] / totalWeight;

    // If the probability of x is zero, then all remaining probabilities will also be 0
    if (0.0 == weight)
    {
      xpDistribution[C_OFS_BOUND] = (double) (x - 1u);
      break;
    }

    xpDistribution[x] = weight;
  }
}

/**************************************************************************************************/
/** \implements initCenteredUniformDistribution
 * 
 **************************************************************************************************/
void initCenteredUniformDistribution
(
  const size_t xBound,
  TPDistribution xpDistribution
)
{
  double probability = 1.0 / (2.0 * ((double) xBound) + 1.0);

  for (size_t x = 0u; xBound >= x; x++)
  {
    xpDistribution[x] = probability;
  }

  xpDistribution[C_OFS_BOUND] = (double) xBound;
}

/**************************************************************************************************/
/** \implements initCompressionErrorDistribution
 * 
 **************************************************************************************************/
void initCompressionErrorDistribution
(
  const size_t xNbCompressedBits,
  TPDistribution xpDistribution
)
{
  size_t x;
  size_t compressed;
  size_t decompressed;
  size_t error;

  // We don't know a bound on the support of the distribution yet
  xpDistribution[C_OFS_BOUND] = 0.0;

  // Initialize the distribution to 0
  for (x = 0u; gModulus / 2u >= x; x++)
  {
    xpDistribution[x] = 0.0;
  }

  // Compress and decompress all the possible values and count the errors
  for (x = 0u; gModulus > x; x++)
  {
    // Compress_q(x, d) = round((2^d / q) * x)
    // round(z) = floor(z + 1/2)
    // So, Compress_q(x, d) = ((x << d) + (q/2)) // q
    compressed = ((x << xNbCompressedBits) + (gModulus >> 1u)) / gModulus;
    // Decompress_q(y, d) = round((q / 2^d) * y)
    // round(z) = floor(z + 1/2)
    // So, Decompress_q(y, d) = ((y * q) + (1 << (d-1))) >> d
    decompressed = ((compressed * gModulus) + (1u << (xNbCompressedBits - 1u)))
                   >> xNbCompressedBits;
    
    error = (decompressed > x) ? (decompressed - x) : (x - decompressed);
    xpDistribution[error]++;
  }

  // Normalize the distribution, and note the last value with non-zero probability
  xpDistribution[0] /= (double) gModulus;

  for (x = 1u; gModulus / 2u >= x; x++)
  {
    xpDistribution[x] /= (double) (2u * gModulus);

    if ((0.0 == xpDistribution[x]) && (0.0 == xpDistribution[C_OFS_BOUND]))
    {
      xpDistribution[C_OFS_BOUND] = (double) (x - 1u);
    }
  }
}

/**************************************************************************************************/
/** \implements loadDistribution
 * 
 **************************************************************************************************/
void loadDistribution
(
  const char* xpFilename,
  TPDistribution xpDistribution
)
{
  FILE* pFile = fopen(xpFilename, "rb");

  size_t nbBitsRead = fread(xpDistribution, sizeof(double), C_LEN_ARRAYS, pFile);

  fclose(pFile);

  if (C_LEN_ARRAYS != nbBitsRead)
  {
    printf("File '%s' is malformed: read %zu values instead of %zu.", xpFilename, nbBitsRead, C_LEN_ARRAYS);
    exit(EXIT_FAILURE);
  }
}

/**************************************************************************************************/
/** \implements printDistribution
 * 
 **************************************************************************************************/
void printDistribution
(
  const TPDistribution xpDistribution
)
{
  printf("Bound: %zu\n", (size_t) xpDistribution[C_OFS_BOUND]);

  printf("------------------------------\n");
  printf("|  Value  |    Probability   |\n");
  printf("------------------------------\n");

  for (size_t value = 0u; gModulus > value; value++)
  {
    size_t index = (gModulus / 2u >= value) ? value : (gModulus - value);
    
    if (index <= xpDistribution[C_OFS_BOUND])
    {
      double probability = xpDistribution[index];

      if (1e-10 <= probability)
      {
        printf("| %7zu | %.14f |\n", value, probability);
      }
      else if (1e-99 <= probability)
      {
        printf("| %7zu | %15.9e  |\n", value, probability);
      }
      else
      {
        printf("| %7zu | %16.9e |\n", value, probability);
      }
    }
  }

  printf("------------------------------\n");
}

/**************************************************************************************************/
/** \implements saveDistribution
 * 
 **************************************************************************************************/
void saveDistribution
(
  const TPDistribution xpDistribution,
  const char* filename
)
{
  FILE* pFile = fopen(filename, "wb");

  fwrite(xpDistribution, sizeof(double), C_LEN_ARRAYS, pFile);

  fclose(pFile);
}

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - IMPLEMENTATION                                                               */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */