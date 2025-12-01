/**************************************************************************************************/
/** \brief    Computations on symmetric distributions
 * 
 *  \author   Julien CAM
 * 
 *  \date     2025/03/26
 *
 *  \file
 **************************************************************************************************/

/* ---------------------------------------------------------------------------------------------- */
/* IMPORTS                                                                                        */
/* ---------------------------------------------------------------------------------------------- */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "distributions.h"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_LEN_ARRAYS  ((gModulus / 2u) + 2u)    //!< Length of the arrays storing distributions
#define C_LEN_THREADS 96u                       //!< Number of threads for parallel operations
#define C_OFS_BOUND   ((gModulus / 2u) + 1u)    //!< Offset to the bound on the distribution

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL VARIABLES                                                                                */
/* ---------------------------------------------------------------------------------------------- */

// The pointers to the input/output distributions of a parallel operation must be global,
// so that every thread can access them without passing them as arguments
static TPDistribution gaDistribution1;
static TPDistribution gaDistribution2;
static TPDistribution gaResults[C_LEN_THREADS];
static TPDistribution gaResult;

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - PROTOTYPE                                                                    */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \brief  Actual implementation of the addition of two symmetric distributions
 *          This function is executed by several threads in parallel,
 *          each of them being in charge of computing a part of the result distribution
 *
 * \param[in]       xpThreadId              Pointer to the number of the current thread
 *                                          Must a size_t, casted to a void*
 *
 **************************************************************************************************/
static void* addDistributionsInternal
(
  void* xpThreadId
);

/**************************************************************************************************/
/** \brief  Actual implementation of the product of two symmetric distributions
 *          This function is executed by several threads in parallel,
 *          each of them being in charge of computing a part of the result distribution
 *
 * \param[in]       xpThreadId              Pointer to the number of the current thread
 *                                          Must a size_t, casted to a void*
 *
 **************************************************************************************************/
static void* multiplyDistributionsInternal
(
  void* xpThreadId
);

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC VARIABLES                                                                               */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* PUBLIC FUNCTIONS - IMPLEMENTATION                                                              */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \implements addDistributions
 * 
 **************************************************************************************************/
void addDistributions
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  TPDistribution xpResult
)
{
  int status;
  pthread_t threads[C_LEN_THREADS];
  size_t t;

  // If |X| <= Bx and |Y| <= By, then |Z| <= Bx + By
  // Moreover, we always have |Z| <= q/2
  if (gModulus / 2u < xpDistribution1[C_OFS_BOUND] + xpDistribution2[C_OFS_BOUND])
  {
    xpResult[C_OFS_BOUND] = (double) (gModulus / 2u);
  }
  else
  {
    xpResult[C_OFS_BOUND] = xpDistribution1[C_OFS_BOUND] + xpDistribution2[C_OFS_BOUND];
  }

  gaDistribution1 = (TPDistribution) xpDistribution1;
  gaDistribution2 = (TPDistribution) xpDistribution2;
  gaResult = xpResult;

  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    status = pthread_create(&threads[t], NULL, addDistributionsInternal, (void*) t);

    if (0 != status)
    {
      fprintf(stderr, "Error: pthread_create failed with status %d\n", status);
      exit(EXIT_FAILURE);
    }
  }

  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    pthread_join(threads[t], NULL);
  }
}

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
/** \implements computeRoundingToOneProbability
 * 
 **************************************************************************************************/
double computeRoundingToOneProbability
(
  const TPDistribution xpDistribution
)
{
  size_t x;
  double result = 0.0;

  // Sum first the lowest probabilities to increase the precision
  // Because gModulus is prime,
  // * q is not a multiple of 4, so: [q/4, 3q/4] = [floor(q/4)+1, q-floor(q/4)-1]
  // * q is not a multiple of 2, so: floor(q/2) =/= q-floor(q/2)
  for (x = gModulus / 2u; gModulus / 4u < x; x--)
  {
    result += 2.0 * xpDistribution[x];
  }

  return result;
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
  size_t k;
  size_t precomputed1 = 2u * xEta;
  double precomputed2 = 1.0 / ((double) (1u << precomputed1));

  // For k between 0 and xEta, the probability of k is [k+xEta among 2*xEta]/2^(2*xEta)
  for (k = 0u; xEta >= k; k++)
  {
    // Initialize the probability to 1/2^(2*xEta)
    xpDistribution[k] = precomputed2;

    // Multiply the probability by [k+xEta among 2*xEta]
    for (size_t i = 0u; xEta + k > i; i++)
    {
      xpDistribution[k] *= ((double) (precomputed1 - i)) / ((double) (i + 1u));
    }
  }

  // Outside [-xEta, xEta], the probability of k is 0
  for (; gModulus / 2u >= k; k++)
  {
    xpDistribution[k] = 0.0;
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
  double totalWeight = 1.0;
  xpDistribution[0u] = 1.0;
  // We don't know a bound on the support of the distribution yet
  xpDistribution[C_OFS_BOUND] = 0.0;

  // The probability of x is proportional to: exp(-(x / xStdDev)^2) / 2)
  for (x = 1u; gModulus / 2u >= x; x++)
  {
    weight = ((double)x) / xStdDev;
    weight = exp(-(weight * weight) / 2.0);
    xpDistribution[x] = weight;
    totalWeight += 2.0 * weight;
  }

  // Normalize the distribution, and note the last value with non-zero probability
  for (x = 0u; gModulus / 2u >= x; x++)
  {
    xpDistribution[x] /= totalWeight;

    if ((0.0 == xpDistribution[x]) && (0.0 == xpDistribution[C_OFS_BOUND]))
    {
      xpDistribution[C_OFS_BOUND] = (double) (x - 1u);
    }
  }

  // If no value has a non-zero probability, the bound is q/2
  if (0.0 == xpDistribution[C_OFS_BOUND])
  {
    xpDistribution[C_OFS_BOUND] = (double) (gModulus / 2u);
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
  size_t x = 0u;
  double probability = 1.0 / (2.0 * ((double)xBound) + 1.0);

  for (; xBound >= x; x++)
  {
    xpDistribution[x] = probability;
  }

  for (; gModulus / 2u >= x; x++)
  {
    xpDistribution[x] = 0.0;
  }

  xpDistribution[C_OFS_BOUND] = xBound;
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
  const char* filename,
  TPDistribution xpDistribution
)
{
  FILE* pFile = fopen(filename, "rb");

  fread(xpDistribution, sizeof(double), C_LEN_ARRAYS, pFile);

  fclose(pFile);
}

/**************************************************************************************************/
/** \implements multiplyDistribution
 * 
 **************************************************************************************************/
void multiplyDistributions
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  TPDistribution xpResult
)
{
  int status;
  size_t x;
  size_t t;
  pthread_t aThreads[C_LEN_THREADS];

  gaDistribution1 = (TPDistribution) xpDistribution1;
  gaDistribution2 = (TPDistribution) xpDistribution2;

  // Each thread will output to a partial result to its own dedicated buffer
  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    gaResults[t] = allocateDistribution();
  }

  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    status = pthread_create(&aThreads[t], NULL, multiplyDistributionsInternal, (void*) t);

    if (0 != status)
    {
      fprintf(stderr, "Error: pthread_create failed with status %d\n", status);
      exit(EXIT_FAILURE);
    }
  }

  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    pthread_join(aThreads[t], NULL);
  }

  for (x = 1u; C_LEN_ARRAYS > x; x++)
  {
    xpResult[x] = 0.0;
  }

  // If |X| < Bx and |Y| < By, then |Z| < Bx * By
  // Moreover, we always have |Z| <= q/2
  if (gModulus / 2u < xpDistribution1[C_OFS_BOUND] * xpDistribution2[C_OFS_BOUND])
  {
    xpResult[C_OFS_BOUND] = (double) (gModulus / 2u);
  }
  else
  {
    xpResult[C_OFS_BOUND] = xpDistribution1[C_OFS_BOUND] * xpDistribution2[C_OFS_BOUND];
  }

  // P(Z=0) = P[(X=0)∪(Y=0)]
  //        = P(X=0) + P(Y=0) - P[(X=0)∩(Y=0)]
  //        = P(X=0) + P(Y=0) - P(X=0) * P(Y=0)
  //        = P(X=0) * [1 - P(Y=0)] + P(Y=0)
  xpResult[0u] = xpDistribution1[0u] * (1.0 - xpDistribution2[0u]) + xpDistribution2[0u];

  // Sum the partial results computed by the threads
  for (t = 0u; C_LEN_THREADS > t; t++)
  {
    for (x = 1u; x <= xpResult[C_OFS_BOUND]; x++)
    {
      xpResult[x] += gaResults[t][x];
    }
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
  double currentValue;
  size_t x;

  printf("Bound: %zu\n", (size_t) xpDistribution[C_OFS_BOUND]);

  printf("------------------------------\n");
  printf("|  Value  |    Probability   |\n");
  printf("------------------------------\n");

  for (x = 0u; gModulus / 2u >= x; x++)
  {
    currentValue = xpDistribution[x];
    
    if (0.0 != currentValue)
    {
      if (1e-10 <= currentValue)
      {
        printf("| %7zu |   %.10f   |\n", x, currentValue);
      }
      else if (1e-99 <= currentValue)
      {
        printf("| %7zu |   %11.5e    |\n", x, currentValue);
      }
      else
      {
        printf("| %7zu |   %12.5e   |\n", x, currentValue);
      }
    }
  }

  for (x = gModulus / 2u; 0u < x; x--)
  {
    currentValue = xpDistribution[x];
    
    if (0.0 != currentValue)
    {
      if (1e-10 <= currentValue)
      {
        printf("| %7zu |   %.10f   |\n", gModulus - x, currentValue);
      }
      else if (1e-99 <= currentValue)
      {
        printf("| %7zu |   %11.5e    |\n", gModulus - x, currentValue);
      }
      else
      {
        printf("| %7zu |   %12.5e   |\n", gModulus - x, currentValue);
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

/**************************************************************************************************/
/** \implements addDistributionsInternal
 * 
 **************************************************************************************************/
static void* addDistributionsInternal
(
  void* xpThreadId
)
{
  size_t threadId = (size_t) xpThreadId;
  size_t x;
  size_t y;
  size_t z;
  size_t bound;
  size_t boundOnX = (size_t) gaDistribution1[C_OFS_BOUND];
  size_t boundOnY = (size_t) gaDistribution2[C_OFS_BOUND];
  size_t boundOnZ = (size_t) gaResult[C_OFS_BOUND];
  double pz;

  for (z = threadId; z <= boundOnZ; z += C_LEN_THREADS)
  {
    pz = 0.0;
    
    // x + y = z        with: | x <= q/2, so P(X=x) is stored in distX[x]
    //                        | y <= q/2, so P(Y=y) is stored in distY[y]
    y = (z < boundOnY) ? z : boundOnY;
    x = z - y;
    bound = (z < boundOnX) ? z : boundOnX;

    for (; x <= bound; x++, y--)
    {
      pz += gaDistribution1[x] * gaDistribution2[y];
    }

    // x + y = z + q    with: | x <= q/2, so P(X=x) is stored in distX[x]
    //                        | y >  q/2, so P(Y=y) is stored in distY[q-y]
    // We apply the change of variable y -> q - y
    x = z + 1u;
    y = 1u;

    for (; (x <= boundOnX) && (y <= boundOnY); x++, y++)
    {
      pz += gaDistribution1[x] * gaDistribution2[y];
    }

    // x + y = z + q    with: | x > q/2, so P(X=x) is stored in distX[q-x]
    //                        | y > q/2, so P(Y=y) is stored in distY[q-y]
    // We apply the change of variable (x, y) -> (q - x, q - y)
    x = boundOnX;
    y = gModulus - z - boundOnX;

    for (; y <= boundOnY; x--, y++)
    {
      pz += gaDistribution1[x] * gaDistribution2[y];
    }

    // x + y = z + q    with: | x >  q/2, so P(X=x) is stored in distX[q-x]
    //                        | y <= q/2, so P(Y=y) is stored in distY[y]
    // We apply the change of variable x -> q - x
    if (boundOnY >= z)
    {
      if (boundOnX < boundOnY - z)
      {
        x = boundOnX;
        y = boundOnX + z;
      }
      else
      {
        x = boundOnY - z;
        y = boundOnY;
      }

      for (; 0u < x; x--, y--)
      {
        pz += gaDistribution1[x] * gaDistribution2[y];
      }
    }

    gaResult[z] = pz;
  }

  // For z greater than the bound on Z, the probability is 0
  for (; gModulus / 2u >= z; z += C_LEN_THREADS)
  {
    gaResult[z] = 0.0;
  }

  pthread_exit(NULL);
}

/**************************************************************************************************/
/** \implements multiplyDistributionsInternal
 * 
 **************************************************************************************************/
static void* multiplyDistributionsInternal
(
  void* xpThreadId
)
{
  size_t threadId = (size_t) xpThreadId;
  size_t x;
  size_t y;
  size_t z;
  size_t boundOnX = (size_t) gaDistribution1[C_OFS_BOUND];
  size_t boundOnY = (size_t) gaDistribution2[C_OFS_BOUND];
  TPDistribution aResult = gaResults[threadId];
  double px;

  for (z = 1u; z <= (gModulus / 2u); z++)
  {
    aResult[z] = 0.0;
  }

  // The x=0u case is treated by calling function
  for (x = threadId + 1u; x <= boundOnX; x += C_LEN_THREADS)
  {
    px = gaDistribution1[x];

    // The y=0u case is treated by calling function
    for (y = 1u; y <= boundOnY; y++)
    {
      // In a distribution array, only probabilities of x between 0 and q/2 are stored
      // When computing the product, we only store the values that lie in this range
      // However, a product can lie in this range without both terms being in it
      // So, we must consider all factors between 0 and q (or -q/2 and q/2)
      // For this reason, we must consider:
      //   z =   x  *   y  mod q
      //  -z =   x  * (-y) mod q
      //  -z = (-x) *   y  mod q
      //   z = (-x) * (-y) mod q
      // So, there are two ways to get the value among {-z,z} that lies in [0, q/2]
      z = (x * y) % gModulus;
      z = (z <= (gModulus / 2u)) ? z : (gModulus - z);
      aResult[z] += 2.0 * px * gaDistribution2[y];
    }
  }

  pthread_exit(NULL);
}

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */