/**************************************************************************************************/
/** \brief    Computations on symmetric distributions leveraging a CPU
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "distributions.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_LEN_ARRAYS  ((gModulus / 2u) + 2u)    //!< Length of the arrays storing distributions
#define C_LEN_THREADS 12u                       //!< Number of threads for parallel operations
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
static void* addDistributions_Internal
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
static void* multiplyDistributions_Internal
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
  size_t bound1 = (size_t) xpDistribution1[C_OFS_BOUND];
  size_t bound2 = (size_t) xpDistribution2[C_OFS_BOUND];
  size_t bound3 = ((bound1 + bound2) <= (gModulus / 2u)) ? (bound1 + bound2) : (gModulus / 2u);

  xpResult[C_OFS_BOUND] = (double) bound3;

  gaDistribution1 = (TPDistribution) xpDistribution1;
  gaDistribution2 = (TPDistribution) xpDistribution2;
  gaResult = xpResult;

  pthread_t threads[C_LEN_THREADS];
  size_t threadNum;

  for (threadNum = 0u; C_LEN_THREADS > threadNum; threadNum++)
  {
    int status = pthread_create(&threads[threadNum], NULL, addDistributions_Internal, (void*) threadNum);

    if (0 != status)
    {
      fprintf(stderr, "Error: pthread_create failed with status %d\n", status);
      exit(EXIT_FAILURE);
    }
  }

  for (threadNum = 0u; C_LEN_THREADS > threadNum; threadNum++)
  {
    pthread_join(threads[threadNum], NULL);
  }
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
  size_t bound1 = (size_t) xpDistribution1[C_OFS_BOUND];
  size_t bound2 = (size_t) xpDistribution2[C_OFS_BOUND];
  size_t bound3 = ((bound1 * bound2) <= (gModulus / 2u)) ? (bound1 * bound2) : (gModulus / 2u);

  xpResult[C_OFS_BOUND] = (double) bound3;

  gaDistribution1 = (TPDistribution) xpDistribution1;
  gaDistribution2 = (TPDistribution) xpDistribution2;

  size_t threadNum;
  pthread_t aThreads[C_LEN_THREADS];

  for (threadNum = 0u; C_LEN_THREADS > threadNum; threadNum++)
  {
    // Each thread will output a partial result to its own dedicated buffer
    gaResults[threadNum] = allocateDistribution();
    gaResults[threadNum][C_OFS_BOUND] = bound3;

    int status = pthread_create(&aThreads[threadNum], NULL, multiplyDistributions_Internal, (void*) threadNum);

    if (0 != status)
    {
      fprintf(stderr, "Error: pthread_create failed with status %d\n", status);
      exit(EXIT_FAILURE);
    }
  }

  for (threadNum = 0u; C_LEN_THREADS > threadNum; threadNum++)
  {
    pthread_join(aThreads[threadNum], NULL);
  }

  // P([X = 0] OR [Y = 0]) = P(X = 0) + P(Y = 0) - P([X = 0] AND [Y = 0])
  xpResult[0u] = xpDistribution1[0u] + xpDistribution2[0u] - xpDistribution1[0u] * xpDistribution2[0u];

  size_t x;

  for (x = 1u; x <= bound3; x++)
  {
    xpResult[x] = 0.0;
  }

  // Sum the partial results computed by the threads
  for (threadNum = 0u; C_LEN_THREADS > threadNum; threadNum++)
  {
    for (x = 1u; x <= bound3; x++)
    {
      xpResult[x] += gaResults[threadNum][x];
    }
  }
}

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - IMPLEMENTATION                                                               */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \implements addDistributions_Internal
 * 
 **************************************************************************************************/
static void* addDistributions_Internal
(
  void* xpThreadId
)
{
  size_t bound1 = (size_t) gaDistribution1[C_OFS_BOUND];
  size_t bound2 = (size_t) gaDistribution2[C_OFS_BOUND];
  size_t bound3 = (size_t) gaResult[C_OFS_BOUND];

  for (size_t z = (size_t) xpThreadId; z <= bound3; z += C_LEN_THREADS)
  {
    double probabilityOfZ = 0.0;

    // y = z - x
    size_t x = (z < bound2) ? 0u : (z - bound2);
    size_t y = (z < bound2) ? z : bound2;
    size_t end = (z < bound1) ? z : bound1;

    for (; x <= end; x++, y--)
    {
      probabilityOfZ += gaDistribution1[x] * gaDistribution2[y];
    }

    // y = x - z
    x = z + 1u;
    y = 1u;
    end = (bound1 < (bound2 + z)) ? bound1 : (bound2 + z);

    for (; x <= end; x++, y++)
    {
      probabilityOfZ += gaDistribution1[x] * gaDistribution2[y];
    }

    // y = x + z
    x = 1u;
    y = z + 1u;
    end = (bound2 > z) ? (bound2 - z) : 0u;
    end = (bound1 < end) ? bound1 : end;

    for (; x <= end; x++, y++)
    {
      probabilityOfZ += gaDistribution1[x] * gaDistribution2[y];
    }

    // y = q - x - z
    x = gModulus - z - bound2;
    y = bound2;
    end = bound1;

    for (; x <= end; x++, y--)
    {
      probabilityOfZ += gaDistribution1[x] * gaDistribution2[y];
    }

    gaResult[z] = probabilityOfZ;
  }

  pthread_exit(NULL);
}

/**************************************************************************************************/
/** \implements multiplyDistributions_Internal
 * 
 **************************************************************************************************/
static void* multiplyDistributions_Internal
(
  void* xpThreadId
)
{
  size_t threadId = (size_t) xpThreadId;
  TPDistribution aResult = gaResults[threadId];
  size_t bound1 = (size_t) gaDistribution1[C_OFS_BOUND];
  size_t bound2 = (size_t) gaDistribution2[C_OFS_BOUND];
  size_t bound3 = (size_t) aResult[C_OFS_BOUND];
  size_t z;

  for (z = 1u; z <= bound3; z++)
  {
    aResult[z] = 0.0;
  }

  // The [x = 0u] case is treated by the calling function
  for (size_t x = threadId + 1u; x <= bound1; x += C_LEN_THREADS)
  {
    double probabilityOfX = gaDistribution1[x];

    // The [y = 0u] case is treated by the calling function
    for (size_t y = 1u; y <= bound2; y++)
    {
      z = (x * y) % gModulus;
      z = (z <= (gModulus / 2u)) ? z : (gModulus - z);
      // If z = x * y, then we also have: z = (-x) * (-y)
      aResult[z] += 2.0 * probabilityOfX * gaDistribution2[y];
    }
  }

  pthread_exit(NULL);
}

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */