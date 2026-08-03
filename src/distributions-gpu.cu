/**************************************************************************************************/
/** \brief    Computations on symmetric distributions leveraging GPUs
 *
 *  \author   Julien CAM
 *
 *  \date     2026/07/21
 *
 *  \file
 **************************************************************************************************/

/* ---------------------------------------------------------------------------------------------- */
/* IMPORTS                                                                                        */
/* ---------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include "distributions.hpp"

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL CONSTANTS, TYPES, ENUM                                                                   */
/* ---------------------------------------------------------------------------------------------- */

#define C_LEN_ARRAYS  ((gModulus / 2u) + 2u)    //!< Length of the arrays storing distributions
#define C_OFS_BOUND   ((gModulus / 2u) + 1u)    //!< Offset to the bound on the distribution

#define C_NB_DEVICES     1u                     //!< Number of available GPUs
#define C_LEN_BLOCK   1024u                     //!< Number of threads per block

#define M_CHECK_STATUS(xStatus)           \
  if (cudaSuccess != xStatus)             \
  {                                       \
    printf("%s in '%s' at line '%d'.\n",  \
      cudaGetErrorString(xStatus),        \
      __FILE__,                           \
      __LINE__                            \
    );                                    \
    exit(EXIT_FAILURE);                   \
  }

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL VARIABLES                                                                                */
/* ---------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------- */
/* LOCAL FUNCTIONS - PROTOTYPE                                                                    */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \brief  Compute one element of the sum of two symmetric distributions
 *
 * \param[in]       xpDistribution1         Pointer to the first symmetric distribution Dx
 * \param[in]       xpDistribution2         Pointer to the second symmetric distribution Dy
 * \param[in]       xBound1                 The last integer with non-zero probability for Dx
 * \param[in]       xBound2                 The last integer with non-zero probability for Dy
 * \param[in]       xBound3                 The last integer with non-zero probability for Dx + Dy
 * \param[in]       xModulus                The modulus to apply to the sum
 * \param[in]       xStart                  The first element to handle on the current GPU
 * \param[out]      xpResult                Pointer to where the result Dz must be stored
 *
 **************************************************************************************************/
__global__ static void addDistributions_Internal
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  const size_t xBound1,
  const size_t xBound2,
  const size_t xBound3,
  const size_t xModulus,
  const size_t xStart,
  TPDistribution xpResult
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

  // The (bound3+1) threads are organized into ceiling((bound3+1) / C_LEN_BLOCK) blocks
  size_t nbBlocks = (bound3 + C_LEN_BLOCK) / C_LEN_BLOCK;

  TPDistribution apgInput1[C_NB_DEVICES];
  TPDistribution apgInput2[C_NB_DEVICES];
  TPDistribution apgOutput[C_NB_DEVICES];
  size_t aOutputLength[C_NB_DEVICES];

  size_t start = 0u;
  size_t deviceNum = 0u;

  while (start <= bound3)
  {
    M_CHECK_STATUS(cudaSetDevice(deviceNum));

    // Each GPU will execute (nbBlocks / C_NB_DEVICES) blocks,
    size_t nbBlocksOnThisDevice = nbBlocks / C_NB_DEVICES;
    // except the first (nbBlocks % C_NB_DEVICES) devices, which will execute one more
    nbBlocksOnThisDevice += (deviceNum < (nbBlocks % C_NB_DEVICES)) ? 1u : 0u;

    // The current GPU handles the output indexes in [start, end]
    size_t end = start + (nbBlocksOnThisDevice * C_LEN_BLOCK) - 1u;
    end = (end > bound3) ? bound3 : end;
    aOutputLength[deviceNum] = end - start + 1u;

    M_CHECK_STATUS(cudaMalloc(&apgInput1[deviceNum], C_LEN_ARRAYS * sizeof(double)));
    M_CHECK_STATUS(cudaMalloc(&apgInput2[deviceNum], C_LEN_ARRAYS * sizeof(double)));
    M_CHECK_STATUS(cudaMalloc(&apgOutput[deviceNum], aOutputLength[deviceNum] * sizeof(double)));

    M_CHECK_STATUS(cudaMemcpy(apgInput1[deviceNum], xpDistribution1, C_LEN_ARRAYS * sizeof(double), cudaMemcpyHostToDevice));
    M_CHECK_STATUS(cudaMemcpy(apgInput2[deviceNum], xpDistribution2, C_LEN_ARRAYS * sizeof(double), cudaMemcpyHostToDevice));

    addDistributions_Internal<<<nbBlocksOnThisDevice,C_LEN_BLOCK>>>(
      apgInput1[deviceNum], apgInput2[deviceNum],
      bound1, bound2, bound3,
      gModulus,
      start,
      apgOutput[deviceNum]
    );

    start = end + 1u;
    deviceNum++;
  }

  start = 0u;
  deviceNum = 0u;

  while (start <= bound3)
  {
    M_CHECK_STATUS(cudaSetDevice(deviceNum));

    M_CHECK_STATUS(cudaMemcpy(&xpResult[start], apgOutput[deviceNum], aOutputLength[deviceNum] * sizeof(double), cudaMemcpyDeviceToHost));

    M_CHECK_STATUS(cudaFree(apgInput1[deviceNum]));
    M_CHECK_STATUS(cudaFree(apgInput2[deviceNum]));
    M_CHECK_STATUS(cudaFree(apgOutput[deviceNum]));

    start += aOutputLength[deviceNum];
    deviceNum++;
  }

  xpResult[C_OFS_BOUND] = (double) bound3;
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
  size_t z;

  xpResult[C_OFS_BOUND] = (double) bound3;

  // P([X = 0] OR [Y = 0]) = P(X = 0) + P(Y = 0) - P([X = 0] AND [Y = 0])
  xpResult[0u] = xpDistribution1[0u] + xpDistribution2[0u] - xpDistribution1[0u] * xpDistribution2[0u];

  for (z = 1u; z <= bound3; z++)
  {
    xpResult[z] = 0.0;
  }

  // The [x = 0u] case is treated above
  for (size_t x = 1u; x <= bound1; x++)
  {
    double probabilityOfX = xpDistribution1[x];

    // The [y = 0u] case is treated above
    for (size_t y = 1u; y <= bound2; y++)
    {
      z = (x * y) % gModulus;
      z = (z <= (gModulus / 2u)) ? z : (gModulus - z);
      // If z = x * y, then we also have: z = (-x) * (-y)
      xpResult[z] += 2.0 * probabilityOfX * xpDistribution2[y];
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
__global__ static void addDistributions_Internal
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  const size_t xBound1,
  const size_t xBound2,
  const size_t xBound3,
  const size_t xModulus,
  const size_t xStart,
  TPDistribution xpResult
)
{
  size_t i = blockIdx.x * C_LEN_BLOCK + threadIdx.x;
  size_t z = xStart + i;

  if (xBound3 < z)
  {
    return;
  }

  double probability = 0.0;

  // y = z - x
  size_t x = (z < xBound2) ? 0u : (z - xBound2);
  size_t y = (z < xBound2) ? z : xBound2;
  size_t end = (z < xBound1) ? z : xBound1;

  for (; x <= end; x++, y--)
  {
    probability += xpDistribution1[x] * xpDistribution2[y];
  }

  // y = x - z
  x = z + 1u;
  y = 1u;
  end = (xBound1 < (xBound2 + z)) ? xBound1 : (xBound2 + z);

  for (; x <= end; x++, y++)
  {
    probability += xpDistribution1[x] * xpDistribution2[y];
  }

  // y = x + z
  x = 1u;
  y = z + 1u;
  end = (xBound2 > z) ? (xBound2 - z) : 0u;
  end = (xBound1 < end) ? xBound1 : end;

  for (; x <= end; x++, y++)
  {
    probability += xpDistribution1[x] * xpDistribution2[y];
  }

  // y = q - x - z
  x = xModulus - z - xBound2;
  y = xBound2;
  end = xBound1;

  for (; x <= end; x++, y--)
  {
    probability += xpDistribution1[x] * xpDistribution2[y];
  }

  xpResult[i] = probability;
}

/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */