/**************************************************************************************************/
/** \brief    Computations on symmetric distributions
 * 
 *  \author   Julien CAM
 * 
 *  \date     2025/03/26
 *
 *  \file
 **************************************************************************************************/
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

/* ---------------------------------------------------------------------------------------------- */
/* IMPORTS                                                                                        */
/* ---------------------------------------------------------------------------------------------- */

#include <stddef.h>

/* ---------------------------------------------------------------------------------------------- */
/* CONSTANTS, TYPES, ENUM                                                                         */
/* ---------------------------------------------------------------------------------------------- */

typedef double* TPDistribution;

/* ---------------------------------------------------------------------------------------------- */
/* VARIABLES                                                                                      */
/* ---------------------------------------------------------------------------------------------- */

// The parameter for modular arithmetic. Must be assigned a value in the main program.
extern size_t gModulus;

/* ---------------------------------------------------------------------------------------------- */
/* FUNCTIONS                                                                                      */
/* ---------------------------------------------------------------------------------------------- */

/**************************************************************************************************/
/** \brief  Add two symmetric distributions:
 *          If X is a random variable with symmetric distribution Dx
 *          and Y is a random variable with symmetric distribution Dy,
 *          compute the distribution Dz of Z = X + Y mod gModulus
 *
 * \param[in]       xpDistribution1         Pointer to the first symmetric distribution Dx
 * \param[in]       xpDistribution2         Pointer to the second symmetric distribution Dy
 * \param[out]      xpResult                Pointer to where the result Dz must be stored
 *
 **************************************************************************************************/
void addDistributions
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  TPDistribution xpResult
);

/**************************************************************************************************/
/** \brief  Compute the probability that a random variable X is rounded to 1 when reduced to 0 or 1.
 *          This corresponds to the probability that X is in the interval [q/4, 3q/4]
 *
 * \param[in]       xpDistribution          Pointer to the symmetric distribution of X
 *
 * \return          The probability that X is rounded to 1
 *
 **************************************************************************************************/
double computeRoundingToOneProbability
(
  const TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Copy a distribution from a buffer into another one
 *
 * \param[in]       xpInput                 The distribution to copy
 * \param[out]      xpOutput                The destination buffer
 *
 **************************************************************************************************/
void copyDistribution
(
  TPDistribution xpInput,
  TPDistribution xpOutput
);

/**************************************************************************************************/
/** \brief  Initialize a centered binomial distribution
 *
 * \param[in]       xEta                    The parameter of the distribution
 * \param[out]      xpDistribution          Pointer to the output symmetric distribution
 *
 **************************************************************************************************/
void initCenteredBinomialDistribution
(
  const size_t xEta,
  TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Initialize a centered normal distribution
 *
 * \param[in]       xStdDev                 The standard deviation of the normal distribution
 * \param[out]      xpDistribution          Pointer to the output symmetric distribution
 *
 **************************************************************************************************/
void initCenteredNormalDistribution
(
  const double xStdDev,
  TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Initialize a centered uniform distribution
 *
 * \param[in]       xBound                  The upper bound of the distribution
 * \param[out]      xpDistribution          Pointer to the output symmetric distribution
 *
 **************************************************************************************************/
void initCenteredUniformDistribution
(
  const size_t xBound,
  TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Initialize a distribution of compression/decompression error
 *
 * \param[in]       xNbCompressedBits       The number of bits used for compression
 * \param[out]      xpDistribution          Pointer to the output symmetric distribution
 *
 **************************************************************************************************/
void initCompressionErrorDistribution
(
  const size_t xNbCompressedBits,
  TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Allocate a buffer that can be used to store a distribution
 *
 * \return  Pointer to the output symmetric distribution
 *
 **************************************************************************************************/
TPDistribution allocateDistribution
(
  void
);

/**************************************************************************************************/
/** \brief  Load a symmetric distribution from a file
 *
 * \param[in]       filename                Name of the file where the distribution is stored
 * \param[out]      xpDistribution          Pointer to where the distribution must be stored
 *
 **************************************************************************************************/
void loadDistribution
(
  const char* filename,
  TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Multiply two symmetric distributions:
 *          If X is a random variable with symmetric distribution Dx
 *          and Y is a random variable with symmetric distribution Dy,
 *          compute the distribution Dz of Z = X * Y mod gModulus
 *
 * \param[in]       xpDistribution1         Pointer to the first symmetric distribution Dx
 * \param[in]       xpDistribution2         Pointer to the second symmetric distribution Dy
 * \param[out]      xpResult                Pointer to where the result Dz must be stored
 *
 **************************************************************************************************/
void multiplyDistributions
(
  const TPDistribution xpDistribution1,
  const TPDistribution xpDistribution2,
  TPDistribution xpResult
);

/**************************************************************************************************/
/** \brief  Print a symmetric distribution (only display elements with non-zero probability)
 *
 * \param[in]       xpDistribution          Pointer to the symmetric distribution to print
 *
 **************************************************************************************************/
void printDistribution
(
  const TPDistribution xpDistribution
);

/**************************************************************************************************/
/** \brief  Save a symmetric distribution to a file
 *
 * \param[in]       xpDistribution          Pointer to the symmetric distribution to save
 * \param[in]       filename                Name of the file where to save the distribution
 *
 **************************************************************************************************/
void saveDistribution
(
  const TPDistribution xpDistribution,
  const char* filename
);

#endif // DISTRIBUTIONS_H
/* ---------------------------------------------------------------------------------------------- */
/* END OF FILE                                                                                    */
/* ---------------------------------------------------------------------------------------------- */