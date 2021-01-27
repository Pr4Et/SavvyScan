/** @file preLR.cpp
 *
 * \brief Pre-LR filter adapted for IMOD.
 *
 * \mainpage Pre-LR filter in IMOD.
 * preNID is pre-reconstruction non-conservative non-linear isotropic diffusion (NID) filter that automatically identifies and reduces local irregularities in the tilt projections (see Maiorca M., Millet C., Hanssen E., Abbey B., Kazmierczak E., Tilley L., J Struct Biol. 2014 Apr;186(1):28-37). In particular, the filter targets highly electron dense regions, such as gold nanoparticles, that can hide proximal biological structures and degrade the overall quality of the reconstructed tomograms.
 *
 * \tableofcontents
 *
 * @param inputFile
 * @param "-input OR -InputStack"   Input Aligned tilt stack file
 * @param "-output OR -OutputFileName"    Output Aligned tilt stack file
 * @param "-angles OR -AnglesFile"	  tlt file (file with the tilt angles)
 * @param "-s OR -sigma"  Float \f$\sigma\f$ sigma
 * @param "-ite OR -Iterations"   number of iterations
 * @param "-views OR -ViewsToProcess"  list of view to process (optional)
 *
 *
 *  @author   Mauro Maiorca, of the Biochemistry & Molecular Biology Department, Bio21 institute, University of Melbourne, Australia, contributed the preNAD program. 
 *  It uses recursive line filter routines from Gregoire Malandain, covered by version 3 of the GPL (see GPL-3.0.txt). 
 *  This work was supported by funding from the Australian Research Council and the National Health and Medical Research Council.
 *
 */
#include <limits>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

#include "b3dutil.h"
#include "iimage.h"
#include "sliceproc.h"
#include "ctfutils.h"
#include "parse_params.h"
#include "nrutil.c"
#include "dsyev3.h"
#include "lapackc.h"
#include "recline.h"

#define MIN_ANGLE 1.0e-6	//tilt Angle less than this is treated as 0.0;
#define MY_PI 3.141592653589
#define MAXLINE 160

/**
 * PreLRTiltProjectionType: support structure for iterating through slices.
 *
 */
typedef struct
{

  /// angle in radiants 
  double angleRadiant;

  /// angle  in degrees
  double angleDegrees;

  /// \f$\sigma\f$ for the current tilt projection 
  double sigma;

  /// The current tilt projection won't be processed if skipMe is set to true
  bool skipMe;

  /// Current tilt projection number
  int TiltNumber;

  /// pointer to the processed tilt projection data at the current iteration
  float **Data;

  /// pointer to the original tilt projection data (prior to processing)
  float **DataOriginal;

} PreLRTiltProjectionType;


/// Accessory function: parsing a comma separated string
int
separateStringCommaValues (std::string InputS,	//!< input string
			   std::vector < float >&values	//!< output vector where to put the values
  )
/**
 * This function is used for parsing multiple sigma values, comma separated.
 */
{

  values.erase (values.begin (), values.end ());
  std::vector < std::string > strings;
  std::istringstream f (InputS);
  std::string s;
  while (std::getline (f, s, ','))
    {
      values.push_back (atof (s.c_str ()));
    }

  return 0;

}



/// Accessory function: recompute the frame of an image as a mirror would do.
void
FillingPadding (float **I,	//!<  INPUT/OUTPUT image, the dimension of this image is unchanged
		unsigned int nx,	//!<  size INPUT in x direction in pixels
		unsigned int ny,	//!<  size INPUT in y direction in pixels
		unsigned int paddingSize	//!<  padding size in pixels
  )
{

/**
 * 
 * This function it ignores previous pixel values in the outer frame of the image, and recompute them as a mirror would do.
 *
 */
  unsigned int p, i, j;		// loop variables
  if (nx > paddingSize + 1 && ny > paddingSize + 1)
    for (p = 0; p < paddingSize; p++)
      {
	for (i = 0; i < nx + 2 * paddingSize; i++)
	  I[i][paddingSize + ny + p] = I[i][ny - 1 - p];
	for (i = 0; i < nx + 2 * paddingSize; i++)
	  I[i][p] = I[i][2 * paddingSize - p];
	for (j = 0; j < ny + 2 * paddingSize; j++)
	  I[paddingSize + nx + p][j] = I[nx - 1 - p][j];
	for (j = 0; j < ny + 2 * paddingSize; j++)
	  I[p][j] = I[2 * paddingSize - p][j];
      }

  return;
}







/// Accessory function for computing 1D zero order and first order gaussian derivatives of a 2D image. It uses recursive line filter routines from Gregoire Malandain, covered by version 3 of the GPL (see GPL-3.0.txt). 
void gaussRecursiveDerivatives1D (double sigma,	//!< standard deviation of Gaussian
				  unsigned int nx,	//!< image dimension in x direction
				  unsigned int ny,	//!< image dimension in y direction
				  unsigned int padding,	//!< padding size (in pixels)
				  double hx,	//!< pixel size in x direction
				  double hy,	//!< pixel size in y direction
				  unsigned int direction,	//!< 0 for X, 1 for Y
				  unsigned int DerivativeOrder,	//!< 0 smoothing, 1 first derivative, 2 second derivative
				  float **f,	//!< input: original image ;
				  float **fo = NULL	//!<output: smoothed (if null, input overwritten)
  )				// input: original image ;  output: smoothed
{

  unsigned int i, j, k;		// loop variables

  FillingPadding (f, nx, ny, padding);

  if (fo == NULL)
    {
      fo = f;
    }


  //if less then a certain threshold, then just don't smooth, copy the image
  const double sigmaThreshold = 0.1;
  if (sigma < sigmaThreshold)
    {
      for (i = padding; i < nx + padding; i++)
	for (j = padding; j < ny + padding; j++)
	  fo[i][j] = f[i][j];
      FillingPadding (fo, nx, ny, padding);
      return;
    }


  //for each line, put the line in a buffer and process it.
  double *bufferIn = NULL;
  double *bufferOut = NULL;
  double *bufferTmp0 = NULL;
  double *bufferTmp1 = NULL;
  //float filterCoefs;

  //recursiveFilterType recFilter = ALPHA_DERICHE;
  recursiveFilterType recFilter = GAUSSIAN_DERICHE;
  derivativeOrder derivOrder = NODERIVATIVE;
  if (DerivativeOrder == 0)
    {
      derivOrder = SMOOTHING;
    }
  else if (DerivativeOrder == 1)
    {
      derivOrder = DERIVATIVE_1;
    }
  else if (DerivativeOrder == 2)
    {
      derivOrder = DERIVATIVE_2;
    }
  else if (DerivativeOrder == 3)
    {
      derivOrder = DERIVATIVE_3;
    }
  RFcoefficientType *rfc =
    InitRecursiveCoefficients (sigma, recFilter, derivOrder);
  if (!rfc)
    exitError ("Allocation structure for recursive filter");

  //fill the buffer
  if (direction == 0)		//along X
    {


      bufferIn = (double *) malloc ((nx + 2 * padding) * sizeof (double));
      bufferOut = (double *) malloc ((nx + 2 * padding) * sizeof (double));
      bufferTmp0 = (double *) malloc ((nx + 2 * padding) * sizeof (double));
      bufferTmp1 = (double *) malloc ((nx + 2 * padding) * sizeof (double));

      if (bufferIn == NULL || bufferOut == NULL || bufferTmp0 == NULL
	  || bufferTmp1 == NULL)
	exitError ("error allocating memory");

      for (j = 0; j < ny + 2 * padding; j++)
	{
	  for (k = 0; k < nx + 2 * padding; k++)
	    {
	      bufferIn[k] = (double) f[k][j];
	      bufferOut[k] = (double) f[k][j];
	      bufferTmp0[k] = (double) f[k][j];
	      bufferTmp1[k] = (double) f[k][j];
	    }
	  RecursiveFilter1D (rfc, bufferIn, bufferOut, bufferTmp0, bufferTmp1,
			     nx + 2 * padding);

	  for (k = 0; k < nx + 2 * padding; k++)
	    {
	      fo[k][j] = (float) bufferOut[k];
	    }

	}

      free (bufferIn);
      free (bufferOut);
      free (bufferTmp0);
      free (bufferTmp1);



    }				//along Y
  if (direction == 1)
    {


      bufferIn = (double *) malloc ((ny + 2 * padding) * sizeof (double));
      bufferOut = (double *) malloc ((ny + 2 * padding) * sizeof (double));
      bufferTmp0 = (double *) malloc ((ny + 2 * padding) * sizeof (double));
      bufferTmp1 = (double *) malloc ((ny + 2 * padding) * sizeof (double));
      if (bufferIn == NULL || bufferOut == NULL || bufferTmp0 == NULL
	  || bufferTmp1 == NULL)
	exitError ("error allocating memory");


      for (j = 0; j < nx + 2 * padding; j++)
	{
	  for (k = 0; k < ny + 2 * padding; k++)
	    {
	      bufferIn[k] = (double) f[j][k];
	      bufferOut[k] = (double) f[j][k];
	      bufferTmp0[k] = (double) f[j][k];
	      bufferTmp1[k] = (double) f[j][k];
	    }
	  RecursiveFilter1D (rfc, bufferIn, bufferOut, bufferTmp0, bufferTmp1,
			     ny + 2 * padding);

	  for (k = 0; k < ny + 2 * padding; k++)
	    {
	      fo[j][k] = (float) bufferOut[k];
	    }

	}
      free (bufferIn);
      free (bufferOut);
      free (bufferTmp0);
      free (bufferTmp1);




    }

  FillingPadding (fo, nx, ny, padding);
  free (rfc);

  return;

}				//gaussRecursiveDerivatives1D






///Masks high irregular region on the image, replace irregular regions with a smoothed version of it, grade the irregularity of each pixel in a [0,1] interval
double
CreateMaskedLocalSmooth (unsigned int nx,	/*!< Dimension (in pixels) along the axis \f$x\f$ */
			 unsigned int ny,	/*!< Dimension (in pixels) along the axis \f$y\f$ */
			 unsigned int padding,	//!< padding size (in pixels)
			 double spacingX,	/*!< Size (in nm) of a pixel along the axis \f$x\f$ */
			 double spacingY,	/*!< Size (in nm) of a pixel along the axis \f$y\f$ */
			 double sigmaMask,	/*!< \f$\sigma\f$ of the Mask */
			 double alpha_,	/*!< \f$\alpha\f$ value */
			 double beta_,	/*!< \f$\beta\f$ value */
			 double tau,	/*!< \f$\tau_1\f$ value */
			 double angleRadiant,	/*!< angle of the ongoing tilt projection (in radiants) */
			 float **I,	/*!< Input image */
			 float **IMask,	/*!< Output mask image */
			 float **IGradeIrragularity,	/*!< Image with irregularities graded: 0 none, 1 max */
			 float **MaskedLocalSmoothedImage,	/*!< Output Masked Local Smoothed Image */
			 float **NormalizedDoG	/*!< Output Masked Local Smoothed Image */
  )
{
//      double sigmaMin = 0.1;

  //output mask file
  float **mask = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);

  //gaussian blurred images
  float **IBlurred = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **IBlurredA = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);

  //Difference of Gaussians
  float **DoG = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **relevantDoG = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **dDoGdx = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **dDoGdy = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);


  //gaussian blurred image
  gaussRecursiveDerivatives1D (sigmaMask, nx, ny, padding, spacingX, spacingY,
			       0, 0, I, IBlurred);
  gaussRecursiveDerivatives1D (sigmaMask, nx, ny, padding, spacingX, spacingY,
			       1, 0, IBlurred, IBlurred);

  double meanDoG = I[padding][padding] - IBlurred[padding][padding];
  double maxDoG = meanDoG;
  double minDoG = meanDoG;
  const double counter = (nx) * (ny);

  double nrot = 0.0;


		/**
		 *
		 * In electron microscopy, the most important source of intensity irregularities is provided by local spatial differences in electron
		 * absorption. In scale space theory, those local spatial irregularities can be retrieved by measuring the space variancy at a defined scale of an image.
		 * This can be achieved by subtracting from the original image a blurred image with a gaussian of a variance \f$\sigma\f$. 
		 *
		 *  \f{eqnarray*}{
		 *  \Gamma_{(\theta,\sigma)} = I - I_{\sigma}
		 *  \f}
		 */
  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {
	DoG[i][j] = I[i][j] - IBlurred[i][j];

      }

		/**
		 *
		 * In order to reduce the contribuition of isolated small regions (compared to \f$\sigma\f$), we blur the resulting image with a gaussian of size \f$ \sigma_\alpha = \alpha \sigma\f$,
		 * with a user defined parameter \f$ \alpha \f$, \f$ \alpha > 0\f$,
		 *  \f{eqnarray*}{
		 *  ({\Gamma_{(\theta,\sigma)}})_{\sigma_{\alpha}} = (I - I_{\sigma})_{\sigma_{\alpha}} 
		 *  \f}
		 *
		 */
  double sigmaMean = sigmaMask * alpha_;
  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       0, 0, DoG, DoG);
  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       1, 0, DoG, DoG);

  maxDoG = 0.0;
  minDoG = 0.0;
  for (unsigned int i = padding; i < nx + padding; i++)
    for (unsigned int j = padding; j < ny + padding; j++)
      {

	NormalizedDoG[i][j] = -DoG[i][j];
	if (NormalizedDoG[i][j] < 0)
	  NormalizedDoG[i][j] = 0;
	if (NormalizedDoG[i][j] > maxDoG)
	  maxDoG = NormalizedDoG[i][j];
      }
  //meanDoG = meanDoG/(counter + 1.0);
  printf ("     minDoG=%f; maxDoG=%f\n", minDoG, maxDoG);

  for (unsigned int i = padding; i < nx + padding; i++)
    for (unsigned int j = padding; j < ny + padding; j++)
      {

	NormalizedDoG[i][j] = (NormalizedDoG[i][j]) / maxDoG;

	if (NormalizedDoG[i][j] > 0.0)
	  {
	    mask[i][j] = 1.0;
	  }
	else
	  {
	    mask[i][j] = 0.0;
	  }

      }



	/**
	 *  In order to obtain a maximum where those values are zero and a minimum in correspondence to local invariants, we compute the squared magnitude of divergence of \f$(I - I_{\sigma})_{\sigma_\alpha} \f$, and again we blur it by convolving with a Gaussian of sigma \f$\sigma_\alpha\f$.
	 *
	 *  \f{eqnarray*}{
	 *    \Psi=|div ((I - I_{\sigma})_{\sigma_\alpha}) |^2_{\sigma_\alpha} =  ((\partial_{x}((I - I_{\sigma})_{\sigma_\alpha}))^2 + (\partial_{y}((I - I_{\sigma})_{\sigma_\alpha}))^2)_{\sigma_\alpha}
	 *  \f}
 	 *
	 */


  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       0, 1, DoG, dDoGdx);
  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       1, 1, DoG, dDoGdy);


  meanDoG =
    (pow ((double)dDoGdx[padding][padding], 2.0) +
     pow ((double)dDoGdy[padding][padding], 2.0));;
  maxDoG = meanDoG;
  minDoG = meanDoG;
  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {
	DoG[i][j] = (pow ((double)dDoGdx[i][j], 2.0) + pow ((double)dDoGdy[i][j], 2.0));
	if (DoG[i][j] > maxDoG)
	  maxDoG = DoG[i][j];
	if (DoG[i][j] < minDoG)
	  minDoG = DoG[i][j];
      }
  meanDoG = meanDoG / (counter + 1.0);

  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       0, 0, DoG, DoG);
  gaussRecursiveDerivatives1D (sigmaMean, nx, ny, padding, spacingX, spacingY,
			       1, 0, DoG, DoG);

  meanDoG = DoG[padding][padding];
  maxDoG = meanDoG;
  minDoG = meanDoG;
  double counterMask = 0;
  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {
	if (mask[i][j] > 0.1)
	  {
	    if (DoG[i][j] > maxDoG)
	      maxDoG = DoG[i][j];
	    if (DoG[i][j] < minDoG)
	      minDoG = DoG[i][j];
	    meanDoG += DoG[i][j];
	    counterMask++;
	  }
      }
  meanDoG = meanDoG / (counterMask + 1.0);


  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {

	//  \hat{\Gamma}=2\frac{(\Gamma (\theta,\sigma))_{\sigma_t}-\overline{(\Gamma (\theta,\sigma))_{\sigma_t}}}{max(\Gamma (\theta,\sigma)_{\sigma_t})-min(\Gamma (\theta,\sigma)_{\sigma_t})}
		/**
		 *  We are interested in regions with higher electron absorption compared to the their local neighborhood, thus all the values \f$ I - I_{\sigma} < 0\f$. 
		 *  For all the values in that range, we are interested in regions with the higher electron absorption irregularities. Thus, we are interested in regions
		 *  \f{eqnarray*}{
			\frac{\Psi-\overline{\Psi}}{max(\Psi)-min(\Psi)} > \tau
		 *  \f}
		 * where \f$\tau\f$ is a value in the interval [0,1], higher \f$\tau\f$, higher \f$\tau\f$ gives high specificity and low sensitivity, we fixed a neutral value of 0.01;
		 *  \f$\overline{\Psi}\f$, \f$max(\Psi)\f$ and \f$min(\Psi)\f$ respectively the mean, the maximum and the minimum vale of \f$\Psi\f$, 
		 * restricted to the regions with higher electron absorption described before (\f$ I - I_{\sigma} < 0\f$).
		 * 
		 */
	DoG[i][j] = (DoG[i][j] - meanDoG) / (maxDoG - minDoG);
	if (DoG[i][j] < 0.0)
	  DoG[i][j] = 0.0;

      }


	 /**
	 *  Finally, we mark with 0 all the values \f$ I - I_{\sigma} < 0\f$. We blur with Gaussian of sigma \f$\sigma_t\f$, and we linearly resample the mask image in the interval [0,1]
	 *  \f{eqnarray*}{
	 *    (Image-Min)/(max-min)
	 *  \f}
	 *  
	 * Where max and min are the maximum and minumum of the function. We refer to this image as a mask.
 	 */

	/**
	 * A final image is obtained by replacing the values of the mask with the blurred image.
	 * 
	 */

  //gaussian blurred image with adapted sigma
  gaussRecursiveDerivatives1D (beta_, nx, ny, padding, spacingX, spacingY, 0,
			       0, I, IBlurredA);
  gaussRecursiveDerivatives1D (beta_, nx, ny, padding, spacingX, spacingY, 1,
			       0, IBlurredA, IBlurredA);

  double IrragularityMax = 0.0;
  double IrragularityMin = 1.0;
  double IrragularityThreshold = 0.005;
  double newMaskValue;
  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {

	//if ( mask[i][j] > 0.1 ) IMask[i][j] = 0;
	//NO if ( mask[i][j] > 0.1 && DoG[i][j] > tau) IMask[i][j] = 0;
	if (mask[i][j] > 0.01 && DoG[i][j] > tau)
	  {
	    IMask[i][j] = DoG[i][j];
	    mask[i][j] = 1.0;
	    MaskedLocalSmoothedImage[i][j] = IBlurredA[i][j];
	  }
	else
	  {
	    MaskedLocalSmoothedImage[i][j] = I[i][j];
	    IMask[i][j] = 0.0;
	    mask[i][j] = 0.0;
	  }
	if (IMask[i][j] < IrragularityThreshold)
	  {
	    IMask[i][j] = 0;
	  }
	if (IMask[i][j] > IrragularityMax)
	  IrragularityMax = IMask[i][j];
	if (IMask[i][j] < IrragularityMin)
	  IrragularityMin = IMask[i][j];

	IMask[i][j] = mask[i][j];
      }



  const double rangeIrregularitues = IrragularityMax - IrragularityMin;
  for (unsigned int i = 0; i < nx + 2 * padding; i++)
    for (unsigned int j = 0; j < ny + 2 * padding; j++)
      {
	IGradeIrragularity[i][j] =
	  (IMask[i][j] - IrragularityMin) / rangeIrregularitues;
      }


  //clear buffers
  free_matrix (IBlurredA, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (IBlurred, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (mask, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (DoG, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (relevantDoG, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (dDoGdx, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (dDoGdy, 0, nx + 2 * padding, 0, ny + 2 * padding);



  return 0;
}


//*****************************
//
//  end of create MASK
//







///Linearity Enhancing Diffusion
double
LinearityEnhancingDiffusion (unsigned int nx,	/*!< Dimension (in pixels) along the axis \f$x\f$ */
			     unsigned int ny,	/*!< Dimension (in pixels) along the axis \f$y\f$ */
			     unsigned int padding,	//padding size
			     double spacingX,	/*!< Size (in nm) of a pixel along the axis \f$x\f$ */
			     double spacingY,	/*!< Size (in nm) of a pixel along the axis \f$y\f$ */
			     double sigma,	/*!< \f$\sigma\f$ */
			     double angleRadiant,	/*!< angle of the ongoing tilt projection (in radiants) */
			     int iterations,	/*!< number of iterations */
			     float **I,	/*!< Input image */
			     float **Iout,	/*!< Output image */
			     float **mask,	/*!< Mask with the graded edges */
			     float **maskInvariants	/*!< Mask with the Invariants */
  )
/**
 *  is the core function of the algorithm
 */
{

  if (Iout == NULL)
    {
      Iout = I;
    }






  double ContrastParameterLambdaLED = 30;
  double ThresholdParameterC = 3.31488;
  double zeroValueTolerance = 1e-15;
  double Alpha = 0.001;
  double timeStep = 0.1;
  double rxx = timeStep / (1.0 * spacingX * spacingX);
  double ryy = timeStep / (1.0 * spacingY * spacingY);
  double rxy = timeStep / (1.4142 * spacingX * spacingY);

  //structure tensor enhancers
  const unsigned int convolutionPrecision = sigma + 1.0;
  const double magnifyGradientRatio =
    (1.0 / convolutionPrecision) * log ((double) convolutionPrecision);



  for (unsigned int i = padding; i < nx + padding; i++)
    for (unsigned int j = padding; j < ny + padding; j++)
      {


	//EVOLUTION
	double wE = rxx * mask[i][j];
	double wW = rxx * mask[i][j];
	double wS = ryy * mask[i][j];
	double wN = ryy * mask[i][j];
	double wSE = rxy * mask[i][j];
	double wNW = rxy * mask[i][j];
	double wNE = rxy * mask[i][j];
	double wSW = rxy * mask[i][j];

	Iout[i][j] = I[i][j]
	  + wE * (I[i + 1][j] - I[i][j])
	  + wW * (I[i - 1][j] - I[i][j])
	  + wS * (I[i][j + 1] - I[i][j])
	  + wN * (I[i][j - 1] - I[i][j])
	  + wSE * (I[i + 1][j + 1] - I[i][j])
	  + wNW * (I[i - 1][j - 1] - I[i][j])
	  + wSW * (I[i - 1][j + 1] - I[i][j])
	  + wNE * (I[i + 1][j - 1] - I[i][j]);

      }



  return 0;


}





/**
 *  Automatic Stop condition and iteration through the whole dataset
 *
 */
void
automaticPreLR (
/// stack in
		 PreLRTiltProjectionType * stackIn,
/// stack out
		 PreLRTiltProjectionType * stackOut,
/// stack out
		 PreLRTiltProjectionType * maskIn,
/// dimension along the axis \f$x\f$
		 unsigned int nx,
/// dimension along the axis \f$y\f$
		 unsigned int ny,
/// number of tilt projections
		 unsigned int nz,
/// the size of border padding
		 unsigned int padding,
///list of \f$\sigma\f$ values
		 std::vector < float >sigmaList,
///list of \f$\alpha\f$ values
		 std::vector < float >alphaSigma_List,
///list of \f$\beta\f$ values
		 std::vector < float >betaSigma_List,
///list of \f$\tau_1\f$ values
		 std::vector < float >tau_List,
/// minimum number of iterations
		 std::vector < float >Iterations_List,
/// spacing (in nm) along  \f$x\f$
		 double spacingX,
/// spacing (in nm) along  \f$y\f$
		 double spacingY,
/// \f$\lambda_c\f$
		 double lambdaC,
/// prints the mask instead of the value
		 bool maskOutput = false,
/// prints the mask instead of the value
		 bool abemusMask = false)
{

  unsigned int ii = 0;

  printf ("START (t=tilt, a=angle, i=iteration)\n");


  float **J0 = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **J = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **J1 = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **J_smooth = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **mask = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **ResidualImage = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **maskDilated = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **IMaskDoG = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);
  float **IGradedIrragularities = matrix (0, nx + 2 * padding, 0, ny + 2 * padding);	// Image with irregularities graded: 0 none, 1 max



  //We have to reach the min number of iterations for each tilt angles
  for (ii = 0; ii < nz; ii++)
    {

      for (unsigned int i = 0; i < nx + 2 * padding; i++)
	for (unsigned int j = 0; j < ny + 2 * padding; j++)
	  {
	    J[i][j] = (stackIn[ii].Data)[i][j];
	    J0[i][j] = (stackIn[ii].Data)[i][j];
	    J1[i][j] = (stackIn[ii].Data)[i][j];
	  }

      std::vector < float >::iterator it = sigmaList.begin ();
      std::vector < float >::iterator alphaIt = alphaSigma_List.begin ();
      std::vector < float >::iterator betaIt = betaSigma_List.begin ();
      std::vector < float >::iterator tauIt = tau_List.begin ();
      std::vector < float >::iterator IterationsIt = Iterations_List.begin ();



      for (; it != sigmaList.end ();
	   ++it, ++alphaIt, ++betaIt, ++tauIt, ++IterationsIt)
	{

	  double sigma0 = *it;
	  double sigmaAdapted = sigma0;
	  double absAngle = abs (stackOut[ii].angleRadiant);


	  printf ("[ITERATE sigma=%f, alpha=%f, beta=%f, tau=%f,] ", sigma0,
		  *alphaIt, *betaIt, *tauIt);


	  for (unsigned int i = 0; i < nx + 2 * padding; i++)
	    for (unsigned int j = 0; j < ny + 2 * padding; j++)
	      {
		J[i][j] = J1[i][j];
	      }




	  if (!abemusMask)
	    {
	      //create the mask
	      //double Specificity = 0.99;
	      //double SpecificitySensitivityBalance = 0.5;
	      CreateMaskedLocalSmooth (nx, ny, padding, spacingX, spacingY,
				       sigmaAdapted, *alphaIt, *betaIt,
				       *tauIt, absAngle, J, mask,
				       IGradedIrragularities, J1, IMaskDoG);
	    }
	  else
	    {
	      //import the mask
	      for (unsigned int i = 0; i < nx + 2 * padding; i++)
		for (unsigned int j = 0; j < ny + 2 * padding; j++)
		  {
		    mask[i][j] = (maskIn[ii].Data)[i][j];
		  }
	    }

//              gaussRecursiveDerivatives1D (1.0, nx, ny, padding, spacingX, spacingY, 0, 0, IGradedIrragularities, IGradedIrragularities);
//              gaussRecursiveDerivatives1D (1.0, nx, ny, padding, spacingX, spacingY, 1, 0, IGradedIrragularities, IGradedIrragularities);
//              gaussRecursiveDerivatives1D (1.0, nx, ny, padding, spacingX, spacingY, 0, 0, mask, mask);
//              gaussRecursiveDerivatives1D (1.0, nx, ny, padding, spacingX, spacingY, 1, 0, mask, mask);


	  //float ** tmpI = J;
	  for (int iii = 0; iii < *IterationsIt; iii++)
	    {
	      //tmpI = J1;
	      //J1 = J;
	      //J=tmpI;
	      LinearityEnhancingDiffusion (nx, ny, padding, spacingX,
					   spacingY, *it,
					   stackOut[ii].angleRadiant, 1, J1,
					   J, mask, mask);
	      for (unsigned int i = 0; i < nx + 2 * padding; i++)
		for (unsigned int j = 0; j < ny + 2 * padding; j++)
		  {
		    J1[i][j] = J[i][j];
		  }

	    }


	}



      for (unsigned int i = padding; i < nx + padding; i++)
	for (unsigned int j = padding; j < ny + padding; j++)
	  {
	    //(stackOut[ii].Data)[i][j] = 100*IGradedIrragularities[i][j]*IMaskDoG[i][j] ;
	    //(stackOut[ii].Data)[i][j] = IMaskDoG[i][j] ;
	    //(stackOut[ii].Data)[i][j] = 100*mask[i][j] ;

	    if (maskOutput)
	      (stackOut[ii].Data)[i][j] = 100 * mask[i][j];
	    else
	      (stackOut[ii].Data)[i][j] = J1[i][j];

	  }

    }



  free_matrix (IGradedIrragularities, 0, nx + 2 * padding, 0,
	       ny + 2 * padding);
  free_matrix (J_smooth, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (ResidualImage, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (mask, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (J0, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (J, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (J1, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (maskDilated, 0, nx + 2 * padding, 0, ny + 2 * padding);
  free_matrix (IMaskDoG, 0, nx + 2 * padding, 0, ny + 2 * padding);



}


/**
 *  Main()
 *
 */
int
main (int argc, char *argv[])
{

  char *stackFn, *maskFn, *angleFn, *outFn, *ViewsToProcess, *sigma, *alpha,
    *beta, *tau, *Iterations;
  int maskOutput;
  int startingView;

  char *progname = imodProgName (argv[0]);

  int numOptArgs, numNonOptArgs;


  unsigned int padding = 40;
  bool rangeRequested = false;
  bool noTiltAngle = false;
  bool gotMaskOutput = false;
  bool gotMaskInput = false;

  std::vector < float >sigmaList;
  std::vector < float >alphaSigma_List;
  std::vector < float >betaSigma_List;
  std::vector < float >tau_List;
  std::vector < float >Iterations_List;

  //parselist of views
  int *anglesToProcess;
  int nAnglesToProcess = 0;

  // Fallbacks from    ../manpages/autodoc2man 2 1 preNID
  int numOptions = 11;
  const char *options[] = {
    "input:InputStack:FN:", "output:OutputFileName:FN:", "angles:AnglesFile:FN:",
    "s:Sigma:FA:", "a:Alpha:FA:", "b:Beta:FA:", "t:Tau:FA:", "ite:Iterations:IA:",
    "im:InputMask:FN:", "mask:MaskOutput:B:", "views:ViewsToProcess:LI:"};

  printf
    ("\n [Mauro Maiorca, of the Biochemistry & Molecular Biology Department, Bio21 institute, University of Melbourne, Australia, contributed the preNID program (adapted for IMOD). ");
  printf
    ("It uses recursive line filter routines from Gregoire Malandain, covered by version 3 of the GPL (see GPL-3.0.txt). ");
  printf (" examples of use:\n\n");

  PipReadOrParseOptions (argc, argv, options, numOptions, progname, 1, 0, 0,
			 &numOptArgs, &numNonOptArgs, NULL);
  if (!PipGetBoolean("usage", &startingView)) {
    PipPrintHelp(progname, 0, 0, 0);
    exit(0);
  }
  if (PipGetString ("InputStack", &stackFn))
    exitError ("No stack specified");
  if (PipGetString ("InputMask", &maskFn))
    {
      maskFn = (char *) "none";
      gotMaskInput = false;
    }
  else
    {
      gotMaskInput = true;
    }
  if (PipGetString ("OutputFileName", &outFn))
    exitError ("OutputFileName is not specified");
  if (!PipGetString ("ViewsToProcess", &ViewsToProcess))
    {
      rangeRequested = true;
      anglesToProcess = parselist (ViewsToProcess, &nAnglesToProcess);
    }
  if (PipGetString ("AnglesFile", &angleFn))
    {
      angleFn = NULL;
      printf
	("[Warning] No angle file is specified, tilt angles are assumed to be all 0.0 degrees.\n");
      noTiltAngle = true;
    }
  if (!PipGetString ("Sigma", &sigma))
    {
      separateStringCommaValues (sigma, sigmaList);
      if (sigmaList.size () < 1)
	exitError ("No sigma is specified, aborting...");
    }
  else
    {
      exitError ("No sigma is specified, aborting...");
    }
  if (!PipGetString ("Alpha", &alpha))
    {
      separateStringCommaValues (alpha, alphaSigma_List);
      if (alphaSigma_List.size () != sigmaList.size ())
	{
	  exitError
	    ("Error with '--Alpha' option: be sure it is consistent with the '--Sigma' option. Aborting...");
	}
    }
  else
    {
      alphaSigma_List.resize (sigmaList.size ());
      for (std::vector < float >::iterator it = alphaSigma_List.begin ();
	   it != alphaSigma_List.end (); ++it)
	*it = 0.5;
      printf ("No alpha inserted, using default value 0.5\n");
    }
  if (!PipGetBoolean ("MaskOutput", &maskOutput))
    {
      gotMaskOutput = true;
    }
  if (!PipGetString ("Beta", &beta))
    {
      separateStringCommaValues (beta, betaSigma_List);
      if (betaSigma_List.size () != sigmaList.size ())
	{
	  exitError
	    ("Error with '--Beta' option: be sure it is consistent with the '--Sigma' option. Aborting...");
	}
    }
  else
    {
      betaSigma_List.resize (sigmaList.size ());
      for (std::vector < float >::iterator it = betaSigma_List.begin ();
	   it != betaSigma_List.end (); ++it)
	*it = 0.5;
      printf ("No beta inserted, using default value 0.5\n");
    }
  if (!PipGetString ("Tau", &tau))
    {
      separateStringCommaValues (tau, tau_List);
      if (tau_List.size () != sigmaList.size ())
	{
	  exitError
	    ("Error with '--Tau' option: be sure it is consistent with the '--Sigma' option. Aborting...");
	}

    }
  else
    {
      tau_List.resize (sigmaList.size ());
      for (std::vector < float >::iterator it = tau_List.begin ();
	   it != tau_List.end (); ++it)
	*it = 0.1;
      printf ("No tau inserted, using default value 0\n");
    }


  if (!PipGetString ("Iterations", &Iterations))
    {
      separateStringCommaValues (Iterations, Iterations_List);
      if (Iterations_List.size () != sigmaList.size ())
	{
	  exitError
	    ("Error with '--Iterations' option: be sure it is consistent with the '--sigma' option. Aborting...");
	}

    }
  else
    {
      Iterations_List.resize (sigmaList.size ());
      for (std::vector < float >::iterator it = Iterations_List.begin ();
	   it != Iterations_List.end (); ++it)
	*it = 1.0;
      printf ("No Iterations inserted, using default value 1.0\n");
    }


  printf ("parameters:\n");
  printf
    ("       stackFn = \"%s\", outFn=\"%s\",  angleFn=\"%s InputMask=\"%s\"",
     stackFn, outFn, angleFn, maskFn);
  printf ("\n       ");

  FILE *fpStack;
  if ((fpStack = iiFOpen (stackFn, "rb")) == 0)
    {
      exitError ("could not open input file %s", stackFn);
      printf ("error 0.5\n");
    }

  FILE *foutput;

  MrcHeader header;
  MrcHeader outHeader;
  int sliceMode;
  /* read header */
  if (mrc_head_read (fpStack, &header))
    {
      printf ("Error 2\n");
      exitError ("reading header of input file %s", stackFn);

    }

  if (mrc_head_read (fpStack, &outHeader))
    {
      printf ("Error 3\n");
      exitError ("reading header of input file %s", stackFn);
    }

  //check the file
  sliceMode = sliceModeIfReal (header.mode);
  if (sliceMode < 0)
    exitError ("File mode is %d; only byte, short integer, or real allowed",
	       header.mode);
  printf ("slice mode=%d\n", sliceMode);

  startingView = 1;
  int numSlices = header.nz;
  printf ("Number of slices=%d\n", numSlices);

  //write the header
  unsigned int nx = header.nx;
  unsigned int ny = header.ny;
  Islice *currSlice;

  float angleSign = 1;
  float minAngle, maxAngle;
  minAngle = -20;
  maxAngle = 20;

  float *tiltAngles = NULL;
  printf ("reading angle file=%s\n", angleFn);
  if (!noTiltAngle)
    tiltAngles =
      readTiltAngles (angleFn, header.nz, angleSign, minAngle, maxAngle);


  // printf("minAngle=%f, maxAngle=%f\n",minAngle,maxAngle);
  foutput = iiFOpen (outFn, "wb");

  // DNM: Set header correcty for a new output file and copy the extended header
  mrcInitOutputHeader(&outHeader);
  outHeader.fp = foutput;
  ImodImageFile *iiFile = iiLookupFileFromFP(fpStack);
  if (iiFile && iiFile->file == IIFILE_MRC && b3dOutputFileType() == OUTPUT_TYPE_MRC) {
    int err = mrcCopyExtraHeader(&header, &outHeader);
    if (err)
      exitError("Copying MRC extended header to output file (error # %d)", err);
  }

  mrc_head_label (&outHeader, "PreNID filtered image");
  //mrc_head_label (&outHeader, "PreNID filtered image, options here");
  mrc_head_write (foutput, &outHeader);



  PreLRTiltProjectionType *stackIn = NULL;
  PreLRTiltProjectionType *stackOut = NULL;
  PreLRTiltProjectionType *stackMask = NULL;
  stackIn =
    (PreLRTiltProjectionType *) malloc (header.nz *
					sizeof (PreLRTiltProjectionType));
  stackOut =
    (PreLRTiltProjectionType *) malloc (header.nz *
					sizeof (PreLRTiltProjectionType));
  stackMask =
    (PreLRTiltProjectionType *) malloc (header.nz *
					sizeof (PreLRTiltProjectionType));

  if (stackIn == NULL || stackOut == NULL || stackMask == NULL)
    exitError ("error allocating memory");


  double spacingX = 1.0;
  double spacingY = 1.0;
  double lambdaE = 30.0, lambdaC = 30.0, lambdaH = 30.0;
  unsigned int ii = 0;

  //set the angles to process
  for (ii = 0; ii < (unsigned int) header.nz && rangeRequested; ii++)
    stackIn[ii].skipMe = true;

  if (rangeRequested)
    {
      for (ii = 0; ii < nAnglesToProcess; ii++)
	if (anglesToProcess[ii] > 0 && anglesToProcess[ii] <= (int) header.nz)
	  {
	    stackIn[(anglesToProcess[ii]) - 1].skipMe = false;
	  }
	else
	  {
	    printf ("[WARNING] projection %d not in the views, ignored!\n",
		    anglesToProcess[ii]);
	  }

    }


  for (ii = 0; ii < (unsigned int) header.nz; ii++)
    {
      float angle;
      if (noTiltAngle)
	{
	  angle = 0.0;
	}
      else
	{
	  angle = tiltAngles[ii];
	}



      currSlice = sliceCreate (nx, ny, sliceMode);
      if (!currSlice)
	exitError ("Creating slice for input");
      //get the type of data we are dealing with
      if (mrc_read_slice (currSlice->data.b, fpStack, &header, ii, 'Z'))
	exitError ("Reading slice %d", ii);

      // Convert slice to floats
      if (sliceMode != SLICE_MODE_FLOAT)
	if (sliceNewMode (currSlice, SLICE_MODE_FLOAT) < 0)
	  exitError ("Converting slice to float");

      stackIn[ii].Data =
	matrix (0, nx + padding + padding, 0, ny + padding + padding);
      stackOut[ii].Data =
	matrix (0, nx + padding + padding, 0, ny + padding + padding);
      stackOut[ii].angleRadiant = (angle * MY_PI) / 180.0;
      stackOut[ii].angleDegrees = angle;
      stackOut[ii].TiltNumber = ii;

      // Copy data into array
      for (unsigned int j = 0; j < ny; j++)
	for (unsigned int i = 0; i < nx; i++)
	  {
	    (stackIn[ii].Data)[i + padding][j + padding] =
	      currSlice->data.f[i + j * nx];
	    (stackOut[ii].Data)[i + padding][j + padding] =
	      currSlice->data.f[i + j * nx];
	  }

      FillingPadding (stackIn[ii].Data, nx, ny, padding);
      FillingPadding (stackOut[ii].Data, nx, ny, padding);

    }



/*
********************************************************************
********************************************************************
********************************************************************
*/

  //Get mask if any
  if (gotMaskInput)
    {
      int sliceMaskMode;
      FILE *fpMask;
      if ((fpMask = iiFOpen (maskFn, "rb")) == 0)
	exitError ("could not open input file %s", maskFn);


      MrcHeader maskHeader;
      if (mrc_head_read (fpMask, &maskHeader))
	exitError ("reading header of input file %s", maskFn);

      sliceMaskMode = sliceModeIfReal (maskHeader.mode);
      if (sliceMaskMode < 0)
	exitError
	  ("File mode is %d; only byte, short integer, or real allowed",
	   maskHeader.mode);

      if (header.nz != maskHeader.nz || header.nx != maskHeader.nx
	  || header.ny != maskHeader.ny)
	exitError ("Mask size not compatible with input image size\n");



      for (ii = 0; ii < (unsigned int) header.nz; ii++)
	{
	  Islice *currMaskSlice;
	  currMaskSlice = sliceCreate (nx, ny, sliceMaskMode);
	  if (!currMaskSlice)
	    exitError ("Creating slice for mask");
	  if (mrc_read_slice
	      (currMaskSlice->data.b, fpMask, &maskHeader, ii, 'Z'))
	    exitError ("Reading mask slide %d", ii);

	  // Convert slice to byte
	  if (sliceMaskMode != SLICE_MODE_FLOAT)
	    if (sliceNewMode (currMaskSlice, SLICE_MODE_FLOAT) < 0)
	      exitError ("Converting slice to float");

	  stackMask[ii].Data =
	    matrix (0, nx + padding + padding, 0, ny + padding + padding);

	  // erasing data
	  for (unsigned int j = 0; j < ny + 2 * padding; j++)
	    for (unsigned int i = 0; i < nx + 2 * padding; i++)
	      (stackMask[ii].Data)[i][j] = 0;

	  // Copy data into array
	  for (unsigned int j = 0; j < ny; j++)
	    for (unsigned int i = 0; i < nx; i++)
	      {

		float value = 0.0;
		if (currMaskSlice->data.f[i + j * nx] > 0.000001)
		  value = 1.0;
		else
		  value = 0.0;

		(stackMask[ii].Data)[i + padding][j + padding] = value;
	      }

	}
      iiFClose (fpMask);

    }




  printf ("start processing\n");
  // ********************************************************************
  // ********************** PRE NON-LINEAR ANISOTROPIC DIFFUSION
  // ********************************************************************

  automaticPreLR (stackIn, stackOut, stackMask, nx, ny, header.nz, padding,
		  sigmaList, alphaSigma_List, betaSigma_List, tau_List,
		  Iterations_List, spacingX, spacingY, lambdaE, gotMaskOutput,
		  gotMaskInput);

  // ********************************************************************
  // ****  END  ***************  PRE NON-LINEAR ANISOTROPIC DIFFUSION
  // ********************************************************************
  printf
    ("\n***************\n write image data and close file \n***************\n");

  // DNM: get the new min/max/mean and write header at the end
  outHeader.amin = 1.e37;
  outHeader.amax = -1.e37;
  outHeader.amean = 0.;
  for (ii = 0; ii < (unsigned int) header.nz; ii++)
    {
      currSlice = sliceCreate (nx, ny, SLICE_MODE_FLOAT);
      if (!currSlice)
	exitError ("Creating slice for output");
      for (unsigned int j = 0; j < ny; j++)
	for (unsigned int i = 0; i < nx; i++)
	  {
	    currSlice->data.f[i + j * nx] =
	      (stackOut[ii].Data)[i + padding][j + padding];
	  }

      // Convert if necessary and write slice
      if (sliceMode != SLICE_MODE_FLOAT)
	if (sliceNewMode (currSlice, sliceMode) < 0)
	  exitError ("Converting slice to short");
      if (mrc_write_slice (currSlice->data.b, foutput, &outHeader, ii, 'Z'))
	exitError ("Writing slice %d", ii);

      sliceMMM(currSlice);
      ACCUM_MIN(outHeader.amin, currSlice->min);
      ACCUM_MAX(outHeader.amax, currSlice->max);
      outHeader.amean += currSlice->mean / outHeader.nz;

      sliceFree (currSlice);
      free_matrix (stackIn[ii].Data, 0, nx + padding + padding, 0,
		   ny + padding + padding);
      free_matrix (stackOut[ii].Data, 0, nx + padding + padding, 0,
		   ny + padding + padding);
      if (gotMaskInput)
	free_matrix (stackMask[ii].Data, 0, nx + padding + padding, 0,
		     ny + padding + padding);
    }


  mrc_head_write(foutput, &outHeader);
  free (stackIn);
  free (stackOut);
  if (gotMaskInput)
    free (stackMask);

  iiFClose (fpStack);
  iiFClose (foutput);


}
