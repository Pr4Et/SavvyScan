/*
 * ctffind.cpp - Function extracted from ctffind program of Rohou and Grigorieff
 *
 * $Id$
 */

#include "core_headers.h"
// IMOD: add include
#include "ctffind.h"

//#define threshold_spectrum

const std::string ctffind_version = "4.1.9";
/*
 * Changelog
 * - 4.1.9
 * -- FRC is between a re-normalized version of the amplitude spectrum, to emphasize phase of the Thon rings over their relative amplitudes
 * -- tweaked criteria for "fit resolution"
 * -- tweaked FRC computation
 * -- astigmatism restraint is off by default
 * -- fixed bug affecting astigmatism 
 * -- tweaked background subtraction (thanks Niko!) - helps with noisy VPP spectra
 */

// IMOD: remove class CtffindApp

class ImageCTFComparison
{
public:
	ImageCTFComparison(int wanted_number_of_images, CTF wanted_ctf, float wanted_pixel_size, bool should_find_phase_shift, bool wanted_astigmatism_is_known, float wanted_known_astigmatism, float wanted_known_astigmatism_angle, bool should_fit_defocus_sweep);
	~ImageCTFComparison();
	void SetImage(int wanted_image_number, Image *new_image);
	void SetCTF(CTF new_ctf);
	CTF ReturnCTF();
	bool AstigmatismIsKnown();
	float ReturnKnownAstigmatism();
	float ReturnKnownAstigmatismAngle();
	bool FindPhaseShift();
  void SetupQuickCorrelation();

	int 	number_of_images;
	Image 	*img;		// Usually an amplitude spectrum, or an array of amplitude spectra
  int number_to_correlate;
  double norm_image;
  double image_mean;
  float *azimuths;
  float *spatial_frequency_squared;
  int *addresses;

private:
	CTF		ctf;
	float	pixel_size;
	bool	find_phase_shift;
	bool	astigmatism_is_known;
	float	known_astigmatism;
	float 	known_astigmatism_angle;
	bool 	fit_defocus_sweep;
};

class CurveCTFComparison
{
public:
	float	*curve;	// Usually the 1D rotational average of the amplitude spectrum of an image
	int		number_of_bins;
	float	reciprocal_pixel_size; // In reciprocal pixels
	CTF		ctf;
	bool 	find_phase_shift;
};

ImageCTFComparison::ImageCTFComparison(int wanted_number_of_images, CTF wanted_ctf, float wanted_pixel_size, bool should_find_phase_shift, bool wanted_astigmatism_is_known, float wanted_known_astigmatism, float wanted_known_astigmatism_angle, bool should_fit_defocus_sweep)
{
	MyDebugAssertTrue(wanted_number_of_images >= 0, "Bad wanted number of images: %i\n",wanted_number_of_images);
	number_of_images = wanted_number_of_images;
	img = new Image [wanted_number_of_images];

	ctf = wanted_ctf;
	pixel_size = wanted_pixel_size;
	find_phase_shift = should_find_phase_shift;
	astigmatism_is_known = wanted_astigmatism_is_known;
	known_astigmatism = wanted_known_astigmatism;
	known_astigmatism_angle = wanted_known_astigmatism_angle;
	fit_defocus_sweep = should_fit_defocus_sweep;
  azimuths = NULL;
  spatial_frequency_squared = NULL;
  addresses = NULL;
  number_to_correlate = 0;
}

ImageCTFComparison::~ImageCTFComparison()
{
	for (int image_counter = 0; image_counter < number_of_images; image_counter++)
	{
		img[image_counter].Deallocate();
	}
	delete [] img;
  delete [] azimuths;
  delete [] spatial_frequency_squared;
  delete [] addresses;
}

void ImageCTFComparison::SetImage(int wanted_image_number, Image *new_image)
{
	MyDebugAssertTrue(wanted_image_number >= 0 && wanted_image_number < number_of_images, "Wanted image number (%i) is out of bounds", wanted_image_number);
	img[wanted_image_number].CopyFrom(new_image);
}

void ImageCTFComparison::SetCTF(CTF new_ctf)
{
	ctf = new_ctf;
}

void ImageCTFComparison::SetupQuickCorrelation()
{
  img[0].SetupQuickCorrelationWithCTF(ctf, number_to_correlate, norm_image, image_mean, NULL, NULL, NULL);
  azimuths = new float[number_to_correlate];
  spatial_frequency_squared = new float[number_to_correlate];
  addresses = new int[number_to_correlate];
  img[0].SetupQuickCorrelationWithCTF(ctf, number_to_correlate, norm_image, image_mean, addresses, spatial_frequency_squared, azimuths);
}

CTF ImageCTFComparison::ReturnCTF() { return ctf; }
bool ImageCTFComparison::AstigmatismIsKnown() { return astigmatism_is_known; }
float ImageCTFComparison::ReturnKnownAstigmatism() { return known_astigmatism; }
float ImageCTFComparison::ReturnKnownAstigmatismAngle() { return known_astigmatism_angle; }
bool ImageCTFComparison::FindPhaseShift() { return find_phase_shift; }

static bool inConjGrad = false;

// This is the function which will be minimised
float CtffindObjectiveFunction(void *scoring_parameters, float array_of_values[] )
{
	ImageCTFComparison *comparison_object = reinterpret_cast < ImageCTFComparison *> (scoring_parameters);

	CTF my_ctf = comparison_object->ReturnCTF();
	if (comparison_object->AstigmatismIsKnown())
	{
		MyDebugAssertTrue(comparison_object->ReturnKnownAstigmatism() >= 0.0,"Known asitgmatism must be >= 0.0");
		my_ctf.SetDefocus(array_of_values[0],array_of_values[0] - comparison_object->ReturnKnownAstigmatism(), comparison_object->ReturnKnownAstigmatismAngle());
	}
	else
	{
		my_ctf.SetDefocus(array_of_values[0],array_of_values[1],array_of_values[2]);
	}
	if (comparison_object->FindPhaseShift())
	{
		if (comparison_object->AstigmatismIsKnown())
		{
			my_ctf.SetAdditionalPhaseShift(array_of_values[1]);
		}
		else
		{
			my_ctf.SetAdditionalPhaseShift(array_of_values[3]);
		}
	}

  /*if (inConjGrad)
    MyDebugPrint("(CtffindObjectiveFunction) D1 = %6.2f pxl D2 = %6.2f pxl, PhaseShift = %6.3f rad, Ast = %5.3f rad, Score = %g\n",my_ctf.GetDefocus1(),my_ctf.GetDefocus2(),my_ctf.GetAdditionalPhaseShift(), my_ctf.GetAstigmatismAzimuth(),- - comparison_object->img[0].GetCorrelationWithCTF(my_ctf));*/

	// Evaluate the function
  if (comparison_object->number_to_correlate) 
    return - comparison_object->img[0].QuickCorrelationWithCTF(my_ctf, comparison_object->number_to_correlate, comparison_object->norm_image, comparison_object->image_mean, comparison_object->addresses, 
                                                               comparison_object->spatial_frequency_squared, comparison_object->azimuths);
  else
    return - comparison_object->img[0].GetCorrelationWithCTF(my_ctf);
}

//#pragma GCC push_options
//#pragma GCC optimize ("O0")

// This is the function which will be minimised when dealing with 1D fitting
float CtffindCurveObjectiveFunction(void *scoring_parameters, float array_of_values[] )
{
	CurveCTFComparison *comparison_object = reinterpret_cast < CurveCTFComparison *> (scoring_parameters);

	CTF my_ctf = comparison_object->ctf;
	my_ctf.SetDefocus(array_of_values[0],array_of_values[0],0.0);
	if (comparison_object->find_phase_shift)
	{
		my_ctf.SetAdditionalPhaseShift(array_of_values[1]);
	}

	// Compute the cross-correlation
	double cross_product = 0.0;
	double norm_curve = 0.0;
	double norm_ctf = 0.0;
	int number_of_values = 0;
	int bin_counter;
	float current_spatial_frequency_squared;
	const float lowest_freq = pow(my_ctf.GetLowestFrequencyForFitting(),2);
	const float highest_freq = pow(my_ctf.GetHighestFrequencyForFitting(),2);
	float current_ctf_value;

	for ( bin_counter = 0 ; bin_counter < comparison_object->number_of_bins; bin_counter ++ )
	{
		current_spatial_frequency_squared = pow(float(bin_counter)*comparison_object->reciprocal_pixel_size,2);
		if (current_spatial_frequency_squared > lowest_freq && current_spatial_frequency_squared < highest_freq)
		{
			current_ctf_value = fabsf(my_ctf.Evaluate(current_spatial_frequency_squared,0.0));
			MyDebugAssertTrue(current_ctf_value >= -1.0 && current_ctf_value <= 1.0,"Bad ctf value: %f",current_ctf_value);
			number_of_values++;
			cross_product += comparison_object->curve[bin_counter] * current_ctf_value;
			norm_curve += pow(comparison_object->curve[bin_counter],2);
			norm_ctf += pow(current_ctf_value,2);
		}
	}

	MyDebugAssertTrue(norm_ctf > 0.0,"Bad norm_ctf: %f\n", norm_ctf);
	MyDebugAssertTrue(norm_curve > 0.0,"Bad norm_curve: %f\n", norm_curve);

	//MyDebugPrint("(CtffindCurveObjectiveFunction) D1 = %6.2f , PhaseShift = %6.3f , Score = %g",array_of_values[0], array_of_values[1], - cross_product / sqrtf(norm_ctf * norm_curve));

	// Note, we are not properly normalizing the cross correlation coefficient. For our
	// purposes this should be OK, since the average power of the theoretical CTF should not
	// change much with defocus. At least I hope so.
	return - cross_product / sqrtf(norm_ctf * norm_curve);


}

//#pragma GCC pop_options
// IMOD: change two to bool for error handling
float FindRotationalAlignmentBetweenTwoStacksOfImages(Image *self, Image *other_image, int number_of_images, float search_half_range, float search_step_size, float minimum_radius, float maximum_radius);
void ComputeImagesWithNumberOfExtremaAndCTFValues(CTF *ctf, Image *number_of_extrema, Image *ctf_values);
int ReturnSpectrumBinNumber(int number_of_bins, float number_of_extrema_profile[], Image *number_of_extrema, long address, Image *ctf_values, float ctf_values_profile[]);
bool ComputeRotationalAverageOfPowerSpectrum( Image *spectrum, CTF *ctf, Image *number_of_extrema, Image *ctf_values, int number_of_bins, double spatial_frequency[], double average[], double average_fit[], double average_renormalized[], float number_of_extrema_profile[], float ctf_values_profile[]);
void OverlayCTF( Image *spectrum, CTF *ctf);
void ComputeFRCBetween1DSpectrumAndFit( int number_of_bins, double average[], double fit[], float number_of_extrema_profile[], double frc[], double frc_sigma[], int first_fit_bin);
bool RescaleSpectrumAndRotationalAverage( Image *spectrum, Image *number_of_extrema, Image *ctf_values, int number_of_bins, double spatial_frequency[], double average[], double average_fit[], float number_of_extrema_profile[], float ctf_values_profile[], int last_bin_without_aliasing, int last_bin_with_good_fit );
void Renormalize1DSpectrumForFRC( int number_of_bins, double average[], double fit[], float number_of_extrema_profile[]);

// IMOD: Exported functions for setting the print function and slice write function
void ctffindSetSliceWriteFunc(WriteSliceType func)
{
  internalSetWriteSliceFunc(func);
}

void ctffindSetPrintFunc(CharArgType func)
{
  internalSetPrintFunc(func);
}

// IMOD: Replace CtffindApp::DoInteractiveUserInput() with this function, remove all user input and
// assignment of input values to variables, replace with assignment of param values to variables   

bool ctffind(CtffindParams &params, float *spectrumArray, int nxDimIn, 
             float *results_array, float **rotationalAvgOut, float **normalizedAvgOut,
             float **fitCurveOut, int &numPointsOut, float &lastBinFreqOut)
{

	// Arguments for this job

  float         pixel_size_of_input_image = params.pixel_size_of_input_image;
  const float     acceleration_voltage = params.acceleration_voltage;
  const float         spherical_aberration = params.spherical_aberration;
  const float     amplitude_contrast = params.amplitude_contrast;
  const int           box_size = params.box_size;
  const float     minimum_resolution = params.minimum_resolution;
  const float         maximum_resolution = params.maximum_resolution;
  const float         minimum_defocus = params.minimum_defocus;
  const float         maximum_defocus = params.maximum_defocus;
  const float         defocus_search_step = params.defocus_search_step;
  const bool      slower_search = params.slower_search;
  const float         astigmatism_tolerance = params.astigmatism_tolerance;
  const bool        find_additional_phase_shift = params.find_additional_phase_shift;
  const float     minimum_additional_phase_shift = params.minimum_additional_phase_shift;
  const float     maximum_additional_phase_shift = params.maximum_additional_phase_shift;
  const float     additional_phase_shift_search_step = params.additional_phase_shift_search_step;
  const bool      astigmatism_is_known = params.astigmatism_is_known;
  const float     known_astigmatism = params.known_astigmatism;
  const float     known_astigmatism_angle = params.known_astigmatism_angle;
  const bool      compute_extra_stats = params.compute_extra_stats;

	float 				pixel_size_for_fitting = params.pixel_size_of_input_image;


	// Maybe the user wants to hold the phase shift value (which they can do by giving the same value for min and max)
	const bool			fixed_additional_phase_shift = fabs(maximum_additional_phase_shift - minimum_additional_phase_shift) < 0.01;

	// This could become a user-supplied parameter later - for now only for developers / expert users
	const bool			follow_1d_search_with_local_2D_brute_force = false;

	/*
	 *  Scoring function
	 */
	float MyFunction(float []);

  // IMOD: remove various unneeded variables for movies and multiple micrographs
	// Other variables
	Image				*average_spectrum = new Image();
	Image				*average_spectrum_masked = new Image();
	Image				*current_power_spectrum = new Image();
	Image				*temp_image = new Image();
	CTF					current_ctf;
	float				average, sigma;
	int					convolution_box_size;
	ImageCTFComparison	*comparison_object_2D;
	CurveCTFComparison	comparison_object_1D;
	float 				estimated_astigmatism_angle;
	float				bf_halfrange[4];
	float				bf_midpoint[4];
	float				bf_stepsize[4];
	float				cg_starting_point[4];
	float				cg_accuracy[4];
	int 				number_of_search_dimensions;
	BruteForceSearch   	*brute_force_search;
	int					counter;
	ConjugateGradient   *conjugate_gradient_minimizer;
	int					number_of_bins_in_1d_spectra;
	Curve				number_of_averaged_pixels;
	Curve				rotational_average;
	Image				*number_of_extrema_image = new Image();
	Image				*ctf_values_image = new Image();
	double				*rotational_average_astig = NULL;
	double				*rotational_average_astig_renormalized = NULL;
	double				*spatial_frequency = NULL;
	double				*spatial_frequency_in_reciprocal_angstroms = NULL;
	double				*rotational_average_astig_fit = NULL;
	float				*number_of_extrema_profile = NULL;
	float				*ctf_values_profile = NULL;
	double				*fit_frc = NULL;
	double				*fit_frc_sigma = NULL;
	int					last_bin_with_good_fit;
	float				best_score_after_initial_phase;
	int					last_bin_without_aliasing;
	float				intermediate_resolution = 5.0;

  // IMOD: set these variables
  bool is_running_locally = true;
  int ix, iy;

  // IMOD: print the messages
	// Some argument checking
	if (minimum_resolution < maximum_resolution)
    {
      wxPrintf("Error: Minimum resolution (%f) higher than maximum resolution (%f).", minimum_resolution,maximum_resolution);
      return false;
    }
	if (minimum_defocus > maximum_defocus)
    {
      wxPrintf("Error: Minimum defocus must be less than maximum defocus.");
      return false;
    }

  // IMOD: remove setting up loops, preparation of output files, gain reference
  double wallStart = ctfWallTime();

	// Prepare the average spectrum image
	average_spectrum->Allocate(box_size,box_size,true);

  // IMOD: Copy the input spectrum into the image array
  for (iy = 0; iy < box_size; iy++)
    memcpy(average_spectrum->real_values + iy * (box_size + average_spectrum->padding_jump_value),
           spectrumArray + iy * nxDimIn, box_size * sizeof(float));


  // IMOD: remove loop on input images, change this from current_power_spectrum to average_spectrum
  // Set origin of amplitude spectrum to 0.0
  average_spectrum->real_values[average_spectrum->ReturnReal1DAddressFromPhysicalCoord(average_spectrum->physical_address_of_box_center_x,average_spectrum->physical_address_of_box_center_y,average_spectrum->physical_address_of_box_center_z)] = 0.0;


  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_before_bg_sub.mrc",1);
  // IMOD: remove all processing to get clipped power spectrum

  // Filter the amplitude spectrum, remove background
  // IMOD: do this unconditionally
  // Try to weaken cross artefacts
  average_spectrum->ComputeAverageAndSigmaOfValuesInSpectrum(float(average_spectrum->logical_x_dimension)*pixel_size_for_fitting/minimum_resolution,float(average_spectrum->logical_x_dimension),average,sigma,12);
  average_spectrum->DivideByConstant(sigma);
  average_spectrum->SetMaximumValueOnCentralCross(average/sigma+10.0);
  //			average_spectrum_masked->CopyFrom(average_spectrum);

  //average_spectrum->QuickAndDirtyWriteSlice("dbg_average_spectrum_before_conv.mrc",1);
  //wxPrintf("max central cross time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  // Compute low-pass filtered version of the spectrum
  convolution_box_size = int( float(average_spectrum->logical_x_dimension) * pixel_size_for_fitting / minimum_resolution * sqrt(2.0) );
  if (IsEven(convolution_box_size)) convolution_box_size++;
  current_power_spectrum->Allocate(average_spectrum->logical_x_dimension,average_spectrum->logical_y_dimension,true);
  current_power_spectrum->SetToConstant(0.0); // According to valgrind, this avoid potential problems later on.
  average_spectrum->SpectrumBoxConvolution(current_power_spectrum,convolution_box_size,float(average_spectrum->logical_x_dimension)*pixel_size_for_fitting/minimum_resolution);
  //wxPrintf("convolution box %d time = %.3f\n", convolution_box_size, ctfWallTime() - wallStart); wallStart = ctfWallTime();

  //current_power_spectrum->QuickAndDirtyWriteSlice("dbg_spec_convoluted.mrc",1);

  // Subtract low-pass-filtered spectrum from the spectrum. This should remove the background slope.
  average_spectrum->SubtractImage(current_power_spectrum);

  //			average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_before_thresh.mrc",1);

  // Threshold high values
  average_spectrum->SetMaximumValue(average_spectrum->ReturnMaximumValue(3,3));

  //			convolution_box_size = int( float(average_spectrum->logical_x_dimension) * pixel_size_for_fitting / minimum_resolution / sqrt(2.0) );
  //			if (IsEven(convolution_box_size)) convolution_box_size++;
  //			current_power_spectrum->SetToConstant(0.0); // According to valgrind, this avoid potential problems later on.
  //			average_spectrum_masked->SpectrumBoxConvolution(current_power_spectrum,convolution_box_size,float(average_spectrum_masked->logical_x_dimension)*pixel_size_for_fitting/minimum_resolution);
  //			average_spectrum_masked->SubtractImage(current_power_spectrum);
  //			average_spectrum_masked->SetMaximumValue(average_spectrum_masked->ReturnMaximumValue(3,3));
  //wxPrintf("subtract & max time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  average_spectrum_masked->CopyFrom(average_spectrum);
  average_spectrum_masked->CosineMask(float(average_spectrum_masked->logical_x_dimension)*pixel_size_for_fitting/std::max(maximum_resolution, 8.0f),float(average_spectrum_masked->logical_x_dimension)*pixel_size_for_fitting/std::max(maximum_resolution, 4.0f), true);
  //      average_spectrum_masked->QuickAndDirtyWriteSlice("dbg_spec_before_thresh.mrc",1);
  //			average_spectrum_masked->CorrectSinc();
  //			average_spectrum_masked->CorrectSinc(float(average_spectrum_masked->logical_x_dimension)*pixel_size_for_fitting/std::max(maximum_resolution, 8.0f), 0.5, true, 0.0);
  //wxPrintf("cosine mask time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  // IMOD: that was end of input loop producing one spectrum
  /*
   *
   *
   * We now have a spectrum which we can use to fit CTFs
   *
   *
   */

  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec.mrc",1);


#ifdef threshold_spectrum
  wxPrintf("DEBUG: thresholding spectrum\n");
  for (counter = 0; counter < average_spectrum->real_memory_allocated; counter ++ )
		{
			average_spectrum->real_values[counter] = std::max(average_spectrum->real_values[counter], -0.0f);
			average_spectrum->real_values[counter] = std::min(average_spectrum->real_values[counter], 1.0f);
		}
  average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_thr.mrc",1);
#endif

  // Set up the CTF object
  current_ctf.Init(acceleration_voltage,spherical_aberration,amplitude_contrast,minimum_defocus,minimum_defocus,0.0,1.0/minimum_resolution,1.0/std::max(maximum_resolution,intermediate_resolution),astigmatism_tolerance,pixel_size_for_fitting,minimum_additional_phase_shift);
  current_ctf.SetDefocus(minimum_defocus/pixel_size_for_fitting,minimum_defocus/pixel_size_for_fitting,0.0);
  current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);

  // Set up the comparison object
  // DNM: Do not tell it to find phase if it is fixed
  comparison_object_2D = new ImageCTFComparison(1,current_ctf,pixel_size_for_fitting,find_additional_phase_shift && ! fixed_additional_phase_shift, astigmatism_is_known, known_astigmatism / pixel_size_for_fitting, known_astigmatism_angle / 180.0 * PI, false);
  comparison_object_2D->SetImage(0,average_spectrum_masked);
  comparison_object_2D->SetupQuickCorrelation();

  // IMOD: remove print

  // Let's look for the astigmatism angle first
  if (astigmatism_is_known)
		{
			estimated_astigmatism_angle = known_astigmatism_angle;
		}
  else
		{
			temp_image->CopyFrom(average_spectrum);
			temp_image->ApplyMirrorAlongY();
			//temp_image->QuickAndDirtyWriteSlice("dbg_spec_y.mrc",1);
			estimated_astigmatism_angle = 0.5 * FindRotationalAlignmentBetweenTwoStacksOfImages(average_spectrum,temp_image,1,90.0,5.0,pixel_size_for_fitting/minimum_resolution,pixel_size_for_fitting/std::max(maximum_resolution,intermediate_resolution));
		}

  //wxPrintf("Estimated astigmatism angle = %f degrees\n", estimated_astigmatism_angle);
  //wxPrintf("astig est time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();



  /*
   * Initial brute-force search, in 1D (fast, but not as accurate)
   */
  if (!slower_search)
		{

			// 1D rotational average
			number_of_bins_in_1d_spectra = int(ceil(average_spectrum_masked->ReturnMaximumDiagonalRadius()));
			rotational_average.SetupXAxis(0.0,float(number_of_bins_in_1d_spectra) * average_spectrum_masked->fourier_voxel_size_x,number_of_bins_in_1d_spectra);
			number_of_averaged_pixels = rotational_average;
			average_spectrum_masked->Compute1DRotationalAverage(rotational_average,number_of_averaged_pixels,true);



			comparison_object_1D.ctf = current_ctf;
			comparison_object_1D.curve = new float[number_of_bins_in_1d_spectra];
			for (counter=0; counter < number_of_bins_in_1d_spectra; counter++)
        {
          comparison_object_1D.curve[counter] = rotational_average.data_y[counter];
        }

      // DNM: Do not find phase if it is fixed
			comparison_object_1D.find_phase_shift = find_additional_phase_shift && ! fixed_additional_phase_shift;
			comparison_object_1D.number_of_bins = number_of_bins_in_1d_spectra;
			comparison_object_1D.reciprocal_pixel_size = average_spectrum_masked->fourier_voxel_size_x;

			// We can now look for the defocus value
			bf_halfrange[0] = 0.5 * (maximum_defocus - minimum_defocus) / pixel_size_for_fitting;
			bf_halfrange[1] = 0.5 * (maximum_additional_phase_shift - minimum_additional_phase_shift);

			bf_midpoint[0] = minimum_defocus / pixel_size_for_fitting + bf_halfrange[0];
			bf_midpoint[1] = minimum_additional_phase_shift + bf_halfrange[1];

			bf_stepsize[0] = defocus_search_step / pixel_size_for_fitting;
			bf_stepsize[1] = additional_phase_shift_search_step;

			if (find_additional_phase_shift && ! fixed_additional_phase_shift)
        {
          number_of_search_dimensions = 2;
        }
			else
        {
          number_of_search_dimensions = 1;
        }

			// Actually run the BF search
			brute_force_search = new BruteForceSearch();
			brute_force_search->Init(&CtffindCurveObjectiveFunction,&comparison_object_1D,number_of_search_dimensions,bf_midpoint,bf_halfrange,bf_stepsize,false,false);
			brute_force_search->Run();

			/*			
              wxPrintf("After 1D brute\n");
              wxPrintf("      DFMID1      DFMID2      ANGAST          CC\n");
              wxPrintf("%12.2f%12.2f%12.2f%12.5f\n",brute_force_search->GetBestValue(0),brute_force_search->GetBestValue(0),0.0,brute_force_search->GetBestScore());
              wxPrintf("%12.2f%12.2f%12.2f%12.5f\n",brute_force_search->GetBestValue(0)*pixel_size_for_fitting,brute_force_search->GetBestValue(0)*pixel_size_for_fitting,0.0,brute_force_search->GetBestScore());
			*/

			// We can now do a local optimization
			// The end point of the BF search is the beginning of the CG search
			for (counter=0;counter<number_of_search_dimensions;counter++)
        {
          cg_starting_point[counter] = brute_force_search->GetBestValue(counter);
        }
			cg_accuracy[0] = 100.0;
			cg_accuracy[1] = 0.05;
			conjugate_gradient_minimizer = new ConjugateGradient();
			conjugate_gradient_minimizer->Init(&CtffindCurveObjectiveFunction,&comparison_object_1D,number_of_search_dimensions,cg_starting_point,cg_accuracy);
			conjugate_gradient_minimizer->Run();
			for (counter=0;counter<number_of_search_dimensions;counter++)
        {
          cg_starting_point[counter] = conjugate_gradient_minimizer->GetBestValue(counter);
        }
			current_ctf.SetDefocus(cg_starting_point[0], cg_starting_point[0],estimated_astigmatism_angle / 180.0 * PI);
			if (find_additional_phase_shift)
        {
          if (fixed_additional_phase_shift)
            {
              current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
            }
          else
            {
              current_ctf.SetAdditionalPhaseShift(cg_starting_point[1]);
            }
        }

			// Remember the best score so far
			best_score_after_initial_phase = - conjugate_gradient_minimizer->GetBestScore();


			// Cleanup
			delete conjugate_gradient_minimizer;
			delete brute_force_search;
			delete [] comparison_object_1D.curve;

		} // end of the fast search over the 1D function


  /*
   * Brute-force search over the 2D scoring function.
   * This is either the first search we are doing, or just a refinement
   * starting from the result of the 1D search
   */
  if (slower_search || (!slower_search && follow_1d_search_with_local_2D_brute_force))
		{
			// Setup the parameters for the brute force search

			if (slower_search) // This is the first search we are doing - scan the entire range the user specified
        {
          if (astigmatism_is_known)
            {
              bf_halfrange[0] = 0.5 * (maximum_defocus-minimum_defocus)/pixel_size_for_fitting;
              bf_halfrange[1] = 0.5 * (maximum_additional_phase_shift-minimum_additional_phase_shift);

              bf_midpoint[0] = minimum_defocus/pixel_size_for_fitting + bf_halfrange[0];
              bf_midpoint[1] = minimum_additional_phase_shift + bf_halfrange[3];

              bf_stepsize[0] = defocus_search_step/pixel_size_for_fitting;
              bf_stepsize[1] = additional_phase_shift_search_step;

              if (find_additional_phase_shift && ! fixed_additional_phase_shift)
                {
                  number_of_search_dimensions = 2;
                }
              else
                {
                  number_of_search_dimensions = 1;
                }
            }
          else
            {
              bf_halfrange[0] = 0.5 * (maximum_defocus-minimum_defocus)/pixel_size_for_fitting;
              bf_halfrange[1] = bf_halfrange[0];
              bf_halfrange[2] = 0.0;
              bf_halfrange[3] = 0.5 * (maximum_additional_phase_shift-minimum_additional_phase_shift);

              bf_midpoint[0] = minimum_defocus/pixel_size_for_fitting + bf_halfrange[0];
              bf_midpoint[1] = bf_midpoint[0];
              bf_midpoint[2] = estimated_astigmatism_angle / 180.0 * PI;
              bf_midpoint[3] = minimum_additional_phase_shift + bf_halfrange[3];

              bf_stepsize[0] = defocus_search_step/pixel_size_for_fitting;
              bf_stepsize[1] = bf_stepsize[0];
              bf_stepsize[2] = 0.0;
              bf_stepsize[3] = additional_phase_shift_search_step;

              if (find_additional_phase_shift && ! fixed_additional_phase_shift)
                {
                  number_of_search_dimensions = 4;
                }
              else
                {
                  number_of_search_dimensions = 3;
                }
            }
        }
			else // we will do a brute-force search near the result of the search over the 1D objective function
        {
          if (astigmatism_is_known)
            {

              bf_midpoint[0] = current_ctf.GetDefocus1();
              bf_midpoint[1] = current_ctf.GetAdditionalPhaseShift();

              bf_stepsize[0] = defocus_search_step/pixel_size_for_fitting;
              bf_stepsize[1] = additional_phase_shift_search_step;

              bf_halfrange[0] = 2.0 * defocus_search_step/pixel_size_for_fitting + 0.1;
              bf_halfrange[1] = 2.0 * additional_phase_shift_search_step + 0.01;


              if (find_additional_phase_shift && ! fixed_additional_phase_shift)
                {
                  number_of_search_dimensions = 2;
                }
              else
                {
                  number_of_search_dimensions = 1;
                }
            }
          else
            {

              bf_midpoint[0] = current_ctf.GetDefocus1();
              bf_midpoint[1] = current_ctf.GetDefocus2();
              bf_midpoint[2] = current_ctf.GetAstigmatismAzimuth();
              bf_midpoint[3] = minimum_additional_phase_shift + bf_halfrange[3];

              bf_stepsize[0] = defocus_search_step/pixel_size_for_fitting;
              bf_stepsize[1] = bf_stepsize[0];
              bf_stepsize[2] = 0.0;
              bf_stepsize[3] = additional_phase_shift_search_step;

              if (astigmatism_tolerance > 0)
                {
                  bf_halfrange[0] = 2.0 * astigmatism_tolerance/pixel_size_for_fitting + 0.1;
                }
              else
                {
                  bf_halfrange[0] = 2.0 * defocus_search_step/pixel_size_for_fitting + 0.1;
                }
              bf_halfrange[1] = bf_halfrange[0];
              bf_halfrange[2] = 0.0;
              bf_halfrange[3] = 2.0 * additional_phase_shift_search_step + 0.01;

              if (find_additional_phase_shift && ! fixed_additional_phase_shift)
                {
                  number_of_search_dimensions = 4;
                }
              else
                {
                  number_of_search_dimensions = 3;
                }
            }
        }

      // DNM: Do one-time set of phase shift for fixed value
      if (find_additional_phase_shift && fixed_additional_phase_shift)
        {
          current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
        }
      
			// Actually run the BF search (we run a local minimizer at every grid point only if this is a refinement search following 1D search (otherwise the full brute-force search would get too long)
			brute_force_search = new BruteForceSearch();
			brute_force_search->Init(&CtffindObjectiveFunction,comparison_object_2D,number_of_search_dimensions,bf_midpoint,bf_halfrange,bf_stepsize,!slower_search,is_running_locally);
			brute_force_search->Run();

			// The end point of the BF search is the beginning of the CG search
			for (counter=0;counter<number_of_search_dimensions;counter++)
        {
          cg_starting_point[counter] = brute_force_search->GetBestValue(counter);
        }

			//
			if (astigmatism_is_known)
        {
          current_ctf.SetDefocus(cg_starting_point[0],cg_starting_point[0] - known_astigmatism / pixel_size_for_fitting, known_astigmatism_angle / 180.0 * PI);
          if (find_additional_phase_shift)
            {
              if (fixed_additional_phase_shift)
                {
                  current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
                }
              else
                {
                  current_ctf.SetAdditionalPhaseShift(cg_starting_point[1]);
                }
            }
        }
			else
        {
          current_ctf.SetDefocus(cg_starting_point[0],cg_starting_point[1],cg_starting_point[2]);
          if (find_additional_phase_shift)
            {
              if (fixed_additional_phase_shift)
                {
                  current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
                }
              else
                {
                  current_ctf.SetAdditionalPhaseShift(cg_starting_point[3]);
                }
            }
        }

			current_ctf.EnforceConvention();

			// Remember the best score so far
			best_score_after_initial_phase = - brute_force_search->GetBestScore();

			delete brute_force_search;
		}

  // IMOD: remove old-school output
  /*wxPrintf("      DFMID1      DFMID2      ANGAST          CC\n");
    wxPrintf("%12.2f%12.2f%12.2f%12.5f\n",current_ctf.GetDefocus1()*pixel_size_for_fitting,current_ctf.GetDefocus2()*pixel_size_for_fitting,current_ctf.GetAstigmatismAzimuth()*180.0/PI,best_score_after_initial_phase);*/

  //wxPrintf("bf time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();
  /*
   * Set up the conjugate gradient minimization of the 2D scoring function
   */
  if (astigmatism_is_known)
		{
			cg_starting_point[0] = current_ctf.GetDefocus1();
			if (find_additional_phase_shift) cg_starting_point[1] = current_ctf.GetAdditionalPhaseShift();
			if (find_additional_phase_shift && ! fixed_additional_phase_shift)
        {
          number_of_search_dimensions = 2;
        }
			else
        {
          number_of_search_dimensions = 1;
        }
			cg_accuracy[0] = 100.0;
			cg_accuracy[1] = 0.05;
		}
  else
		{
			cg_accuracy[0] = 100.0;
			cg_accuracy[1] = 100.0;
			cg_accuracy[2] = 0.025;
			cg_accuracy[3] = 0.05;
			cg_starting_point[0] = current_ctf.GetDefocus1();
			cg_starting_point[1] = current_ctf.GetDefocus2();
      if (slower_search || (!slower_search && follow_1d_search_with_local_2D_brute_force))
        {
          // we did a search against the 2D power spectrum so we have a better estimate
          // of the astigmatism angle in the CTF object
          cg_starting_point[2] = current_ctf.GetAstigmatismAzimuth();
        }
      else
        {
          // all we have right now is the guessed astigmatism angle from the mirror
          // trick before any CTF fitting was even tried
          cg_starting_point[2] = estimated_astigmatism_angle / 180.0 * PI;
        }

			if (find_additional_phase_shift) cg_starting_point[3] = current_ctf.GetAdditionalPhaseShift();
			if (find_additional_phase_shift && ! fixed_additional_phase_shift)
        {
          number_of_search_dimensions = 4;
        }
			else
        {
          number_of_search_dimensions = 3;
        }
		}
  // CG minimization
  inConjGrad = true;
  comparison_object_2D->SetCTF(current_ctf);
  conjugate_gradient_minimizer = new ConjugateGradient();
  conjugate_gradient_minimizer->Init(&CtffindObjectiveFunction,comparison_object_2D,number_of_search_dimensions,cg_starting_point,cg_accuracy);
  current_ctf.Init(acceleration_voltage,spherical_aberration,amplitude_contrast,minimum_defocus,minimum_defocus,0.0,1.0/minimum_resolution,1.0/maximum_resolution,astigmatism_tolerance,pixel_size_for_fitting,minimum_additional_phase_shift);
  conjugate_gradient_minimizer->Run();
  //wxPrintf("CG time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  // Remember the results of the refinement
  for (counter=0;counter<number_of_search_dimensions;counter++)
		{
			cg_starting_point[counter] = conjugate_gradient_minimizer->GetBestValue(counter);
		}
  if (astigmatism_is_known)
		{
			current_ctf.SetDefocus(cg_starting_point[0],cg_starting_point[0] - known_astigmatism / pixel_size_for_fitting, known_astigmatism_angle / 180.0 * PI);
			if (find_additional_phase_shift)
        {
          if (fixed_additional_phase_shift)
            {
              current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
            }
          else
            {
              current_ctf.SetAdditionalPhaseShift(cg_starting_point[1]);
            }
        }
		}
  else
		{
			current_ctf.SetDefocus(cg_starting_point[0],cg_starting_point[1],cg_starting_point[2]);
			if (find_additional_phase_shift)
        {
          if (fixed_additional_phase_shift)
            {
              current_ctf.SetAdditionalPhaseShift(minimum_additional_phase_shift);
            }
          else
            {
              current_ctf.SetAdditionalPhaseShift(cg_starting_point[3]);
            }
        }
		}
  current_ctf.EnforceConvention();
  /*wxPrintf("%12.2f%12.2f%12.2f%12.5f   Final Values\n",current_ctf.GetDefocus1()*pixel_size_for_fitting,current_ctf.GetDefocus2()*pixel_size_for_fitting,current_ctf.GetAstigmatismAzimuth()*180.0/PI,-conjugate_gradient_minimizer->GetBestScore());*/

  // IMOD: remove old-school output and setting current_output_location

  // Generate diagnostic image
  //average_spectrum.QuickAndDirtyWriteSlice("dbg_spec_diag_start.mrc",1);
  average_spectrum->AddConstant(-1.0 * average_spectrum->ReturnAverageOfRealValuesOnEdges());

  /*
   *  Attempt some renormalisations - we want to do this over a range not affected by the central peak or strong Thon rings,
   *  so as to emphasize the "regular" Thon rings
   */
  float start_zero = sqrtf(current_ctf.ReturnSquaredSpatialFrequencyOfAZero(3,0.0));
  float finish_zero = sqrtf(current_ctf.ReturnSquaredSpatialFrequencyOfAZero(4,0.0));
  float normalization_radius_min = start_zero * average_spectrum->logical_x_dimension;
  float normalization_radius_max = finish_zero * average_spectrum->logical_x_dimension;

  if (start_zero > current_ctf.GetHighestFrequencyForFitting() || start_zero < current_ctf.GetLowestFrequencyForFitting() || finish_zero > current_ctf.GetHighestFrequencyForFitting() || finish_zero < current_ctf.GetLowestFrequencyForFitting())
		{
			normalization_radius_max = current_ctf.GetHighestFrequencyForFitting() * average_spectrum->logical_x_dimension;
			normalization_radius_min = std::max(0.5f * normalization_radius_max, current_ctf.GetLowestFrequencyForFitting() * average_spectrum->logical_x_dimension);
		}

  MyDebugAssertTrue(normalization_radius_max > normalization_radius_min,"Bad values for min (%f) and max (%f) normalization radii\n");

  if (normalization_radius_max - normalization_radius_min > 2.0)
		{
			average_spectrum->ComputeAverageAndSigmaOfValuesInSpectrum(	normalization_radius_min,
                                                                  normalization_radius_max,
                                                                  average,sigma);
			average_spectrum->CircleMask(5.0,true);
			average_spectrum->SetMaximumValueOnCentralCross(average);
			average_spectrum->SetMinimumAndMaximumValues(average - 4.0 * sigma, average + 4.0 * sigma);
			average_spectrum->ComputeAverageAndSigmaOfValuesInSpectrum(	normalization_radius_min,
                                                                  normalization_radius_max,
                                                                  average,sigma);
			average_spectrum->AddConstant(-1.0 * average);
			average_spectrum->MultiplyByConstant(1.0 / sigma);
			average_spectrum->AddConstant(average);
		}
  //wxPrintf("renorm time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_diag_1.mrc",1);

  // 1D rotational average
  number_of_bins_in_1d_spectra = int(ceil(average_spectrum->ReturnMaximumDiagonalRadius()));
  rotational_average.SetupXAxis(0.0,float(number_of_bins_in_1d_spectra) * average_spectrum->fourier_voxel_size_x ,number_of_bins_in_1d_spectra);
  rotational_average.ZeroYData();
  //number_of_averaged_pixels.ZeroYData();
  number_of_averaged_pixels = rotational_average;
  average_spectrum->Compute1DRotationalAverage(rotational_average,number_of_averaged_pixels,true);
  //wxPrintf("regular rotavg time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  // Rotational average, taking astigmatism into account
  if (compute_extra_stats)
		{
			number_of_extrema_image->Allocate(average_spectrum->logical_x_dimension,average_spectrum->logical_y_dimension,true);
			ctf_values_image->Allocate(average_spectrum->logical_x_dimension,average_spectrum->logical_y_dimension,true);
			spatial_frequency 						= new double[number_of_bins_in_1d_spectra];
			rotational_average_astig 				= new double[number_of_bins_in_1d_spectra];
			rotational_average_astig_renormalized	= new double[number_of_bins_in_1d_spectra];
			rotational_average_astig_fit			= new double[number_of_bins_in_1d_spectra];
			number_of_extrema_profile 				= new float[number_of_bins_in_1d_spectra];
			ctf_values_profile 						= new float[number_of_bins_in_1d_spectra];
			fit_frc									= new double[number_of_bins_in_1d_spectra];
			fit_frc_sigma							= new double[number_of_bins_in_1d_spectra];
			ComputeImagesWithNumberOfExtremaAndCTFValues(&current_ctf, number_of_extrema_image, ctf_values_image);
			//number_of_extrema_image->QuickAndDirtyWriteSlice("dbg_num_extrema.mrc",1);
			//ctf_values_image->QuickAndDirtyWriteSlice("dbg_ctf_values.mrc",1);
      // IMOD: returnif error
			if (!ComputeRotationalAverageOfPowerSpectrum(average_spectrum, &current_ctf, number_of_extrema_image, ctf_values_image, number_of_bins_in_1d_spectra, spatial_frequency, rotational_average_astig, rotational_average_astig_fit, rotational_average_astig_renormalized, number_of_extrema_profile, ctf_values_profile))
        return false;

			// Here, do FRC
			int first_fit_bin = 0;
			for (int bin_counter = number_of_bins_in_1d_spectra - 1; bin_counter >= 0; bin_counter -- )
        {
          if (spatial_frequency[bin_counter] >= current_ctf.GetLowestFrequencyForFitting()) first_fit_bin = bin_counter;
        }
			ComputeFRCBetween1DSpectrumAndFit(number_of_bins_in_1d_spectra,rotational_average_astig_renormalized,rotational_average_astig_fit,number_of_extrema_profile,fit_frc,fit_frc_sigma,first_fit_bin);

			// At what bin does CTF aliasing become problematic?
			last_bin_without_aliasing = 0;
			int location_of_previous_extremum = 0;
			for (counter=1;counter<number_of_bins_in_1d_spectra;counter++)
        {
          if (number_of_extrema_profile[counter]-number_of_extrema_profile[counter-1] >= 0.9)
            {
              // We just reached a new extremum
              if (counter-location_of_previous_extremum < 4)
                {
                  last_bin_without_aliasing = location_of_previous_extremum;
                  break;
                }
              location_of_previous_extremum = counter;
            }
        }
			/* IMOD
         if (is_running_locally && old_school_input && last_bin_without_aliasing != 0)
         {
         wxPrintf("CTF aliasing apparent from %0.1f Angstroms\n",pixel_size_for_fitting / spatial_frequency[last_bin_without_aliasing]);
         }*/
		}
  //wxPrintf("extra stat time = %.3f\n", ctfWallTime() - wallStart); wallStart = ctfWallTime();

  //average_spectrum.QuickAndDirtyWriteSlice("dbg_spec_diag_2.mrc",1);

  // Until what frequency were CTF rings detected?
  if (compute_extra_stats)
		{
			static float low_threshold = 0.1;
			static float frc_significance_threshold = 0.5; // In analogy to the usual criterion when comparing experimental results to the atomic model
			static float high_threshold = 0.66;
			bool at_last_bin_with_good_fit;
			int number_of_bins_above_low_threshold = 0;
			int number_of_bins_above_significance_threshold = 0;
			int number_of_bins_above_high_threshold = 0;

      // Fix for IMOD: keep track of last unique value and fix tests at end
      int last_bin_with_unique_value = 0;
			int first_bin_to_check = int(sqrtf(current_ctf.ReturnSquaredSpatialFrequencyOfAZero(1,0.0))*average_spectrum->logical_x_dimension);
			//wxPrintf("Will only check from bin %i of %i onwards\n", first_bin_to_check, number_of_bins_in_1d_spectra);
			last_bin_with_good_fit = -1;
			for (counter=first_bin_to_check;counter<number_of_bins_in_1d_spectra;counter++)
        {
          //wxPrintf("On bin %i, fit_frc = %f, rot averate astig = %f\n", counter, fit_frc[counter], rotational_average_astig[counter]);
          at_last_bin_with_good_fit = ( (number_of_bins_above_low_threshold          > 3) &&   ( fit_frc[counter] < low_threshold)  )
            ||
            ( (number_of_bins_above_high_threshold         > 3) &&   ( fit_frc[counter] < frc_significance_threshold) );
          if (at_last_bin_with_good_fit)
            {
              last_bin_with_good_fit = counter;
              break;
            }
          // Count number of bins above given thresholds
          if (fit_frc[counter] > low_threshold) number_of_bins_above_low_threshold++;
          if (fit_frc[counter] > frc_significance_threshold) number_of_bins_above_significance_threshold++;
          if (fit_frc[counter] > high_threshold) number_of_bins_above_high_threshold++;
          if (counter && (fit_frc[counter] != fit_frc[counter - 1] || fit_frc[counter] == 1.0))
            last_bin_with_unique_value = counter;
        }
			//wxPrintf("%i bins out of %i checked were above significance threshold\n",number_of_bins_above_significance_threshold,number_of_bins_in_1d_spectra-first_bin_to_check);
			if ( number_of_bins_above_significance_threshold == number_of_bins_in_1d_spectra-first_bin_to_check || (last_bin_with_good_fit < 0 && number_of_bins_above_high_threshold > 3))
        last_bin_with_good_fit = std::min(last_bin_with_unique_value, number_of_bins_in_1d_spectra - 1);
			if ( number_of_bins_above_significance_threshold == 0 ) last_bin_with_good_fit = 1;
			last_bin_with_good_fit = std::min(last_bin_with_good_fit,number_of_bins_in_1d_spectra - 1);
		}
  else
		{
			last_bin_with_good_fit = 1;
		}
#ifdef DEBUG
  //MyDebugAssertTrue(last_bin_with_good_fit >= 0 && last_bin_with_good_fit < number_of_bins_in_1d_spectra,"Did not find last bin with good fit: %i", last_bin_with_good_fit);
  if (! (last_bin_with_good_fit >= 0 && last_bin_with_good_fit < number_of_bins_in_1d_spectra) )
		{
			wxPrintf("WARNING: Did not find last bin with good fit: %i\n", last_bin_with_good_fit);
		}
#else
  if (last_bin_with_good_fit < 1 && last_bin_with_good_fit >= number_of_bins_in_1d_spectra)
		{
			last_bin_with_good_fit = 1;
		}
#endif

  // Prepare output diagnostic image
  //average_spectrum->AddConstant(- average_spectrum->ReturnAverageOfRealValuesOnEdges()); // this used to be done in OverlayCTF / CTFOperation in the Fortran code
  //average_spectrum.QuickAndDirtyWriteSlice("dbg_spec_diag_3.mrc",1);
  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_before_rescaling.mrc",1);
  if (compute_extra_stats) {
    if (!RescaleSpectrumAndRotationalAverage(average_spectrum,number_of_extrema_image,ctf_values_image,number_of_bins_in_1d_spectra,spatial_frequency,rotational_average_astig,rotational_average_astig_fit,number_of_extrema_profile,ctf_values_profile,last_bin_without_aliasing,last_bin_with_good_fit))
      return false;
  }
  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_before_thresholding.mrc",1);

  /* IMOD: not needed
     average_spectrum->ComputeAverageAndSigmaOfValuesInSpectrum(	normalization_radius_min, normalization_radius_max, average,sigma);
     average_spectrum->SetMinimumAndMaximumValues(average - sigma, average + 2.0 * sigma );*/

  //average_spectrum->QuickAndDirtyWriteSlice("dbg_spec_before_overlay.mrc",1);

  // IMOD: Remove detailed results output and aliasing warning and lots of output to files

  delete comparison_object_2D;

  // IMOD: Removed terminal output on other outputs

	// Send results back
	results_array[0] = current_ctf.GetDefocus1() * pixel_size_for_fitting;				// Defocus 1 (Angstroms)
	results_array[1] = current_ctf.GetDefocus2() * pixel_size_for_fitting;				// Defocus 2 (Angstroms)
	results_array[2] = current_ctf.GetAstigmatismAzimuth() * 180.0 / PI;	// Astigmatism angle (degrees)
	results_array[3] = current_ctf.GetAdditionalPhaseShift();				// Additional phase shift (e.g. from phase plate) (radians)
	results_array[4] = - conjugate_gradient_minimizer->GetBestScore();		// CTFFIND score
  // IMOD: add || !compute_extra_stats to avoid uncomputed items
	if (last_bin_with_good_fit == 0 || !compute_extra_stats)
    {
      results_array[5] = 0.0;															//	A value of 0.0 indicates that the calculation to determine the goodness of fit failed for some reason
    }
	else
    {
      results_array[5] = pixel_size_for_fitting / spatial_frequency[last_bin_with_good_fit];		//	The resolution (Angstroms) up to which Thon rings are well fit by the CTF
    }
	if (last_bin_without_aliasing == 0 || !compute_extra_stats)
    {
      results_array[6] = 0.0;															// 	A value of 0.0 indicates that no aliasing was detected
    }
	else
    {
      results_array[6] = pixel_size_for_fitting / spatial_frequency[last_bin_without_aliasing]; 	//	The resolution (Angstroms) at which aliasing was just detected

    }

  // IMOD: Send back rotational average and other parameters
  if (rotationalAvgOut) {
    *rotationalAvgOut = (float *)malloc(number_of_bins_in_1d_spectra * sizeof(float));
    if (*rotationalAvgOut)
      for (counter = 0; counter < number_of_bins_in_1d_spectra; counter++)
        (*rotationalAvgOut)[counter] = rotational_average.data_x[counter];
  }
  if (normalizedAvgOut)
    *normalizedAvgOut = NULL;
  if (fitCurveOut)
    *fitCurveOut = NULL;
  numPointsOut = number_of_bins_in_1d_spectra;
  lastBinFreqOut = (number_of_bins_in_1d_spectra - 1.) / box_size;

  // And use the astigmatism result instead if extra-stats, plus fit and normalized
  if (compute_extra_stats && normalizedAvgOut && fitCurveOut) {
    if (*rotationalAvgOut) 
      for (counter = 0; counter < number_of_bins_in_1d_spectra; counter++)
        (*rotationalAvgOut)[counter] = rotational_average_astig[counter];
    *normalizedAvgOut = (float *)malloc(number_of_bins_in_1d_spectra * sizeof(float));
    *fitCurveOut = (float *)malloc(number_of_bins_in_1d_spectra * sizeof(float));
    if (*normalizedAvgOut)
      for (counter = 0; counter < number_of_bins_in_1d_spectra; counter++)
        (*normalizedAvgOut)[counter] = rotational_average_astig_renormalized[counter];
    if (*fitCurveOut)
      for (counter = 0; counter < number_of_bins_in_1d_spectra; counter++)
        (*fitCurveOut)[counter] = rotational_average_astig_fit[counter];
  }
  //wxPrintf("final time = %.3f\n", ctfWallTime() - wallStart);

  // IMOD: remove items here that were removed at start
	// Cleanup
	delete average_spectrum;
	delete average_spectrum_masked;
	delete current_power_spectrum;
	delete temp_image;
	delete number_of_extrema_image;
	delete ctf_values_image;
	if (compute_extra_stats)
    {
      delete [] spatial_frequency;
      delete [] rotational_average_astig;
      delete [] rotational_average_astig_renormalized;
      delete [] rotational_average_astig_fit;
      delete [] number_of_extrema_profile;
      delete [] ctf_values_profile;
      delete [] fit_frc;
      delete [] fit_frc_sigma;
    }
	delete conjugate_gradient_minimizer;



	// Return
	return true;
}

/*
 * Go from an experimental radial average with decaying Thon rings to a function between 0.0 and 1.0 for every oscillation.
 * This is done by treating each interval between a zero and an extremum of the CTF separately, and for each of them,
 * sorting and ranking the values in the radial average.
 * Each value is then replaced by its rank, modified to make it looks like a |CTF| signal.
 * This makes sense as a preparation for evaluating the quality of fit of a CTF when we want to ignore the amplitude of the Thon
 * rings and just focus on whether the fit agrees in terms of the positions of the zeros and extrema.
 * Without this, a very good fit doesn't always have a great FRC for regions where the experimental radial average is decaying rapidly.
 */
void Renormalize1DSpectrumForFRC( int number_of_bins, double average[], double fit[], float number_of_extrema_profile[])
{
	int bin_counter;
	int bin_of_previous_extremum;
	int bin_of_current_extremum;
	int i;
	int bin_of_zero;
	std::vector<float> temp_vector;
	std::vector<size_t> temp_ranks;
	//
	bin_of_previous_extremum = 0;
	bin_of_current_extremum = 0;
	for (bin_counter = 1; bin_counter < number_of_bins; bin_counter ++ )
	{
		if (number_of_extrema_profile[bin_counter]-number_of_extrema_profile[bin_counter-1] >= 0.9)
		{
			// We just passed an extremum, at bin_counter-1
			// (number_of_extrema_profile keeps track of the count of extrema before the spatial frequency corresponding to this bin)
			bin_of_current_extremum = bin_counter - 1;
			if (bin_of_previous_extremum > 0)
			{
				if ((bin_of_current_extremum - bin_of_previous_extremum >= 4 && false) || (number_of_extrema_profile[bin_counter] < 7))
				{
					// Loop from the previous extremum to the one we just found
					// (there is a zero in between, let's find it)
					// TODO: redefine the zero as the lowest point between the two extrema?
					bin_of_zero = (bin_of_current_extremum - bin_of_previous_extremum)/2 + bin_of_previous_extremum;
					for (i=bin_of_previous_extremum;i<bin_of_current_extremum;i++)
					{
						if (fit[i] < fit[i-1] && fit[i] < fit[i+1]) bin_of_zero = i;
					}
					//wxPrintf("bin zero = %i\n",bin_of_zero);

					// Now we can rank before the zero (the downslope)
					//wxPrintf("downslope (including zero)...\n");
					temp_vector.clear();
					for (i=bin_of_previous_extremum; i<=bin_of_zero;i++)
					{
						//wxPrintf("about to push back %f\n",float(average[i]));
						temp_vector.push_back(float(average[i]));
					}
					temp_ranks = rankSort(temp_vector);
					for (i=bin_of_previous_extremum; i<=bin_of_zero;i++)
					{
						//wxPrintf("replaced %f",average[i]);
						average[i] = double(float(temp_ranks.at(i-bin_of_previous_extremum))/float(temp_vector.size()-1));
						average[i] = sin(average[i] * PI * 0.5);
						//wxPrintf(" with %f\n",average[i]);
					}

					// Now we can rank after the zero (upslope)
					//wxPrintf("upslope...\n");
					temp_vector.clear();
					for (i=bin_of_zero+1; i<bin_of_current_extremum;i++)
					{
						//wxPrintf("about to push back %f\n",float(average[i]));
						temp_vector.push_back(float(average[i]));
					}
					temp_ranks = rankSort(temp_vector);
					for (i=bin_of_zero+1; i<bin_of_current_extremum;i++)
					{
						//wxPrintf("replaced %f",average[i]);
						average[i] = double(float(temp_ranks.at(i-bin_of_zero-1)+1)/float(temp_vector.size()+1));
						average[i] = sin(average[i] * PI * 0.5);
						//wxPrintf(" with %f\n",average[i]);
					}
					//MyDebugAssertTrue(abs(average[bin_of_zero]) < 0.01,"Zero bin (%i) isn't set to zero: %f\n", bin_of_zero, average[bin_of_zero]);

				}
				else
				{
					// A simpler way, without ranking, is just normalize
					// between 0.0 and 1.0 (this usually works quite well when Thon rings are on a flat background anyway)
					float min_value = 1.0;
					float max_value = 0.0;
					for (i=bin_of_previous_extremum;i<bin_of_current_extremum;i++)
					{
						if (average[i] > max_value) max_value = average[i];
						if (average[i] < min_value) min_value = average[i];
					}
					for (i=bin_of_previous_extremum;i<bin_of_current_extremum;i++)
					{
						average[i] -= min_value;
						if (max_value - min_value > 0.0001) average[i] /= (max_value - min_value);
					}
				}
			}
			bin_of_previous_extremum = bin_of_current_extremum;
		}
	}
}

//
void ComputeFRCBetween1DSpectrumAndFit( int number_of_bins, double average[], double fit[], float number_of_extrema_profile[], double frc[], double frc_sigma[], int first_fit_bin)
{
	int bin_counter;
  /* IMOD icl 11 switch to allocation
  int half_window_width[number_of_bins]; */
  int *half_window_width = new int[number_of_bins];
	int bin_of_previous_extremum;
	int i;
	int first_bin, last_bin;
	double spectrum_mean, fit_mean;
	double spectrum_sigma, fit_sigma;
	double cross_product;
	float number_of_bins_in_window;

	const int minimum_window_half_width = number_of_bins / 40;

	// First, work out the size of the window over which we'll compute the FRC value
	bin_of_previous_extremum = 0;
	for (bin_counter=1; bin_counter < number_of_bins; bin_counter++)
	{
		if (number_of_extrema_profile[bin_counter] != number_of_extrema_profile[bin_counter-1])
		{
			for (i=bin_of_previous_extremum;i<bin_counter;i++)
			{
				half_window_width[i] = std::max(minimum_window_half_width,int((1.0 + 0.1 * float(number_of_extrema_profile[bin_counter]))  * float(bin_counter - bin_of_previous_extremum + 1)));
				half_window_width[i] = std::min(half_window_width[i],number_of_bins/2 - 1);
				MyDebugAssertTrue(half_window_width[i] < number_of_bins/2,"Bad half window width: %i. Number of bins: %i\n",half_window_width[i],number_of_bins);
			}
			bin_of_previous_extremum = bin_counter;
		}
	}
	half_window_width[0] = half_window_width[1];
	for (bin_counter=bin_of_previous_extremum; bin_counter < number_of_bins; bin_counter++)
	{
		half_window_width[bin_counter] = half_window_width[bin_of_previous_extremum-1];
	}

	// Now compute the FRC for each bin
	for (bin_counter=0; bin_counter < number_of_bins; bin_counter++)
	{
		if (bin_counter < first_fit_bin)
		{
			frc[bin_counter] = 1.0;
		}
		else
		{
			spectrum_mean = 0.0;
			fit_mean = 0.0;
			spectrum_sigma = 0.0;
			fit_sigma = 0.0;
			cross_product = 0.0;
			// Work out the boundaries
			first_bin = bin_counter - half_window_width[bin_counter];
			last_bin = bin_counter + half_window_width[bin_counter];
			if (first_bin < first_fit_bin)
			{
				first_bin = first_fit_bin;
				last_bin = first_bin + 2 * half_window_width[bin_counter] + 1;
			}
			if (last_bin >= number_of_bins)
			{
				last_bin = number_of_bins - 1;
				first_bin = last_bin - 2 * half_window_width[bin_counter] - 1;
			}
			MyDebugAssertTrue(first_bin >=0 && first_bin < number_of_bins,"Bad first_bin: %i",first_bin);
			MyDebugAssertTrue(last_bin >=0 && last_bin < number_of_bins,"Bad last_bin: %i",last_bin);
			// First pass
			for (i=first_bin;i<=last_bin;i++)
			{
				spectrum_mean += average[i];
				fit_mean += fit[i];
			}
			number_of_bins_in_window = float(2 * half_window_width[bin_counter] + 1);
			//wxPrintf("bin %03i, number of extrema: %f, number of bins in window: %f\n", bin_counter, number_of_extrema_profile[bin_counter], number_of_bins_in_window);
			spectrum_mean /= number_of_bins_in_window;
			fit_mean      /= number_of_bins_in_window;
			// Second pass
			for (i=first_bin;i<=last_bin;i++)
			{
				cross_product += (average[i] - spectrum_mean) * (fit[i] - fit_mean);
				spectrum_sigma += pow(average[i] - spectrum_mean,2);
				fit_sigma += pow(fit[i] - fit_mean,2);
			}
			MyDebugAssertTrue(spectrum_sigma > 0.0 && spectrum_sigma < 10000.0,"Bad spectrum_sigma: %f\n",spectrum_sigma);
			MyDebugAssertTrue(fit_sigma > 0.0 && fit_sigma < 10000.0,"Bad fit sigma: %f\n",fit_sigma);
			if (spectrum_sigma > 0.0 && fit_sigma > 0.0)
			{
				frc[bin_counter] = cross_product / (sqrtf(spectrum_sigma/number_of_bins_in_window) * sqrtf(fit_sigma/number_of_bins_in_window)) / number_of_bins_in_window;
			}
			else
			{
				frc[bin_counter] = 0.0;
			}
			frc_sigma[bin_counter] = 2.0 / sqrtf(number_of_bins_in_window);
		}
		//wxPrintf("First fit bin: %i\n", first_fit_bin);
		MyDebugAssertTrue(frc[bin_counter] > -1.01 && frc[bin_counter] < 1.01, "Bad FRC value: %f", frc[bin_counter]);
	}
  /*IMOD icl 11 : free */
  delete [] half_window_width;
}



//
void OverlayCTF( Image *spectrum, CTF *ctf)
{
	MyDebugAssertTrue(spectrum->is_in_memory, "Spectrum memory not allocated");

	//
	EmpiricalDistribution values_in_rings;
	EmpiricalDistribution values_in_fitting_range;
	int i;
	int j;
	long address;
	float i_logi, i_logi_sq;
	float j_logi, j_logi_sq;
	float current_spatial_frequency_squared;
	float current_azimuth;
	const float lowest_freq  = pow(ctf->GetLowestFrequencyForFitting(),2);
	const float highest_freq = pow(ctf->GetHighestFrequencyForFitting(),2);
	float current_ctf_value;
	float target_sigma;

	//spectrum->QuickAndDirtyWriteSlice("dbg_spec_overlay_entry.mrc",1);

	//
	address = 0;
	for (j=0;j < spectrum->logical_y_dimension;j++)
	{
		j_logi = float(j-spectrum->physical_address_of_box_center_y) * spectrum->fourier_voxel_size_y;
		j_logi_sq = powf(j_logi,2);
		for (i=0 ;i < spectrum->logical_x_dimension; i++)
		{
			i_logi = float(i-spectrum->physical_address_of_box_center_x) * spectrum->fourier_voxel_size_x;
			i_logi_sq = powf(i_logi,2);
			//
			current_spatial_frequency_squared = j_logi_sq + i_logi_sq;
			//
			if (current_spatial_frequency_squared > lowest_freq && current_spatial_frequency_squared <= highest_freq)
			{
				current_azimuth = atan2(j_logi,i_logi);
				current_ctf_value = fabs(ctf->Evaluate(current_spatial_frequency_squared,current_azimuth));
				if (current_ctf_value > 0.5) values_in_rings.AddSampleValue(spectrum->real_values[address]);
				values_in_fitting_range.AddSampleValue(spectrum->real_values[address]);
				//if (current_azimuth <= ctf->GetAstigmatismAzimuth()  && current_azimuth >= ctf->GetAstigmatismAzimuth() - 3.1415*0.5) spectrum->real_values[address] = current_ctf_value;
				if (j < spectrum->physical_address_of_box_center_y && i < spectrum->physical_address_of_box_center_x) spectrum->real_values[address] = current_ctf_value;
			}
			if (current_spatial_frequency_squared <= lowest_freq)
			{
				spectrum->real_values[address] = 0.0;
			}
			//
			address++;
		}
		address += spectrum->padding_jump_value;
	}

	//spectrum->QuickAndDirtyWriteSlice("dbg_spec_overlay_1.mrc",1);

	/*

	// We will renormalize the experimental part of the diagnostic image
	target_sigma = sqrtf(values_in_rings.GetSampleVariance()) ;


	if (target_sigma > 0.0)
	{
		address = 0;
		for (j=0;j < spectrum->logical_y_dimension;j++)
		{
			j_logi = float(j-spectrum->physical_address_of_box_center_y) * spectrum->fourier_voxel_size_y;
			j_logi_sq = powf(j_logi,2);
			for (i=0 ;i < spectrum->logical_x_dimension; i++)
			{
				i_logi = float(i-spectrum->physical_address_of_box_center_x) * spectrum->fourier_voxel_size_x;
				i_logi_sq = powf(i_logi,2);
				//
				current_spatial_frequency_squared = j_logi_sq + i_logi_sq;
				// Normalize the experimental part of the diagnostic image
				if (i > spectrum->physical_address_of_box_center_x || j > spectrum->physical_address_of_box_center_y)
				{
					spectrum->real_values[address] /= target_sigma;
				}
				else
				{
					// Normalize the outside of the theoretical part of the diagnostic image
					if (current_spatial_frequency_squared > highest_freq) spectrum->real_values[address] /= target_sigma;
				}

				address++;
			}
			address += spectrum->padding_jump_value;
		}
	}
	*/

	//spectrum->QuickAndDirtyWriteSlice("dbg_spec_overlay_final.mrc",1);
}


// Rescale the spectrum and its 1D rotational avereage so that the peaks and troughs are at 0.0 and 1.0. The location of peaks and troughs are worked out
// by parsing the suppilied 1D average_fit array
// IMOD: switch to bool to allow error return
bool RescaleSpectrumAndRotationalAverage( Image *spectrum, Image *number_of_extrema, Image *ctf_values, int number_of_bins, double spatial_frequency[], double average[], double average_fit[], float number_of_extrema_profile[], float ctf_values_profile[], int last_bin_without_aliasing, int last_bin_with_good_fit )
{
	MyDebugAssertTrue(spectrum->is_in_memory, "Spectrum memory not allocated");
	MyDebugAssertTrue(number_of_bins > 1,"Bad number of bins: %i\n",number_of_bins);

	//
	const bool spectrum_is_blank = spectrum->IsConstant();
	const int rescale_based_on_maximum_number = 2; // This peak will be used as a renormalization.
	const int sg_width = 7;
	const int sg_order = 2;
	const bool rescale_peaks = false; // if this is false, only the background will be subtracted, the Thon rings "heights" will be unaffected
  /* IMOD icl 11 switch to allocation
	float background[number_of_bins];
	float peak[number_of_bins]; */
	float *background = new float[number_of_bins];
	float *peak = new float[number_of_bins];
	int bin_counter;
	bool at_a_maximum, at_a_minimum, maximum_at_previous_bin, minimum_at_previous_bin;
	int location_of_previous_maximum, location_of_previous_minimum;
	int current_maximum_number = 0;
	int normalisation_bin_number;
	int i;
	int j;
	bool actually_do_rescaling;
	int chosen_bin;
	long address;
	int last_bin_to_rescale;
	float min_scale_factor;
	float scale_factor;
	float rescale_peaks_to;

	Curve *minima_curve = new Curve;
	Curve *maxima_curve = new Curve;

	// Initialise arrays and variables
	for (bin_counter=0; bin_counter < number_of_bins; bin_counter++)
	{
		background[bin_counter] = 0.0;
		peak[bin_counter] = 0.0;
	}
	location_of_previous_maximum = 0;
	location_of_previous_minimum = 0;
	current_maximum_number = 0;
	at_a_maximum = false;
	at_a_minimum = true; // Note, this may not be true if we have the perfect phase plate

	//
	if ( ! spectrum_is_blank )
	{
		for (bin_counter=1; bin_counter < number_of_bins - 1; bin_counter ++)
		{
			// Remember where we were before - minimum, maximum or neither
			maximum_at_previous_bin = at_a_maximum;
			minimum_at_previous_bin = at_a_minimum;
			// Are we at a CTF min or max?
			at_a_minimum = (average_fit[bin_counter] <= average_fit[bin_counter-1]) && (average_fit[bin_counter] <= average_fit[bin_counter+1]);
			at_a_maximum = (average_fit[bin_counter] >= average_fit[bin_counter-1]) && (average_fit[bin_counter] >= average_fit[bin_counter+1]);
			// It could be that the CTF is constant in this region, in which case we stay at a minimum if we were there
			if (at_a_maximum && at_a_minimum)
			{
				at_a_minimum = minimum_at_previous_bin;
				at_a_maximum = maximum_at_previous_bin;
			}
			// Fill in values for the background or peak by linear interpolation
			if (at_a_minimum)
			{
				for (i=location_of_previous_minimum+1;i<=bin_counter;i++)
				{
					// Linear interpolation of average values at the peaks and troughs of the CTF
					background[i] = average[location_of_previous_minimum] * float(bin_counter-i) / float(bin_counter-location_of_previous_minimum) + average[bin_counter] * float(i-location_of_previous_minimum) / float(bin_counter-location_of_previous_minimum);
				}
				location_of_previous_minimum = bin_counter;
				minima_curve->AddPoint(spatial_frequency[bin_counter],average[bin_counter]);
			}
			if (at_a_maximum)
			{
				if ((! maximum_at_previous_bin) && (average_fit[bin_counter] > 0.7)) current_maximum_number = current_maximum_number + 1;
				for (i=location_of_previous_maximum+1;i<=bin_counter;i++)
				{
					// Linear interpolation of average values at the peaks and troughs of the CTF
					peak[i]       = average[location_of_previous_maximum] * float(bin_counter-i) / float(bin_counter-location_of_previous_maximum) + average[bin_counter] * float(i-location_of_previous_maximum) / float(bin_counter-location_of_previous_maximum);
					//
					if (current_maximum_number == rescale_based_on_maximum_number) normalisation_bin_number = bin_counter;
				}
				location_of_previous_maximum = bin_counter;
				maxima_curve->AddPoint(spatial_frequency[bin_counter],average[bin_counter]);
			}
			if (at_a_maximum && at_a_minimum)
			{
				MyPrintfRed("Rescale spectrum: Error. At a minimum and a maximum simultaneously.");
        /* IMOD: do not abort
				abort(); */
        delete minima_curve;
        delete maxima_curve;
        delete [] background;
        delete [] peak;
        return false;
			}
		}

		// Fit the minima and maximum curves using Savitzky-Golay smoothing
		if (maxima_curve->number_of_points > sg_width) maxima_curve->FitSavitzkyGolayToData(sg_width, sg_order);
		if (minima_curve->number_of_points > sg_width) minima_curve->FitSavitzkyGolayToData(sg_width, sg_order);

		// Replace the background and peak envelopes with the smooth min/max curves
		for (bin_counter=0;bin_counter<number_of_bins;bin_counter++)
		{
			if (minima_curve->number_of_points > sg_width) background[bin_counter] =  minima_curve->ReturnSavitzkyGolayInterpolationFromX(spatial_frequency[bin_counter]);
			if (maxima_curve->number_of_points > sg_width) peak[bin_counter]       =  maxima_curve->ReturnSavitzkyGolayInterpolationFromX(spatial_frequency[bin_counter]);
		}

		// Now that we have worked out a background and a peak envelope, let's do the actual rescaling
		actually_do_rescaling = (peak[normalisation_bin_number] - background[normalisation_bin_number]) > 0.0;
		if (last_bin_without_aliasing != 0)
		{
			last_bin_to_rescale = std::min(last_bin_with_good_fit,last_bin_without_aliasing);
		}
		else
		{
			last_bin_to_rescale = last_bin_with_good_fit;
		}
		/*  IMOD: don't need the spectrum rescaled
    if (actually_do_rescaling)
		{
			min_scale_factor = 0.2;
			rescale_peaks_to = 0.75;
			address = 0;
			for (j=0;j<spectrum->logical_y_dimension;j++)
			{
				for (i=0;i<spectrum->logical_x_dimension;i++)
				{
					chosen_bin = ReturnSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, address, ctf_values, ctf_values_profile);

          // IMOD: Return failure
          if (chosen_bin < 0)
            return false;
					if (chosen_bin <= last_bin_to_rescale)
					{
						spectrum->real_values[address] -= background[chosen_bin]; // This alone makes the spectrum look very nice already
						if (rescale_peaks) spectrum->real_values[address] /= std::min(1.0f,std::max(min_scale_factor,peak[chosen_bin]-background[chosen_bin])) / rescale_peaks_to; // This is supposed to help "boost" weak Thon rings
					}
					else
					{
						spectrum->real_values[address] -= background[last_bin_to_rescale];
						if (rescale_peaks) spectrum->real_values[address] /= std::min(1.0f,std::max(min_scale_factor,peak[last_bin_to_rescale]-background[last_bin_to_rescale])) / rescale_peaks_to;
					}
					//
					address++;
				}
				address += spectrum->padding_jump_value;
			}
		}
		else
		{
			MyDebugPrint("(RescaleSpectrumAndRotationalAverage) Warning: bad peak/background detection");
		}*/

		// Rescale the 1D average
		if (peak[normalisation_bin_number] > background[normalisation_bin_number])
		{
			for (bin_counter=0;bin_counter<number_of_bins;bin_counter++)
			{

				average[bin_counter] = (average[bin_counter] - background[bin_counter]) / (peak[normalisation_bin_number] - background[normalisation_bin_number]) * 0.95;
				// We want peaks to reach at least 0.1
				if ( ((peak[bin_counter] - background[bin_counter]) < 0.1) && (fabs(peak[bin_counter]-background[bin_counter]) > 0.000001) && bin_counter <= last_bin_without_aliasing)
				{
					average[bin_counter] = average[bin_counter] / (peak[bin_counter]-background[bin_counter]) * ( peak[normalisation_bin_number] - background[normalisation_bin_number] ) * 0.1;
				}
			}
		}
		else
		{
			MyDebugPrint("(RescaleSpectrumAndRotationalAverage): unable to rescale 1D average experimental spectrum\n");
		}


	} // end of test of spectrum_is_blank

	// Cleanup
	delete minima_curve;
	delete maxima_curve;
  /* IMOD icl 11 free and error return */
  delete [] background;
  delete [] peak;
  return true;
}

// IMOD: return error
bool ComputeRotationalAverageOfPowerSpectrum( Image *spectrum, CTF *ctf, Image *number_of_extrema, Image *ctf_values, int number_of_bins, double spatial_frequency[], double average[], double average_fit[], double average_rank[], float number_of_extrema_profile[], float ctf_values_profile[])
{
	MyDebugAssertTrue(spectrum->is_in_memory, "Spectrum memory not allocated");
	MyDebugAssertTrue(number_of_extrema->is_in_memory,"Number of extrema image not allocated");
	MyDebugAssertTrue(ctf_values->is_in_memory,"CTF values image not allocated");
	MyDebugAssertTrue(spectrum->HasSameDimensionsAs(number_of_extrema),"Spectrum and number of extrema images do not have same dimensions");
	MyDebugAssertTrue(spectrum->HasSameDimensionsAs(ctf_values),"Spectrum and CTF values images do not have same dimensions");
	//
	const bool spectrum_is_blank = spectrum->IsConstant();
	const float min_angular_distances_from_axes_radians = 10.0 / 180.0 * PI;
	int counter;
	float azimuth_of_mid_defocus;
	float angular_distance_from_axes;
	float current_spatial_frequency_squared;
  /* IMOD icl 11 switch to allocation
  int number_of_values[number_of_bins];*/
	int *number_of_values = new int[number_of_bins];
	int i, j;
	long address;
	float ctf_diff_from_current_bin;
	int chosen_bin;

	// Initialise the output arrays
	for (counter=0; counter<number_of_bins; counter++)
	{
		average[counter] = 0.0;
		average_fit[counter] = 0.0;
		average_rank[counter] = 0.0;
		ctf_values_profile[counter] = 0.0;
		number_of_values[counter] = 0;
	}

	//
	if (! spectrum_is_blank)
	{
		// For each bin of our 1D profile we compute the CTF. We choose the azimuth to be mid way between the two defoci of the astigmatic CTF
		azimuth_of_mid_defocus = ctf->GetAstigmatismAzimuth() + PI * 0.25;
		// We don't want the azimuth too close to the axes, which may have been blanked by the central-cross-artefact-suppression-system (tm)
    /* IMOD icl 11 add (float) */
		angular_distance_from_axes = fmod(azimuth_of_mid_defocus, (float)PI * 0.5f);
		if(fabs(angular_distance_from_axes) < min_angular_distances_from_axes_radians)
		{
			if (angular_distance_from_axes > 0.0)
			{
				azimuth_of_mid_defocus = min_angular_distances_from_axes_radians;
			}
			else
			{
				azimuth_of_mid_defocus = - min_angular_distances_from_axes_radians;
			}
		}
		if (fabs(angular_distance_from_axes) > 0.5 * PI - min_angular_distances_from_axes_radians)
		{
			if (angular_distance_from_axes > 0.0)
			{
				azimuth_of_mid_defocus = PI * 0.5 - min_angular_distances_from_axes_radians;
			}
			else
			{
				azimuth_of_mid_defocus = - PI * 0.5 + min_angular_distances_from_axes_radians;
			}
		}
		// Now that we've chosen an azimuth, we can compute the CTF for each bin of our 1D profile
		for (counter=0;counter < number_of_bins; counter++)
		{
			current_spatial_frequency_squared = powf(float(counter) * spectrum->fourier_voxel_size_y, 2);
			spatial_frequency[counter] = sqrt(current_spatial_frequency_squared);
			ctf_values_profile[counter] = ctf->Evaluate(current_spatial_frequency_squared,azimuth_of_mid_defocus);
			number_of_extrema_profile[counter] = ctf->ReturnNumberOfExtremaBeforeSquaredSpatialFrequency(current_spatial_frequency_squared,azimuth_of_mid_defocus);
		}

		// Now we can loop over the spectrum again and decide to which bin to add each component
		address = 0;
		for (j=0; j<spectrum->logical_y_dimension; j++)
		{
			for (i=0; i < spectrum->logical_x_dimension; i++)
			{
				ctf_diff_from_current_bin = std::numeric_limits<float>::max();
				chosen_bin = ReturnSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, address, ctf_values, ctf_values_profile);
        // IMOD: return failure
        if (chosen_bin < 0) {
          wxPrintf("ReturnSpectrumBinNumber failed to find bin for i = %d  j = %d\n", i, j);
          return false;
        }
				average[chosen_bin] += spectrum->real_values[address];
				number_of_values[chosen_bin]++;
				//
				address++;
			}
			address += spectrum->padding_jump_value;
		}

		// Do the actual averaging
		for (counter = 0; counter < number_of_bins; counter++)
		{
			if (number_of_values[counter] > 0)
			{
				average[counter] = average[counter] / float(number_of_values[counter]);
			}
			else
			{
				average[counter] = 0.0;
			}
			average_fit[counter] = fabs(ctf_values_profile[counter]);
		}

	}

	// Compute the rank version of the rotational average
	for (counter = 0; counter < number_of_bins; counter ++ ) { average_rank[counter] = average[counter]; }
	Renormalize1DSpectrumForFRC(number_of_bins,average_rank,average_fit,number_of_extrema_profile);
  // IMOD icl 11: delete and error return 
  delete [] number_of_values;
  return true;
}


int ReturnSpectrumBinNumber(int number_of_bins, float number_of_extrema_profile[], Image *number_of_extrema, long address, Image *ctf_values, float ctf_values_profile[])
{
	int current_bin;
	float diff_number_of_extrema;
	float diff_number_of_extrema_previous;
	float diff_number_of_extrema_next;
	float ctf_diff_from_current_bin;
	float ctf_diff_from_current_bin_old;
	int chosen_bin;
	//
	//MyDebugPrint("address: %li - number of extrema: %f - ctf_value: %f\n", address, number_of_extrema->real_values[address], ctf_values->real_values[address]);
	MyDebugAssertTrue(address < number_of_extrema->real_memory_allocated,"Oops, bad address: %li\n",address);
	// Let's find the bin which has the same number of preceding extrema and the most similar ctf value
	ctf_diff_from_current_bin = std::numeric_limits<float>::max();
	chosen_bin = -1;
	for (current_bin=0; current_bin < number_of_bins; current_bin++)
	{
		diff_number_of_extrema = fabs(number_of_extrema->real_values[address] - number_of_extrema_profile[current_bin]);
		if (current_bin > 0)
		{
			diff_number_of_extrema_previous = fabs(number_of_extrema->real_values[address]- number_of_extrema_profile[current_bin-1]);
		}
		else
		{
			diff_number_of_extrema_previous = std::numeric_limits<float>::max();
		}
		if (current_bin < number_of_bins - 1)
		{
			diff_number_of_extrema_next = fabs(number_of_extrema->real_values[address] - number_of_extrema_profile[current_bin+1]);
		}
		else
		{
			diff_number_of_extrema_next = std::numeric_limits<float>::max();
		}
		//
		if (number_of_extrema->real_values[address] > number_of_extrema_profile[number_of_bins-1])
		{
			chosen_bin = number_of_bins - 1;
		}
		else
		{
			if ( diff_number_of_extrema <= 0.01 || (  diff_number_of_extrema <  diff_number_of_extrema_previous &&
					                                  diff_number_of_extrema <= diff_number_of_extrema_next &&
													  number_of_extrema_profile[std::max(current_bin-1,0)] != number_of_extrema_profile[std::min(current_bin+1,number_of_bins-1)]  ) )
			{
				// We're nearly there
				// Let's look for the position for the nearest CTF value
				ctf_diff_from_current_bin_old = ctf_diff_from_current_bin;
				ctf_diff_from_current_bin = fabs(ctf_values->real_values[address] - ctf_values_profile[current_bin]);
				if (ctf_diff_from_current_bin < ctf_diff_from_current_bin_old)
				{
					//MyDebugPrint("new chosen bin: %i\n",current_bin);
					chosen_bin = current_bin;
				}
			}
		}
	}
	if (chosen_bin == -1)
	{
		MyPrintfRed("Could not find bin\n");
    /* IMOD: do not abort
		abort(); */
    return -1;
	}
	else
	{
		//MyDebugAssertTrue(chosen_bin > 0 && chosen_bin < number_of_bins,"Oops, bad chosen bin number: %i (number of bins = %i)\n",chosen_bin,number_of_bins);
		//MyDebugPrint("final chosen bin = %i\n", chosen_bin);
		return chosen_bin;
	}
}
/*
integer function ComputePowerSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, &
                                                i,j,ctf_values,ctf_values_profile) result(chosen_bin)
    integer,        intent(in)  ::  number_of_bins
    real,           intent(in)  ::  number_of_extrema_profile(:)
    type(Image),    intent(in)  ::  number_of_extrema
    integer,        intent(in)  ::  i,j                         !<  Physical memory address
    type(Image),    intent(in)  ::  ctf_values
    real,           intent(in)  ::  ctf_values_profile(:)
    ! private variables
    integer     ::  current_bin
    real        ::  diff_number_of_extrema, diff_number_of_extrema_previous, diff_number_of_extrema_next
    real        ::  ctf_diff_from_current_bin
    real        ::  ctf_diff_from_current_bin_old
    ! Let's find the bin which has the same number of preceding extrema and the most similar ctf value
    ctf_diff_from_current_bin = huge(1.0e0)
    chosen_bin = 0
    do current_bin=1,number_of_bins
        diff_number_of_extrema  = abs(number_of_extrema%real_values(i,j,1) - number_of_extrema_profile(current_bin))
        if (current_bin .gt. 1) then
            diff_number_of_extrema_previous = abs(number_of_extrema%real_values(i,j,1) &
                                                - number_of_extrema_profile(current_bin-1))
        else
            diff_number_of_extrema_previous = huge(1.0e0)
        endif
        if (current_bin .lt. number_of_bins) then
            diff_number_of_extrema_next     = abs(number_of_extrema%real_values(i,j,1) &
                                                - number_of_extrema_profile(current_bin+1))
        else
            diff_number_of_extrema_next = huge(1.0e0)
        endif
        if (number_of_extrema%real_values(i,j,1) .gt. number_of_extrema_profile(number_of_bins)) then
            chosen_bin = number_of_bins
        else
            if (        diff_number_of_extrema .le. 0.01 &
                .or.    (     diff_number_of_extrema .lt. diff_number_of_extrema_previous &
                        .and. diff_number_of_extrema .le. diff_number_of_extrema_next &
                        .and. number_of_extrema_profile(max(current_bin-1,1)) &
                            .ne. number_of_extrema_profile(min(current_bin+1,number_of_bins))) &
                ) then
                ! We're nearly there
                ! Let's look for the position of the nearest CTF value
                ctf_diff_from_current_bin_old = ctf_diff_from_current_bin
                ctf_diff_from_current_bin = abs(ctf_values%real_values(i,j,1) - ctf_values_profile(current_bin))
                if (ctf_diff_from_current_bin .lt. ctf_diff_from_current_bin_old) then
                    chosen_bin = current_bin
                endif
            endif
        endif
    enddo
    if (chosen_bin .eq. 0) then
        print *, number_of_extrema_profile
        print *, i, j, number_of_extrema%real_values(i,j,1), ctf_values%real_values(i,j,1)
        print *, diff_number_of_extrema, diff_number_of_extrema_previous, diff_number_of_extrema_next
        call this_program%TerminateWithFatalError('ComputeRotationalAverageOfPowerSpectrum','Could not find bin')
    endif
end function ComputePowerSpectrumBinNumber
*/


// Compute an image where each pixel stores the number of preceding CTF extrema. This is described as image "E" in Rohou & Grigorieff 2015 (see Fig 3)
void ComputeImagesWithNumberOfExtremaAndCTFValues(CTF *ctf, Image *number_of_extrema, Image *ctf_values)
{
	MyDebugAssertTrue(number_of_extrema->is_in_memory,"Memory not allocated");
	MyDebugAssertTrue(ctf_values->is_in_memory,"Memory not allocated");
	MyDebugAssertTrue(ctf_values->HasSameDimensionsAs(number_of_extrema),"Images do not have same dimensions");

	int i, j;
	float i_logi, i_logi_sq;
	float j_logi, j_logi_sq;
	float current_spatial_frequency_squared;
	float current_azimuth;
	long address;

	address = 0;
	for (j=0;j<number_of_extrema->logical_y_dimension;j++)
	{
		j_logi = float(j - number_of_extrema->physical_address_of_box_center_y) * number_of_extrema->fourier_voxel_size_y;
		j_logi_sq = pow(j_logi,2);
		for (i=0;i<number_of_extrema->logical_x_dimension;i++)
		{
			i_logi = float(i - number_of_extrema->physical_address_of_box_center_x) * number_of_extrema->fourier_voxel_size_x;
			i_logi_sq = pow(i_logi,2);
			// Where are we?
			current_spatial_frequency_squared = j_logi_sq + i_logi_sq;
			if (current_spatial_frequency_squared > 0.0)
			{
				current_azimuth = atan2(j_logi,i_logi);
			}
			else
			{
				current_azimuth = 0.0;
			}
			//
			ctf_values->real_values[address] = ctf->Evaluate(current_spatial_frequency_squared,current_azimuth);
			number_of_extrema->real_values[address] = ctf->ReturnNumberOfExtremaBeforeSquaredSpatialFrequency(current_spatial_frequency_squared,current_azimuth);
			//
			address++;
		}
		address += number_of_extrema->padding_jump_value;
	}

	number_of_extrema->is_in_real_space = true;
	ctf_values->is_in_real_space = true;
}

// Align rotationally a (stack) of image(s) against another image. Return the rotation angle that gives the best normalised cross-correlation.
float FindRotationalAlignmentBetweenTwoStacksOfImages(Image *self, Image *other_image, int number_of_images, float search_half_range, float search_step_size, float minimum_radius, float maximum_radius)
{
	MyDebugAssertTrue(self[0].is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(self[0].is_in_real_space, "Not in real space");
	MyDebugAssertTrue(self[0].logical_z_dimension == 1, "Meant for images, not volumes");
	MyDebugAssertTrue(other_image[0].is_in_memory, "Memory not allocated - other_image");
	MyDebugAssertTrue(other_image[0].is_in_real_space, "Not in real space - other_image");
	MyDebugAssertTrue(other_image[0].logical_z_dimension == 1, "Meant for images, not volumes - other_image");
	MyDebugAssertTrue(self[0].HasSameDimensionsAs(&other_image[0]),"Images and reference images do not have same dimensions.");

	// Local variables
	const float minimum_radius_sq = pow(minimum_radius,2);
	const float maximum_radius_sq = pow(maximum_radius,2);
	const float inverse_logical_x_dimension = 1.0 / float(self[0].logical_x_dimension);
	const float inverse_logical_y_dimension = 1.0 / float(self[0].logical_y_dimension);
	float best_cc = - std::numeric_limits<float>::max();
	float best_rotation = - std::numeric_limits<float>::max();
	float current_rotation = - search_half_range;
	float current_rotation_rad;
	EmpiricalDistribution cc_numerator_dist;
	EmpiricalDistribution cc_denom_self_dist;
	EmpiricalDistribution cc_denom_other_dist;
	int current_image;
	int i, i_logi;
	float i_logi_frac, ii_phys;
	int j, j_logi;
	float j_logi_frac, jj_phys;
	float current_interpolated_value;
	long address_in_other_image;
	float current_cc;



	// Loop over possible rotations
	while ( current_rotation < search_half_range + search_step_size )
	{

		current_rotation_rad = current_rotation / 180.0 * PI;
		cc_numerator_dist.Reset();
		cc_denom_self_dist.Reset();
		cc_denom_other_dist.Reset();
		// Loop over the array of images
		for (current_image=0; current_image < number_of_images; current_image++)
		{
			// Loop over the other (reference) image
			address_in_other_image = 0;
			for (j=0; j < other_image[0].logical_y_dimension; j++)
			{
				j_logi = j - other_image[0].physical_address_of_box_center_y;
				j_logi_frac = pow(j_logi * inverse_logical_y_dimension,2);
				for (i=0; i < other_image[0].logical_x_dimension; i++)
				{
					i_logi = i - other_image[0].physical_address_of_box_center_x;
					i_logi_frac = pow(i_logi * inverse_logical_x_dimension,2) + j_logi_frac;

					if (i_logi_frac >= minimum_radius_sq && i_logi_frac <= maximum_radius_sq)
					{
						// We do ccw rotation to go from other_image (reference) to self (input image)
						ii_phys = i_logi * cos(current_rotation_rad) - j_logi * sin(current_rotation_rad) + self[0].physical_address_of_box_center_x ;
						jj_phys = i_logi * sin(current_rotation_rad) + j_logi * cos(current_rotation_rad) + self[0].physical_address_of_box_center_y ;
						//
						if (int(ii_phys) > 0 && int(ii_phys)+1 < self[0].logical_x_dimension && int(jj_phys) > 0 && int(jj_phys)+1 < self[0].logical_y_dimension ) // potential optimization: we have to compute the floor and ceiling in the interpolation routine. Is it not worth doing the bounds checking in the interpolation routine somehow?
						{
							self[0].GetRealValueByLinearInterpolationNoBoundsCheckImage(ii_phys,jj_phys,current_interpolated_value);
							//MyDebugPrint("%g %g\n",current_interpolated_value,other_image[0].real_values[address_in_other_image]);
							cc_numerator_dist.AddSampleValue(current_interpolated_value * other_image[current_image].real_values[address_in_other_image]);
							cc_denom_other_dist.AddSampleValue(pow(other_image[0].real_values[address_in_other_image],2)); // potential optimization: since other_image is not being rotated, we should only need to compute this quantity once, not for every potential rotation
							cc_denom_self_dist.AddSampleValue(pow(current_interpolated_value,2));
						}
					}
					address_in_other_image++;
				} // i
				address_in_other_image += other_image[0].padding_jump_value;
			} // end of loop over other (reference) image
		} // end of loop over array of images

		current_cc = cc_numerator_dist.GetSampleSum() / sqrt(cc_denom_other_dist.GetSampleSum()*cc_denom_self_dist.GetSampleSum());


		if (current_cc > best_cc)
		{
			best_cc = current_cc;
			best_rotation = current_rotation;
		}

		// Increment the rotation
		current_rotation += search_step_size;

	} // end of loop over rotations

	return best_rotation;
}


