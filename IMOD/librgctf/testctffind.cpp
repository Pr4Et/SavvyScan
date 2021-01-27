#include "b3dutil.h"
#include "mrcslice.h"
#include "ctffind.h"
#include "parse_params.h"
#include "cfft.h"

static int writeSlice(const char *filename, float *data, int xsize, int ysize)
{
  Islice slice;
  sliceInit(&slice, xsize, ysize, MRC_MODE_FLOAT, data);
	return sliceWriteMRCfile(filename, &slice);
}


int main( int argc, char *argv[] )
{
  CtffindParams params;
  const char *options[] = 
    {"volt::F:Voltage", "sph::F:Sperical Abb", "box::I:box size", "rmin::F:Minimum res",
     "rmax::F:Maximum res", "dmin::F:Minimum defocus", "dmax::F:Maximum defocus",
     "dstep::F:Defocus step", "fast::B:Faster search", "atol::F:Astigmatism tolerance",
     "known::B:Astigmatism is known", "astig::F:Known astig", "angle::F:Known angle",
     "noex::B:No extra computations", "phase::B:Find phase", "phmin::F:Minimum phase",
     "phmax::F:Maximum phase", "phstep::F:Phase step", 
     "resamp::F:Resample amplitude spectrum to this resolution"};
  int numOptions = sizeof(options) / sizeof(char *);
  int numOptArg, numNonOptArg, val, numPoints, padSize, start, end, useBox;
  char *filename;
  MrcHeader hdata;
  void *imageArray;
  float *spectrum, *rotationalAvg, *normalizedAvg, *fitCurve;
  float results_array[7], lastBinFreq, xscale, yscale, zscale, resampRes = 2.8;
  double wallStart  = wallTime();

  ctffindSetSliceWriteFunc(writeSlice);

  setExitPrefix("testctffind");
  PipExitOnError(0, "testctffind");
  PipParseInput(argc, argv, options, numOptions, &numOptArg, &numNonOptArg);
  if (!numNonOptArg) {
    PipPrintHelp("testctffind", 0, 1, 0);
    exit(0);
  }
  params.acceleration_voltage = 200.;
  params.spherical_aberration = 2.;
  params.amplitude_contrast = 0.07;
  params.box_size = 256;
  params.minimum_resolution = 50.;
  params.maximum_resolution = 10;
  params.minimum_defocus = 5000.;
  params.maximum_defocus = 80000.;
  params.defocus_search_step = 500.;
  params.astigmatism_tolerance = -100.;
  params.additional_phase_shift_search_step = .1;
  params.known_astigmatism = 0.;
  params.known_astigmatism_angle = 0.;
  params.minimum_additional_phase_shift = 0.;
  params.maximum_additional_phase_shift = 0.;
  
  PipGetNonOptionArg(0, &filename);
  PipGetFloat("volt", &params.acceleration_voltage);
  PipGetFloat("sph", &params.spherical_aberration);
  PipGetInteger("box", &params.box_size);
  PipGetFloat("rmin", &params.minimum_resolution);
  PipGetFloat("rmax", &params.maximum_resolution);
  PipGetFloat("dmin", &params.minimum_defocus);
  PipGetFloat("dmax", &params.maximum_defocus);
  PipGetFloat("dstep", &params.defocus_search_step);
  PipGetFloat("atol", &params.astigmatism_tolerance);
  val = 0;
  PipGetBoolean("fast", &val);
  params.slower_search = val == 0;
  val = 0;
  PipGetBoolean("known", &val);
  params.astigmatism_is_known = val != 0;
  val = 0;
  PipGetBoolean("phase", &val);
  params.find_additional_phase_shift = val != 0;
  val = 0;
  PipGetBoolean("noex", &val);
  params.compute_extra_stats = val == 0;
  if (params.astigmatism_is_known)
    PipGetFloat("astig", &params.known_astigmatism);
  PipGetFloat("angle", &params.known_astigmatism_angle);
  if (params.find_additional_phase_shift) {
    PipGetFloat("phmin", &params.minimum_additional_phase_shift);
    PipGetFloat("phmax", &params.maximum_additional_phase_shift);
    PipGetFloat("phstep", &params.additional_phase_shift_search_step);
  }
  params.box_size = 2 * ((params.box_size + 1) / 2);
  PipGetFloat("resamp", &resampRes);

  FILE *fp = fopen(filename, "rb");
  if (!fp)
    exitError("Opening file %s", filename);
  if (mrc_head_read(fp, &hdata))
    exitError("Reading header of file %s - %s", filename, b3dGetError());
  free(filename);

  mrc_get_scale(&hdata, &xscale, &yscale, &zscale);
  params.pixel_size_of_input_image = xscale;

  useBox = params.box_size;
  if (resampRes > params.pixel_size_of_input_image * 2.)
    useBox = 2 * (B3DNINT(1. + 0.5 * params.box_size * resampRes / 
                          params.pixel_size_of_input_image) / 2);

  spectrum = B3DMALLOC(float, useBox * (useBox + 2));
  if (!spectrum)
    exitError("Allocating spectrum array");

  for (int iz = 0; iz < hdata.nz; iz++) {
    imageArray = mrc_mread_slice(fp, &hdata, iz, 'Z');
    if (!imageArray)
      exitError("Reading section 0 - %s", b3dGetError());

    padSize = B3DMAX(hdata.nx, hdata.ny);
    padSize = niceFrame(padSize, 2, niceFFTlimit());
    val = spectrumScaled(imageArray, hdata.mode, hdata.nx, hdata.ny, spectrum, -padSize, 
                         useBox, 0, 0., -1, todfft);
    if (val)
      exitError("Making reduced spectrum, error code %d", val);
    if (useBox > params.box_size) {
      start = (useBox - params.box_size) / 2;
      end = start + params.box_size - 1;
      if (extractAndBinIntoArray(spectrum, MRC_MODE_FLOAT, useBox + 2, start, end, start,
                                 end, 1, spectrum, params.box_size + 2, 0, 0, 0, &val,
                                 &val))
        exitError("Extracting reduced spectrum from larger box to smaller");
      params.pixel_size_of_input_image *= (float)useBox / (float)params.box_size;
    }

    if (!ctffind(params, spectrum, params.box_size + 2, results_array, &rotationalAvg,
                 &normalizedAvg, &fitCurve, numPoints, lastBinFreq))
      exitError("An error occurred in ctffind");

    if (hdata.nz > 1) {
      printf("%1.f  %.1f  %.2f  %.4f  %.4f  %.1f  %.1f\n", results_array[0], 
             results_array[1], results_array[2], results_array[3], results_array[4],
             results_array[5], results_array[6]);
    } else {
      printf("Defocus 1 %.2f  2 %.2f  angle %.2f", results_array[0], results_array[1],
             results_array[2]);
      if (params.find_additional_phase_shift)
        printf("    Added phase shift %.4f", results_array[3]);
      printf("\nScore %.4f\n", results_array[4]);
      if (params.compute_extra_stats) {
        printf("Thon rings well fit to %.1f\n", results_array[5]);
        if (results_array[6])
          printf("CTF aliasing detected at %.1f\n", results_array[6]);
      }
    }
    free(imageArray);
    B3DFREE(rotationalAvg);
    B3DFREE(normalizedAvg);
    B3DFREE(fitCurve);
  }
  printf("ctffind time %.3f\n", wallTime() - wallStart);
  free(spectrum);
  return 0;
}

