// When compiling in SerialEM, the default is to export functions, this is a define
// used in SerialEM which will undefine that export, then dllexport.h will import
#if defined(NO_WARN_MBCS_MFC_DEPRECATION) && defined(DLL_IM_EX)
#undef DLL_IM_EX
#endif
#include "dllexport.h"

struct CtffindParams {
	float 				pixel_size_of_input_image;
	float 		acceleration_voltage;
	float       	spherical_aberration;
	float 		amplitude_contrast;
	int         	box_size;
	float 		minimum_resolution;
	float       	maximum_resolution;
	float       	minimum_defocus;
	float       	maximum_defocus;
	float       	defocus_search_step;
	bool			slower_search;
	float       	astigmatism_tolerance;
	bool       	find_additional_phase_shift;
	float  		minimum_additional_phase_shift;
	float			maximum_additional_phase_shift;
	float			additional_phase_shift_search_step;
	bool  		astigmatism_is_known;
	float 		known_astigmatism;
	float 		known_astigmatism_angle;
	bool 			compute_extra_stats;
};

typedef void (*CharArgType)(const char *);
typedef int (*WriteSliceType)(const char *, float *, int , int);

DLL_IM_EX bool ctffind(CtffindParams &params, float *spectrumArray, int nxDimIn,
                       float *results_array, float **rotationalAvgOut,
                       float **normalizedAvgOut, float **fitCurveOut,
                       int &numPointsOut, float &lastBinFreqOut);

DLL_IM_EX void ctffindSetSliceWriteFunc(WriteSliceType func);
DLL_IM_EX void ctffindSetPrintFunc(CharArgType func);
