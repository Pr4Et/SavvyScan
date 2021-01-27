typedef struct Peak {
  float x;
  float y;
  float z;
  float value;
  long  physical_address_within_image;
} Peak;

typedef struct Kernel2D {
  int   pixel_index[4];
  float pixel_weight[4];
} Kernel2D;

typedef struct CurvePoint {
  int   index_m;
  int   index_n;
  float value_m;
  float value_n;
} CurvePoint;

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <cstdarg>
#include <cfloat>
#include <complex>
const std::complex<float> I(0.0,1.0);
#include <iterator>
#include <utility>
#include <vector>
#include <limits>
#include <math.h>
#include "defines.h"
#include "functions.h"
#include "randomnumbergenerator.h"
#include "matrix.h"
#include "ctf.h"
#include "curve.h"
#include "angles_and_shifts.h"
#include "empirical_distribution.h"
#include "image.h"
#include "brute_force_search.h"
#include "conjugate_gradient.h"

extern RandomNumberGenerator global_random_number_generator;

