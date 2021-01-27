/**

Livewire - Core code for running the Livewire algorithm
Copyright (C) 2011  Jeffrey Bush  jeff@coderforlife.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

**/

#include "Weights.h"

#include "Colors.h"
#include "FilterBase.h"

#include <QRgb>

#include <deque>

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <assert.h>

const double PI_4 = M_PI / 4;

#ifndef IMOD_PLUGIN
//#define SAVE_WEIGHT_IMAGE
#endif

#ifndef IMOD_PLUGIN
#include "BitmapWriter.h"
#endif

#define MIN(a, b)	(((a) < (b)) ? (a) : (b))
#define MAX(a, b)	(((a) > (b)) ? (a) : (b))

using namespace Livewire;

/// <summary>
/// Calculates a sigmoid function for a single value with the given parameters.
/// Equation:
/// f(x) = floor(max_y / (1 + exp(slope * (halfmax_x - x)));
/// </summary>
/// <param name="x">The value of x to calculate the sigmoid for</param>
/// <param name="slope_inv">The inverse of the slope of the function, the slope follows that smaller values create a more linear relationship while larger values tend toward a step function</param>
/// <param name="halfmax_x">The value of x at which the result is half of the max, below this values will become smaller and above they will become larger</param>
/// <param name="max_y">The maximum value of the result</param>
/// <returns>The evaluated sigmoid function value</returns>
template<int min_y, int max_y>
inline static int sigmoid(const int slope_inv, const int halfmax_x, const int x) { return (int)(min_y + (max_y - min_y) / (1 + exp((halfmax_x - x) / (double)slope_inv))); }

#pragma region Coalescing Functions
static void CoalesceGrayByte(const uint w, const uint h, const uint stride, const byte *in, byte *out)
{
	for (const byte *out_end = out + w*h; out < out_end; out += w, in += stride)
		memcpy(out, in, w);
}

static void CoalesceGrayUShort(const uint w, const uint h, const uint stride, const byte *in, byte *out)
{
	for (const byte *out_end = out + w*h; out < out_end; out += w, in += stride)
	{
		const unsigned short *_in = (const unsigned short*)in;
		for (uint x = 0; x < w; ++x)
			out[x] = (_in[x] >> 8) & 0xFF;
	}
}

template <uint channel>
static void CoalesceChannel(const uint w, const uint h, const uint stride, const byte *in, byte *out)
{
	static const uint shift = channel * 8;
	for (const byte *out_end = out + w*h; out < out_end; out += w, in += stride)
	{
		const QRgb *_in = (const QRgb*)in;
		for (uint x = 0; x < w; ++x)
			out[x] = (_in[x] >> shift) & 0xFF;
	}
}
#define CoalesceRed(w, h, stride, in, out)		CoalesceChannel<2>(w, h, stride, in, out)
#define CoalesceGreen(w, h, stride, in, out)	CoalesceChannel<1>(w, h, stride, in, out)
#define CoalesceBlue(w, h, stride, in, out)		CoalesceChannel<0>(w, h, stride, in, out)

template <uint r, uint g, uint b>
static void CoalesceWeightedAvgRGB(const uint w, const uint h, const uint stride, const byte *in, byte *out)
{
	for (const byte *out_end = out + w*h; out < out_end; out += w, in += stride)
	{
		const QRgb *_in = (const QRgb*)in;
		for (uint x = 0; x < w; ++x)
		{
			QRgb rgb = _in[x];
			out[x] = (r * qRed(rgb) + g * qGreen(rgb) + b * qBlue(rgb)) / (r + g + b);
		}
	}
}

#define CoalesceAvgRGB(w, h, stride, in, out)		CoalesceWeightedAvgRGB<   1,   1,   1>(w, h, stride, in, out)
#define CoalesceLuma(w, h, stride, in, out)			CoalesceWeightedAvgRGB<2126,7152, 722>(w, h, stride, in, out)
#define CoalesceLuma601(w, h, stride, in, out)		CoalesceWeightedAvgRGB< 299, 587, 114>(w, h, stride, in, out)
#define CoalesceLumaSMPTE(w, h, stride, in, out)	CoalesceWeightedAvgRGB< 212, 701,  87>(w, h, stride, in, out)

template <void conv(byte R, byte G, byte B, byte& H, byte& S, byte& X)>
static void CoalesceWeightedHSX(const uint w, const uint h, const uint stride, const byte *in, byte *out)
{
	for (const byte *out_end = out + w*h; out < out_end; out += w, in += stride)
	{
		const QRgb *_in = (const QRgb*)in;
		for (uint x = 0; x < w; ++x)
		{
			QRgb rgb = _in[x];
			byte H, S, X;
			conv(qRed(rgb), qGreen(rgb), qBlue(rgb), H, S, X);
			out[x] = 0.6 * X + 0.3 * H + 0.1 * S;
		}
	}
}

#define CoalesceWeightedHSV(w, h, stride, in, out)	CoalesceWeightedHSX<RGBtoHSV>(w, h, stride, in, out)
#define CoalesceWeightedHSL(w, h, stride, in, out)	CoalesceWeightedHSX<RGBtoHSL>(w, h, stride, in, out)
#define CoalesceWeightedHSI(w, h, stride, in, out)	CoalesceWeightedHSX<RGBtoHSI>(w, h, stride, in, out)

#pragma endregion

#pragma region Binning and Noise Reduction Functions
template<uint windowSize>
struct MedianFilter : public RCRSFilter<windowSize>
{
	virtual byte SelectValue(const byte* list, const byte)
	{
		static const uint ws2_2 = windowSize*windowSize / 2;
		return (windowSize & 1) ? list[ws2_2] : (byte)(((uint)list[ws2_2 - 1] + list[ws2_2]) / 2);
	}
	virtual byte SelectValue(const byte* list, uint count, const byte)
	{
		const uint half = count / 2;
		return (count & 1) ? list[half] : (byte)(((uint)list[half - 1] + list[half]) / 2);
	}
};
template<uint windowSize>
struct MedianBinner : public WindowBinner<MedianFilter<windowSize> > { };

template<uint windowSize> struct MeanKernel { static const uint Matrix[windowSize][windowSize], Vert[windowSize], Horz[windowSize], Total = windowSize*windowSize; };
#define MEAN_KERNEL(N, ...) template<> const uint MeanKernel<N>::Vert[N] = {__VA_ARGS__}; template<> const uint MeanKernel<N>::Horz[N] = {__VA_ARGS__}; template<> const uint MeanKernel<N>::Matrix[N][N]
MEAN_KERNEL(2, 1, 1)          = { {1,1}, {1,1} };
MEAN_KERNEL(3, 1, 1, 1)       = { {1,1,1}, {1,1,1}, {1,1,1} };
MEAN_KERNEL(4, 1, 1, 1, 1)    = { {1,1,1,1}, {1,1,1,1}, {1,1,1,1}, {1,1,1,1} };
MEAN_KERNEL(5, 1, 1, 1, 1, 1) = { {1,1,1,1,1}, {1,1,1,1,1}, {1,1,1,1,1}, {1,1,1,1,1}, {1,1,1,1,1} };
template<uint windowSize>
struct MeanFilter : public SepConvFilter<windowSize, MeanKernel<windowSize> >
{
	// This could be done with built-in WindowFilter for convolutions, but since all values are one we optimize a little bit
	// TODO: this needs to also be updated when SepConvFilter finally gets implemented
	virtual byte FilterWindow(const byte** window)
	{
		static const uint ws2 = windowSize * windowSize;
		uint v = 0;
		for (uint y = 0; y < windowSize; ++y)
			for (uint x = 0; x < windowSize; ++x)
				v += window[y][x];
		return (byte)(v / ws2);
	}
	virtual byte FilterWindow(const byte** window, uint w, uint h, uint, uint)
	{
		uint v = 0;
		for (uint y = 0; y < h; ++y)
			for (uint x = 0; x < w; ++x)
				v += window[y][x];
		return (byte)(v / (w * h));
	}
};
template<uint windowSize>
struct MeanBinner : public WindowBinner<MeanFilter<windowSize> > { };

template<uint windowSize> struct GaussianKernel { static const uint Matrix[windowSize][windowSize], Vert[windowSize], Horz[windowSize], Total; };
#define BRACED_INIT_LIST(...) {__VA_ARGS__}
#define GAUSSIAN_KERNEL(N, T, S) template<> const uint GaussianKernel<N>::Total = T; \
	template<> const uint GaussianKernel<N>::Vert[N] = BRACED_INIT_LIST(S); \
	template<> const uint GaussianKernel<N>::Horz[N] = BRACED_INIT_LIST(S); \
	template<> const uint GaussianKernel<N>::Matrix[N][N]

// Gaussian 3x3 with stddev = 1
GAUSSIAN_KERNEL(3, 16, (1u,2u,1u)) =	  {	{ 1,  2,  1},
											{ 2,  4,  2},
											{ 1,  2,  1}, };

// Gaussian 4x4 with stddev = 1
GAUSSIAN_KERNEL(4, 64, (1u,3u,3u,1u)) =	  {	{ 1,  3,  3,  1},
											{ 3,  9,  9,  3},
											{ 3,  9,  9,  3},
											{ 1,  3,  3,  1}, };

// Gaussian 5x5 with stddev = ~1.1
GAUSSIAN_KERNEL(5, 256, (1u,4u,6u,4u,1u)) =	  {	{ 1,  4,  6,  4,  1},
												{ 4, 16, 24, 16,  4},
												{ 6, 24, 36, 24,  6},
												{ 4, 16, 24, 16,  4},
												{ 1,  4,  6,  4,  1}, };

template<uint windowSize>
struct GaussianFilter : public SepConvFilter<windowSize, GaussianKernel<windowSize> > { };

// Gaussian 2x2 with any stddev is equivalent to mean filter
template<>
struct GaussianFilter<2> : public MeanFilter<2> { };

template<uint windowSize>
struct GaussianBinner : public WindowBinner<GaussianFilter<windowSize> > { };

#pragma endregion

#pragma region Edge Detection Functions

template<uint windowSize, typename kernel>
struct EdgeFilter : public WindowFilter<windowSize>
{
	// These are double-convolution filters with weird scaling, so we don't actually extend from SepConvFilter
	typedef kernel Kernel;

	// Since edge filters do so poorly along the image edges simply set edges to a not-very-edge-like-value
	virtual byte FilterWindow(const byte**, uint, uint, uint, uint) { return 224; }
	virtual byte FilterWindow(const byte **window)
	{
		// TODO: implement separability? however due to scaling that might not work
		int Gx = 0, Gy = 0;
		for (uint y = 0; y < windowSize; ++y)
			for (uint x = 0; x < windowSize; ++x)
			{
				byte val = window[y][x];
				Gx += Kernel::Matrix[x][y] * val;
				Gy += Kernel::Matrix[y][x] * val;
			}
		//return ~(byte)sqrt((double)((Gx*Gx + Gy*Gy) / (Kernel::Total / (255 * 255))); // real way to scale the data, but the slope is way too slow for real data
		//return ~(byte)sigmoid<Kernel::Total / 8, 0, -255>(Gx*Gx + Gy*Gy); // only top half of curve (shaped more like the sqrt)
		return ~(byte)sigmoid<0, 255>(Kernel::Total / 50, Kernel::Total / 8, Gx*Gx + Gy*Gy); // TODO: magic numbers 50 and 8
	}
};

template<uint windowSize> struct SobelKernel { static const int Matrix[windowSize][windowSize], Total; };
#define SOBEL_KERNEL(N, T) template<> const int SobelKernel<N>::Total = T; template<> const int SobelKernel<N>::Matrix[N][N]
SOBEL_KERNEL(3, 1300500) = {	{-1, 0, 1},
								{-2, 0, 2},
								{-1, 0, 1}, };
SOBEL_KERNEL(5, 170885700) = {	{-1,  -2,  0,  2, 1}, // 131.4x total from Sobel-3
								{-4,  -8,  0,  8, 4},
								{-6, -12,  0, 12, 6},
								{-4,  -8,  0,  8, 4},
								{-1,  -2,  0,  2, 1}, };
template<uint windowSize>
struct SobelFilter : public EdgeFilter<windowSize, SobelKernel<windowSize> > { };

// Scharr should be slightly more accurate then Sobel (3px) (although I can't really tell the difference)
template<uint windowSize> struct ScharrKernel { static const int Matrix[windowSize][windowSize], Total; };
#define SCHARR_KERNEL(N, T) template<> const int ScharrKernel<N>::Total = T; template<> const int ScharrKernel<N>::Matrix[N][N]
SCHARR_KERNEL(3, 23148900) = {	{ -3, 0,  3},  // 17.8x total from Sobel
								{-10, 0, 10},
								{ -3, 0,  3}, };
template<uint windowSize>
struct ScharrFilter : public EdgeFilter<windowSize, ScharrKernel<windowSize> > { };

template<uint windowSize>
struct CannyFilter : Filter<windowSize>
{
	// Canny is quite advanced and takes longer than any of the other algorithms
	// This implementation of Canny uses a higher accuracy linear interpolation during non-maximal
	// suppression and Otsu's Multi-threshold Method for automatically determining the threshold
	// values. It goes through the image 4-5 times in multiple stages of processing. It doesn't do
	// any noise reduction itself, so it is recommend to either reduce pixels or filter the image
	// also. This uses 2 image-sized double matrices of temporaries (so 8*2*w*h bytes, or 4MB for
	// 512x512 blocks) along with a 256x256 double matrix (512kb) temporary.
	//
	// TODO: the performance for use with livewire of this wasn't as good as was hoped. Some ideas
	// that were never tested were using a 5x5 Sobel (just change the window size), check that the
	// interpretation of angles is correct, messing with Otsu's method to make sure we are getting
	// good values, and not giving a completely binary image but integrating the G_max data with
	// the YES and MAYBE values.
	virtual void Run(const uint w, const uint h, const uint stride, const byte *in, byte *out, bool)
	{
		// We always set the edges to NO (even if wholesOnly) since it makes the algorithm easier to implement

		const int RADIUS = windowSize / 2;

		// Step 1: Gradient Magnitude and Direction
		// We use the Sobel filter for this
		double* G_mag = (double*)BlockPool::Get(w*h*sizeof(double));
		double* G_dir = (double*)BlockPool::Get(w*h*sizeof(double));
		memset(G_mag,                0, w * sizeof(double) * RADIUS);  // zero-out top edge
		memset(G_mag+(h-RADIUS-1)*w, 0, w * sizeof(double) * RADIUS);  // zero-out bottom edge
		double G_mag_max = 0;
		for (uint y = RADIUS; y < h - RADIUS; ++y)
		{
			const uint yw = y*w, ys = y*stride;
			for (uint x = 0; x < RADIUS; ++x) { G_mag[yw+x] = 0; G_mag[yw+w-x-1] = 0; } // zero-out left and right edges
			for (uint x = RADIUS; x < w - RADIUS; ++x)
			{
				const uint I = ys + x;
				int Gx = 0, Gy = 0;
				for (uint i = 0; i < windowSize; ++i)
				{
					const uint in_off = I + (i-RADIUS)*stride - RADIUS;
					for (uint j = 0; j < windowSize; ++j)
					{
						byte val = in[in_off + j];
						Gx += SobelKernel<windowSize>::Matrix[i][j] * val;
						Gy += SobelKernel<windowSize>::Matrix[j][i] * val;
					}
				}

				// TODO: currently the Gx and Gy are both -1* from the Python code, but this doesn't seem to matter due to symmetry
				const double mag = hypot(Gx, Gy); // = sqrt(Gx*Gx + Gy*Gy)
				G_mag[yw + x] = mag;
				G_dir[yw + x] = atan2((double)Gy, (double)Gx); // -> - pi to pi, I believe the circle has 0 to the right and pi / 2 to the bottom of the image   [TODO: check]
				if (mag > G_mag_max) { G_mag_max = mag; }
			}
		}

		// Step 2: Non-maximum Suppression
		const int OFFS[] = { 1, 1+(int)w, (int)w, (int)w-1, -1, }; // the offsets for the different 'quarters'
		const byte YES = 0x00, NO = 0xFF, MAYBE = 0x80;
		memset(out, NO, stride * h * sizeof(byte));
		// We also calculate the histogram at this point
		const static uint NBINS = 256;
		const double bin_size = NBINS / G_mag_max;
		int hist[NBINS];
		memset(hist, 0, sizeof(hist));
		for (uint y = RADIUS; y < h - RADIUS; ++y)
		{
			const uint yw = y*w, ys = y*stride;
			for (uint x = RADIUS; x < w - RADIUS; ++x)
			{
				const uint I = yw + x;
				const double dir = G_dir[I] + M_PI; // we need the direction to always be positive, it doesn't matter that we rotated 180 degrees b/c of symmetry
				double qrtr_dbl;
				const double rem = modf(dir / PI_4, &qrtr_dbl) / PI_4; // how far along we are to the next quarter (clockwise) from the quarter start
				const int quarter = (int)(qrtr_dbl + 0.1) % 4; // start of quarter 0-3
				const double mag = G_mag[I], rem1 = 1 - rem;
				const int d1 = OFFS[quarter], d2 = OFFS[quarter+1];
				if ((mag > (G_mag[I+d1]*rem1 + G_mag[I+d2]*rem)) &&
					(mag > (G_mag[I-d1]*rem1 + G_mag[I-d2]*rem)))
				{
					out[ys + x] = YES;
					// Count the value in the histogram (TODO: or should ALL pixels be considered (adding a bunch of zeros))
					int bin = (int)(mag * bin_size);
					if (bin == NBINS) { --bin; }
					++hist[bin];
				}
			}
		}
		BlockPool::Return(G_dir);

		// Step 3: Calculate and Apply Double Threshold 
		// Calculate the thresholds using Otsu's Multi-threshold Method with 2 thresholds
		uint t1 = 0, t2 = 0;
		{
			int *P = hist; // we just reuse the histogram memory for the P values (which is the cumulative sum of hist)
			int S[NBINS] = { 0 };
			for (uint i = 1; i < NBINS; ++i)
			{
				S[i] = S[i-1] + hist[i] * i;
				P[i] += P[i-1]; //P[i] = P[i-1] + hist[i]; - but P == hist
			}
			double *H = (double*)memset(BlockPool::Get(NBINS*NBINS*sizeof(double)), 0, NBINS*NBINS*sizeof(double));
			for (uint i = 1; i < NBINS; ++i)
			{
				for (uint j = i; j < NBINS; ++j)
				{
					int Pij = P[j] - P[i-1];
					if (Pij != 0)
					{
						int Sij = S[j] - S[i-1];
						H[i*NBINS+j] = ((double)Sij*Sij) / Pij;
					}
				}
			}
			double max_sigma = 0;
			for (uint i = 1; i < NBINS - 3; ++i)
			{
				for (uint j = i + 1; j < NBINS - 2; ++j)
				{
					double sigma = H[1*NBINS + i] + H[(i+1)*NBINS + j] + H[(j+1)*NBINS + NBINS-1]; // TODO: why is it H[1] and not H[0]?
					if (sigma > max_sigma) { t1 = i; t2 = j; max_sigma = sigma; }
				}
			}
			BlockPool::Return(H);
		}
		const double T_YES = t2 / bin_size;
		const double T_NO  = t1 / bin_size;

		// TOOD: possibly do some correction of maybes in this loop so we have less work later
		std::deque<int> q; // queue of YESes
		for (uint y = RADIUS; y < h - RADIUS; ++y)
		{
			const uint yw = y*w, ys = y*stride;
			for (uint x = RADIUS; x < w - RADIUS; ++x)
			{
				const uint I = ys + x;
				if (out[I] != NO)
				{
					const double mag = G_mag[yw + x];
					if (mag <= T_NO) { out[I] = NO; }
					else if (mag < T_YES) { out[I] = MAYBE; }
					else { q.push_back((int)I); }
				}
			}
		}
		BlockPool::Return(G_mag);

		// Step 4: Edge Tracking by Hysteresis
		// Promote MAYBEs to YES if next to a YES
		while (!q.empty())
		{
			int I = q.front();
			q.pop_front();
			for (int i = -RADIUS; i <= RADIUS; ++i)
			{
				int out_off = I + i*stride;
				for (int j = -RADIUS; j <= RADIUS; ++j)
				{
					if (out[out_off + j] == MAYBE)
					{
						out[out_off + j] = YES;
						q.push_back(out_off + j);
					}
				}
			}
		}
		// Set all remaining MAYBEs to NO
		for (uint y = RADIUS; y < h - RADIUS; ++y)
		{
			const uint ys = y*stride;
			for (uint x = RADIUS; x < w - RADIUS; ++x)
			{
				const uint I = ys + x;
				if (out[I] == MAYBE) { out[I] = NO; }
			}
		}
	}
};

#pragma endregion

#pragma region Accentuation Functions
/// <summary>Calculates a sigmoid function for the given value with a slope of 0.05 and half-max at 128.</summary>
/// <param name="x">The value of x to calculate the sigmoid for (0-255)</param>
/// <returns>The evaluated sigmoid function value (0-255)</returns>
inline static byte sigmoid_accentuation(byte x) { return sigmoid<0, 255>(20, 128, x); }

//inline static double sigmoid(float x, double halfmax_x, double max_y, double slope) { return max_y / (1 + exp(slope * (halfmax_x - x))); }
//inline static double sigmoid(float x) { return sigmoid(x, 0.5, 1.0, 10.0); }

struct SigmoidAccentuationFilter : public PixelFilter { virtual byte FilterPixel(byte x) { return sigmoid<0, 255>(20, 128, x); } };
#pragma endregion

struct Invert : public PixelFilter { virtual byte FilterPixel(byte x) { return ~x; } };

inline static void Swap(byte *&a, byte *&b) { byte *x = a; a = b; b = x; }

#define WINDOW_SIZE(m) ((m) & 0xF) // works for PixelReductionMethod, NoiseReductionMethod, and EdgeDetectionMethod

#define STATUS_NOT_DONE	0 // have not done the block and it isn't in the queue
#define STATUS_WILL_DO	1 // have not done the block but it is in the queue
#define STATUS_DOING	2 // actively computing the block and no longer in the queue
#define STATUS_DONE		3 // done computing the block / it is usable
// TODO: status really only takes up 2 bits, can fit 4 in a byte

#define BLOCK_SIZE	0x100
#define BLOCK_SIZE_	LOG2(BLOCK_SIZE)
CASSERT(IS_POWER_OF_TWO(BLOCK_SIZE)); // must be a power of 2, there are numerous mathematical shortcuts taken with this assumption in mind

Weights::Settings::Settings(CoalescingMethod Method, bool Invert, PixelReductionMethod PixelReduction, NoiseReductionMethod NoiseReduction, AccentuationMethod Accentuation, EdgeDetectionMethod EdgeDetection) :
	Method(Method), Invert(Invert), PixelReduction(PixelReduction), NoiseReduction(NoiseReduction), Accentuation(Accentuation), EdgeDetection(EdgeDetection) {}
bool Weights::Settings::operator ==(const Weights::Settings& r) const
{
	return this->Method == r.Method && this->PixelReduction == r.PixelReduction && this->NoiseReduction == r.NoiseReduction && this->EdgeDetection == r.EdgeDetection && this->Accentuation == r.Accentuation && this->Invert == r.Invert;
}

const Weights::Settings Weights::GrayscaleSettings;
const Weights::Settings Weights::ColorSettings(Weights::WeightedHSV, false, Weights::Mean2pxWindow, Weights::NoNoiseReduction, Weights::Sigmoid, Weights::Sobel5);

inline static uint CalcFilterOverflow(const Weights::Settings& settings)
{
	return WINDOW_SIZE(settings.NoiseReduction) / 2 + WINDOW_SIZE(settings.EdgeDetection) / 2;
}

Weights::Weights() : Threaded("Weights"),
	_data_raw(NULL), _data(NULL), _status(NULL), _width_raw(0), _height_raw(0), _stride(0), _width(0), _height(0),
	_scale(WINDOW_SIZE(this->_settings.PixelReduction)), _filter_overflow(CalcFilterOverflow(this->_settings)) { }
Weights::~Weights() { this->Stop(); free(this->_data); free(this->_status); }

inline void Weights::UpdatedScaleOrSize()
{
	if (this->_width_raw || this->_height_raw)
	{
		this->_width = ScaleBack(this->_width_raw, WINDOW_SIZE(this->_scale));
		this->_height = ScaleBack(this->_height_raw, WINDOW_SIZE(this->_scale));
		this->_width_status = ScaleBack<BLOCK_SIZE>(this->_width);
		this->_height_status = ScaleBack<BLOCK_SIZE>(this->_height);
		this->_last_col_width = this->_width - ((this->_width_status-1) << BLOCK_SIZE_);

		this->_block_queue.SetSize(this->_width_status, this->_height_status);

		const uint wh_s = this->_width_status*this->_height_status;
		this->_data   = (byte**)memset(realloc(this->_data, wh_s*sizeof(byte*)), NULL, wh_s*sizeof(byte*));
		this->_status = (byte*)memset(realloc(this->_status, wh_s), STATUS_NOT_DONE, wh_s);
		this->SetTotalProgress(wh_s);
	}
}

void Weights::SetSettings(const Settings& settings)
{
	const byte *data_raw = this->_data_raw;

	this->Stop();

	this->_settings = settings;
	this->_filter_overflow = CalcFilterOverflow(settings);
	bool size_changed = this->_scale != (uint)WINDOW_SIZE(settings.PixelReduction);
	if (size_changed)
	{
		this->_scale = WINDOW_SIZE(settings.PixelReduction);
		this->UpdatedScaleOrSize();
	}

	if (data_raw)
	{
		this->_data_raw = data_raw;
		this->Start();
	}

	emit SettingsChanged(size_changed, false);
}
const Weights::Settings& Weights::GetSettings() const { return this->_settings; }

void Weights::SetImage(const byte* imageData, uint W, uint H, DataFormat format, uint stride)
{
	this->Stop();

	bool size_changed = this->_width_raw != W || this->_height_raw != H;
	if (size_changed)
	{
		this->_width_raw = W;
		this->_height_raw = H;
		this->UpdatedScaleOrSize();
	}

	this->_stride = stride;
	this->_format = format;
	this->_data_raw = imageData;

	this->Start();

	emit SettingsChanged(size_changed, true);
}

bool Weights::Get(uint X, uint Y, byte* w) const
{
	if (!this->_status) // no image set at all
		return false;
	const uint x = X >> BLOCK_SIZE_, i = (Y >> BLOCK_SIZE_) * this->_width_status + x;
	if (this->_status[i] == STATUS_WILL_DO || this->_status[i] == STATUS_DOING)
	{
		this->_status_lock.lock();
		while (this->_status[i] == STATUS_WILL_DO || this->_status[i] == STATUS_DOING)
			this->_block_finished.wait(&this->_status_lock, 250);
		this->_status_lock.unlock();
	}
	if (this->_status[i] == STATUS_NOT_DONE) // block not complete
		return false;
	const uint x_ = X & (BLOCK_SIZE - 1), y_ = Y & (BLOCK_SIZE - 1); // X % d (where d is a power of 2) = X & (d - 1)
	*w = this->_data[i][((x == this->_width_status-1) ? (y_ * this->_last_col_width) : (y_ << BLOCK_SIZE_)) + x_];
	//*w = this->_data[i][(y_ + 1) * (2 + ((x == this->_width_status - 1) ? this->_last_col_width : BLOCK_SIZE)) + x_ + 1];
	return true;
}

void Weights::Stop()
{
	if (!this->_status) // no image set at all so already 'stopped' 
		return;

	emit this->Stopping();

	const uint wh = this->_width_status*this->_height_status;
	
	if (this->IsExecuting())
	{
		QMutexLocker(&this->_status_lock);
		QMutexLocker(&this->_queue_lock);

		Threaded::Stop(false);

		this->_block_queue.Clear();
		this->_blocks_queued.wakeAll();

		memset(this->_status, STATUS_NOT_DONE, wh);
		this->_block_finished.wakeAll();
	}

	this->_data_raw = NULL;
	for (uint I = 0; I < wh; ++I)
		if (this->_data[I])
		{
			BlockPool::Return(this->_data[I]);
			this->_data[I] = NULL;
		}

	this->wait();
}

inline static uint CalcScore(uint x, uint y, double X, double Y)
{
	uint _X = (uint)X, _Y = (uint)Y;
	if (x == _X) return (y == _Y) ? 0 : (((y < _Y) ? (Y - y - 1) : (y - Y)) * BLOCK_SIZE);
	if (y == _Y) return                  ((x < _X) ? (X - x - 1) : (x - X)) * BLOCK_SIZE;

	double dx = (x < _X) ? (X - x - 1) : (x - X);
	double dy = (y < _Y) ? (Y - y - 1) : (y - Y);
	return (uint)(sqrt(dx*dx+dy*dy) * BLOCK_SIZE);
}
typedef struct _dblpoint { double x, y; } dblpoint;
static uint CalcScoreCB(uint x, uint y, uint, dblpoint *pt) { return CalcScore(x, y, pt->x, pt->y); }

void Weights::CalculateRegion(uint x, uint y, uint min_room)
{
	//assert(this->_data_raw != NULL);

	const uint l = x > min_room ? (x - min_room) >> BLOCK_SIZE_ : 0;
	const uint t = y > min_room ? (y - min_room) >> BLOCK_SIZE_ : 0;
	const uint r = x + min_room < this->_width  ? (x + min_room) >> BLOCK_SIZE_ : this->_width_status  - 1;
	const uint b = y + min_room < this->_height ? (y + min_room) >> BLOCK_SIZE_ : this->_height_status - 1;
	const double X = x / (double)BLOCK_SIZE, Y = y / (double)BLOCK_SIZE;
	const dblpoint pt = {X, Y};

	uint I = t * this->_width_status + l, i;

	bool added = false;

	QMutexLocker locker1(&this->_status_lock);
	QMutexLocker locker2(&this->_queue_lock);

	this->_block_queue.UpdateAllScores((PointPriorityQueue::CalcScore)CalcScoreCB, &pt);

	for (y = t; y <= b; I += this->_width_status, ++y)
	{
		for (x = l, i = I; x <= r; ++x, ++i)
		{
			if (this->_status[i] == STATUS_NOT_DONE)
			{
				// add to the queue
				this->_block_queue.Enqueue(x, y, i, CalcScore(x, y, X, Y));
				this->_status[i] = STATUS_WILL_DO;
				added = true;
			}
		}
	}
		
	if (added)
		this->_blocks_queued.wakeAll();
}

void Weights::Run()
{
	// TODO: make multi-threaded
	// int threadCount = QThread::idealThreadCount(); // equals the number of processor cores

	uint x, y, I, score;
	bool not_empty;

	do
	{
		for (;;)
		{
			this->_queue_lock.lock();
			not_empty = this->_block_queue.Dequeue(x, y, I, score);
			this->_queue_lock.unlock();
			
			if (!this->IsExecuting()) return;
			if (!not_empty) break;

			this->CalcBlock(x, y, I);
			this->IncProgress();
		}

		this->_queue_lock.lock();
		if (this->IsExecuting())
			this->_blocks_queued.wait(&this->_queue_lock);
		this->_queue_lock.unlock();
	}
	while (this->IsExecuting());
}

void Weights::CalcBlock(uint x_s, uint y_s, uint I)
{
	this->_status_lock.lock();
	this->_status[I] = STATUS_DOING;
	this->_status_lock.unlock();

	const uint w = this->_width,     h = this->_height;
	const uint W = this->_width_raw, H = this->_height_raw;
	const uint scale = this->_scale, extra = this->_filter_overflow;

	const bool short_right = x_s == this->_width_status-1, short_bottom = y_s == this->_height_status-1;
	const bool has_l = (bool)x_s, has_t = (bool)y_s, has_r = !short_right, has_b = !short_bottom;
	const bool wholesOnly = has_l && has_t && has_r && has_b;
	const uint l_extra = has_l * extra, t_extra = has_t * extra;

	// post-scaling dimensions (they shrink as we apply filters)
	uint bw = (short_right  ? this->_last_col_width : (BLOCK_SIZE + extra)) + l_extra;
	uint bh = (short_bottom ? h & (BLOCK_SIZE - 1)  : (BLOCK_SIZE + extra)) + t_extra;
	const uint stride = bw;
	uint off = 0; // updated as we apply filters

	// pre-scaling dimensions
	const uint BW = short_right  ? W % (scale << BLOCK_SIZE_) + l_extra * scale : (bw * scale);
	const uint BH = short_bottom ? H % (scale << BLOCK_SIZE_) + t_extra * scale : (bh * scale);

	const byte *in = this->_data_raw + (((y_s << BLOCK_SIZE_) - t_extra) * this->_stride + ((x_s << BLOCK_SIZE_) - l_extra) * this->_format) * scale;
	byte *out  = (byte*)BlockPool::Get(BW*BH);
	byte *temp = (byte*)BlockPool::Get(bw*bh);
	byte *big = out, *small = temp; // if scale > 1, big is the bigger of the two memory allocations (otherwise they are the same size)

	switch (this->_format)
	{
	case GrayscaleByte:   CoalesceGrayByte  (BW, BH, this->_stride, in, out); break;
	case GrayscaleUShort: CoalesceGrayUShort(BW, BH, this->_stride, in, out); break;
	case RGB:
		switch (this->_settings.Method)
		{
		case RedChannel:	CoalesceRed			(BW, BH, this->_stride, in, out); break;
		case GreenChannel:	CoalesceGreen		(BW, BH, this->_stride, in, out); break;
		case BlueChannel:	CoalesceBlue		(BW, BH, this->_stride, in, out); break;
		case AvgRGB:		CoalesceAvgRGB		(BW, BH, this->_stride, in, out); break;
		case Luma:			CoalesceLuma		(BW, BH, this->_stride, in, out); break;
		case Luma601:		CoalesceLuma601		(BW, BH, this->_stride, in, out); break;
		case LumaSMPTE:		CoalesceLumaSMPTE	(BW, BH, this->_stride, in, out); break;
		case WeightedHSV:	CoalesceWeightedHSV	(BW, BH, this->_stride, in, out); break;
		case WeightedHSL:	CoalesceWeightedHSL	(BW, BH, this->_stride, in, out); break;
		case WeightedHSI:	CoalesceWeightedHSI	(BW, BH, this->_stride, in, out); break;
		}
	}

	if (this->_settings.Invert) { Invert().Run(BW, BH, BW, out, out, wholesOnly); }

	switch (this->_settings.PixelReduction)
	{
	case NoPixelReduction: break; // do nothing
	case Median2pxWindow:	MedianBinner<2>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Median3pxWindow:	MedianBinner<3>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Median4pxWindow:	MedianBinner<4>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Median5pxWindow:	MedianBinner<5>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Mean2pxWindow:		MeanBinner<2>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Mean3pxWindow:		MeanBinner<3>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Mean4pxWindow:		MeanBinner<4>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Mean5pxWindow:		MeanBinner<5>().	Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Gaussian3pxWindow:	GaussianBinner<3>().Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Gaussian4pxWindow:	GaussianBinner<4>().Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	case Gaussian5pxWindow:	GaussianBinner<5>().Run(BW, BH, BW, out, temp); Swap(out, temp); break;
	}

	switch (this->_settings.NoiseReduction)
	{
	case NoNoiseReduction: break; // do nothing
	case MedianFilter3pxWindow:		MedianFilter<3>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case MedianFilter5pxWindow:		MedianFilter<5>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case MeanFilter3pxWindow:		MeanFilter<3>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case MeanFilter5pxWindow:		MeanFilter<5>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case GaussianFilter3pxWindow:	GaussianFilter<3>().Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case GaussianFilter5pxWindow:	GaussianFilter<5>().Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	}

	{
		const uint ws = WINDOW_SIZE(this->_settings.NoiseReduction) / 2;
		bw -= (has_l + has_r) * ws;
		bh -= (has_t + has_b) * ws;
		off += (has_l + has_t * stride) * ws;
	}

	switch (this->_settings.Accentuation)
	{
	case NoAccentuation: break; // do nothing
	case Sigmoid:	SigmoidAccentuationFilter().Run(bw, bh, stride, out + off, out + off, wholesOnly); break;
	}

	//{
	//	const uint ws = WINDOW_SIZE(this->_settings.Accentuation) / 2;
	//	bw -= (has_l + has_r) * ws;
	//	bh -= (has_t + has_b) * ws;
	//	off += (has_l + has_t * stride) * ws;
	//}

	switch (this->_settings.EdgeDetection)
	{
	case NoEdgeDetection: break; // do nothing
	case Sobel3:		SobelFilter<3>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case Sobel5:		SobelFilter<5>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case Scharr:		ScharrFilter<3>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	case Canny:			CannyFilter<3>().	Run(bw, bh, stride, out + off, temp + off, wholesOnly); Swap(out, temp); break;
	}

	{
		const uint ws = WINDOW_SIZE(this->_settings.EdgeDetection) / 2;
		bw -= (has_l + has_r) * ws;
		bh -= (has_t + has_b) * ws;
		off += (has_l + has_t * stride) * ws;
	}

	// Remove any padding due to stride != bw and move all data into the smaller memory
	for (uint i = 0; i < bh; ++i) { memmove(small + i*bw, out + i*stride + off, bw); }
	this->_data[I] = small;
	BlockPool::Return(big);

	this->_status_lock.lock();
	this->_status[I] = STATUS_DONE;
#ifdef _DEBUG
	qDebug("[%p] %u done", QThread::currentThreadId(), I);
#endif
	this->_block_finished.wakeAll();
	this->_status_lock.unlock();

#ifdef SAVE_WEIGHT_IMAGE
	this->Checkpoint("saving weights image");
	WriteBitmap(this->_data, w, h, BLOCK_SIZE, "weights");
#endif
}

uint Weights::GetBlocksCalculated(QVector<QPoint>& done, QVector<QPoint>& doing, QVector<QPoint>& not_done) const
{
	if (!this->_status) return 0;
	const uint BS = this->_scale << BLOCK_SIZE_;
	const uint w = this->_width_status * BS, h = this->_height_status * BS;
	for (uint y = 0, I = 0; y < h; y += BS)
		for (uint x = 0; x < w; x += BS, ++I)
		{
			if (this->_status[I] == STATUS_DONE)												done.push_back(QPoint(x, y));
			else if (this->_status[I] == STATUS_DOING || this->_status[I] == STATUS_WILL_DO)	doing.push_back(QPoint(x, y));
			else if (this->_status[I] == STATUS_NOT_DONE)										not_done.push_back(QPoint(x, y));
		}
	return BS;
}
#ifndef IMOD_PLUGIN
void Weights::SaveImage(const char *name) const
{
	const uint total = this->_width_status * this->_height_status;
	this->_status_lock.lock();
	for (uint i = 0; i < total; ++i)
	{
		while (this->_status[i] == STATUS_WILL_DO || this->_status[i] == STATUS_DOING)
			this->_block_finished.wait(&this->_status_lock);
	}
	this->_status_lock.unlock();
	WriteBitmap(this->_data, this->_width, this->_height, BLOCK_SIZE, name);
}
#endif
