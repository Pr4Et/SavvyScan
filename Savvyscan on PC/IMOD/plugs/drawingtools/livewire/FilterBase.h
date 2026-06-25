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

// Includes the base filter classes but no actual filters
// Classes:
//   Filter<windowSize>        -  the base of all filters
//   PixelFilter               -  a generic filter that operates on individual pixels
//   WindowFilter<windowSize>  -  a generic filter that operates on windows of pixels 
//   RCRSFilter<windowSize>    -  a rank-conditioned rank-selection filter
//   ConvFilter<windowSize>    -  a convolution-based filter
//   SepConvFilter<windowSize> -  a convolution-based filter where the kernel is separable
//   WindowBinner<windowSize>  -  wraps a WindowFilter and adds the ability to bin

#ifndef FILTERBASE_H
#define FILTERBASE_H

#include "general.h"

template<uint windowSize>
struct Filter
{
	static const uint WindowSize = windowSize;
	virtual void Run(const uint w, const uint h, const uint stride, const byte *in, byte *out, bool wholesOnly) = 0;
	/*static const uint Matrix[windowSize][windowSize];*/
	virtual ~Filter() { }
};

struct PixelFilter : public Filter<1>
{
	virtual void Run(const uint w, const uint h, const uint stride, const byte *in, byte *out, bool)
	{
		for (uint y = 0; y < h; ++y)
		{
			const uint I = y*stride;
			for (uint x = 0; x < w; ++x)
			{
				out[I + x] = FilterPixel(in[I + x]);
			}
		}
	}
	virtual byte FilterPixel(byte p) = 0;
	virtual ~PixelFilter() { }
};

template<uint windowSize> // windowSize must be odd and greater than 1 (except for use with WindowBinner)
struct WindowFilter : public Filter<windowSize>
{
	virtual void Run(const uint w, const uint h, const uint stride, const byte *in, byte *out, bool wholesOnly)
	{
		// Make sure that the window size is >1 and odd
		// TODO: support this check only when this function is actually called (right now binners triggers this error)
		//CASSERT(windowSize != 1 && (windowSize & 1) == 1);

		static const uint WS1 = windowSize - 1, ws = windowSize / 2, ws1 = ws + 1;
		const uint w1 = w - 1, w_ws = w - ws, ws_w = w + ws, /*h1 = h - 1,*/ h_ws = h - ws, ws_h = h + ws;
		const byte *window[windowSize];

		uint X, Y, Yw, I;

#define ADVANCE_X(N) for (uint y = 0; y < N; ++y) ++window[y];
		if (wholesOnly)
		{
			// full neighborhood only
			const uint gap = stride - w + 2 * ws;
			for (Y = ws, I = ws*(stride + 1), Yw = 0; Y < h_ws; ++Y, Yw += stride, I += gap)
			{
				for (uint y = 0, i = Yw; y < windowSize; ++y, i += stride) window[y] = in + i;
				for (X = ws; X < w_ws; ++X) { out[I++] = FilterWindow(window); ADVANCE_X(windowSize); }
			}
		}
		else
		{
			const uint gap = stride - w;
			for (uint y = 0, i = 0; y < WS1; ++y, i += stride) window[y] = in + i;

			// missing left columns and top rows
			for (X = 0; X < ws; ++X) { for (Y = 0, I = X; Y < ws; ++Y, I += stride) out[I] = FilterWindow(window, ws1 + X, ws1 + Y, X, Y); }

			// missing no columns and top rows
			for (/*X = ws*/; X < w_ws; ++X) { for (Y = 0, I = X; Y < ws; ++Y, I += stride) out[I] = FilterWindow(window, windowSize, ws1 + Y, ws, Y);	ADVANCE_X(WS1); }

			// missing right columns and top rows
			for (/*X = w_ws*/; X < w1; ++X) { for (Y = 0, I = X; Y < ws; ++Y, I += stride) out[I] = FilterWindow(window, ws_w - X, ws1 + Y, ws, Y);		ADVANCE_X(WS1); }
			/*X = w1*/ for (Y = 0, I = w1; Y < ws; ++Y, I += stride) out[I] = FilterWindow(window, ws1, ws1 + Y, ws, Y);

			I += gap;
			for (/*Y = ws,*/ I = ws*stride, Yw = 0; Y < h_ws; ++Y, Yw += stride, I += gap)
			{
				for (uint y = 0, i = Yw; y < windowSize; ++y, i += stride) window[y] = in + i;

				// missing left columns and no rows
				for (X = 0; X < ws; ++X) { out[I++] = FilterWindow(window, ws1 + X, windowSize, X, ws); }

				// full neighborhood
				for (/*X = ws*/; X < w_ws; ++X) { out[I++] = FilterWindow(window);									ADVANCE_X(windowSize); }

				// missing right columns and no rows
				for (/*X = w_ws*/; X < w1; ++X) { out[I++] = FilterWindow(window, ws_w - X, windowSize, ws, ws);	ADVANCE_X(windowSize); }
				out[I++] = FilterWindow(window, ws1, windowSize, ws, ws);
			}

			for (uint y = 0, i = Yw; y < WS1; ++y, i += stride) window[y] = in + i;
			Yw = I;

			// missing left columns and bottom rows
			for (X = 0; X < ws; ++X) { for (Y = h_ws, I = X + Yw; Y < h; ++Y, I += stride) out[I] = FilterWindow(window, ws1 + X, ws_h - Y, X, ws); }

			// missing no columns and bottom rows
			for (/*X = ws*/; X < w_ws; ++X) { for (Y = h_ws, I = X + Yw; Y < h; ++Y, I += stride) out[I] = FilterWindow(window, windowSize, ws_h - Y, ws, ws);	ADVANCE_X(WS1); }

			// missing right columns and bottom rows
			for (/*X = w_ws*/; X < w1; ++X) { for (Y = h_ws, I = X + Yw; Y < h; ++Y, I += stride) out[I] = FilterWindow(window, ws_w - X, ws_h - Y, ws, ws);		ADVANCE_X(WS1); }
			/*X = w1*/ for (Y = h_ws, I = w1 + Yw; Y < h; ++Y, I += stride) out[I] = FilterWindow(window, ws1, ws_h - Y, ws, ws);
		}
#undef ADVANCE_X
	}
	virtual byte FilterWindow(const byte **window) = 0;
	virtual byte FilterWindow(const byte **window, uint w, uint h, uint cx, uint cy) = 0;
	virtual ~WindowFilter() { }
};

template<>
struct WindowFilter<3> : public Filter<3>
{
	void Run(const uint w, const uint h, const uint stride, const byte *in, byte *out, bool wholesOnly)
	{
		// Specialize the most common and smallest window filter
	
		const uint w1 = w - 1, h1 = h - 1;
		const byte *window[3];
	
		uint X, Y, Yw, I = 0;
	
		if (wholesOnly)
		{
			// full neighborhood only
			const uint gap = stride - w + 2;
			for (I = stride + 1, Y = 1, Yw = stride; Y < h1; ++Y, Yw += stride, I += gap)
			{
				window[0] = in + Yw - stride;
				window[1] = in + Yw;
				window[2] = in + Yw + stride;
				for (X = 1; X < w1; ++X, ++window[0], ++window[1], ++window[2]) { out[I++] = FilterWindow(window); }
			}
		}
		else
		{
			const uint gap = stride - w;
			window[0] = in;
			window[1] = in + stride;
			out[I++] = FilterWindow(window, 2, 2, 0, 0);																// missing left columns and top rows
			for (X = 1; X < w1; ++X) { out[I++] = FilterWindow(window, 3, 2, 1, 0); ++window[0]; ++window[1]; }			// missing no columns and top rows
			out[I++] = FilterWindow(window, 2, 2, 1, 0);																// missing right columns and top rows
	
			I += gap;
			for (/*I = stride,*/ Y = 1, Yw = stride; Y < h1; ++Y, Yw += stride, I += gap)
			{
				window[0] = in + Yw - stride;
				window[1] = in + Yw;
				window[2] = in + Yw + stride;
				out[I++] = FilterWindow(window, 2, 3, 0, 1);															// missing left columns and no rows
				for (X = 1; X < w1; ++X, ++window[0], ++window[1], ++window[2])	{ out[I++] = FilterWindow(window); }	// full neighborhood
				out[I++] = FilterWindow(window, 2, 3, 1, 1);															// missing right columns and no rows
			}
		
			window[0] = in + Yw - stride;
			window[1] = in + Yw;
			out[I++] = FilterWindow(window, 2, 2, 0, 1);																// missing left columns and bottom rows
			for (X = 1; X < w1; ++X) { out[I++] = FilterWindow(window, 3, 2, 1, 1); ++window[0]; ++window[1]; }			// missing no columns and bottom rows
			out[I  ] = FilterWindow(window, 2, 2, 1, 1);																// missing right columns and bottom rows
		}
	}
	virtual byte FilterWindow(const byte **window) = 0;
	virtual byte FilterWindow(const byte **window, uint w, uint h, uint cx, uint cy) = 0;
	virtual ~WindowFilter() { }
};

template<typename filter> // Filter::WindowSize must be greater than 1
struct WindowBinner : filter
{
	typedef filter Filter;
	static const uint WindowSize = Filter::WindowSize;
	void Run(const uint W, const uint H, const uint stride, const byte *in, byte *out)
	{
		// Make sure that the window size is >1
		CASSERT(WindowSize > 1);
	
		// W and H are the raw width and height (in)
		// w and h are the binned width and height (out)
	
		static const uint ws = WindowSize / 2; // WS1 = WindowSize - 1
	
		const uint w_f =  W        / WindowSize, h_f =  H        / WindowSize; // the number of full bins
		//const uint w   = (W + WS1) / windowSize, h   = (H + WS1) / windowSize; // binned size
		const uint w_r =  W - w_f  * WindowSize, h_r =  H - h_f  * WindowSize; // the remainder that doesn't fit in a full bin

		const byte *window[WindowSize];
	
#define ADVANCE_X(N) for (uint y = 0; y < N; ++y) window[y] += WindowSize;
	
		uint I = 0, Yw = 0;
		for (uint Y = 0; Y < h_f; ++Y)
		{
			for (uint y = 0; y < WindowSize; ++y, Yw += stride) window[y] = in + Yw;
	
			// full neighborhood
			for (uint X = 0; X < w_f; ++X)	{ out[I++] = Filter::FilterWindow(window); ADVANCE_X(WindowSize); }
	
			// missing right columns and no rows
			if (w_r)						{ out[I++] = Filter::FilterWindow(window, w_r, WindowSize, ws, ws); }
		}
	
		if (h_r)
		{
			for (uint y = 0; y < h_r; ++y, Yw += stride) window[y] = in + Yw;
	
			// missing no columns and bottom rows
			for (uint X = 0; X < w_f; ++X)	{ out[I++] = Filter::FilterWindow(window, WindowSize, h_r, ws, ws); ADVANCE_X(h_r); }
			
			// missing right columns and bottom rows
			if (w_r)						{ out[I] = Filter::FilterWindow(window, w_r, h_r, ws, ws); }
		}
	
#undef ADVANCE_X
	}
};

template<uint windowSize>
struct RCRSFilter : public WindowFilter<windowSize>
{
	static int byte_comp(const void *a, const void *b) { return ((int)*(byte*)a) - *(byte*)b; }
	virtual byte FilterWindow(const byte** window)
	{
		static const uint ws2 = windowSize*windowSize, ws_2 = windowSize / 2;
		byte list[ws2];
		for (uint y = 0, i = 0; y < windowSize; ++y, i += windowSize)
			memcpy(list + i, window[y], windowSize);
		qsort(list, ws2, 1, &byte_comp);
		return SelectValue(list, window[ws_2][ws_2]);
	}
	virtual byte FilterWindow(const byte** window, uint w, uint h, uint cx, uint cy)
	{
		const uint count = w*h;
		byte list[windowSize*windowSize];
		for (uint y = 0, i = 0; y < h; ++y, i += w)
			memcpy(list + i, window[y], w);
		qsort(list, w*h, 1, &byte_comp);
		return SelectValue(list, count, window[cx][cy]);
	}
	virtual byte SelectValue(const byte* list, const byte val) = 0;
	virtual byte SelectValue(const byte* list, uint count, const byte val) = 0;
	virtual ~RCRSFilter() { }
};

template<uint windowSize, typename kernel>
struct ConvFilter : public WindowFilter<windowSize>
{
	typedef kernel Kernel;
	virtual byte FilterWindow(const byte** window)
	{
		uint v = 0;
		for (uint y = 0; y < windowSize; ++y)
			for (uint x = 0; x < windowSize; ++x)
				v += Kernel::Matrix[x][y] * window[y][x];
		return (byte)(v / Kernel::Total);
	}
	virtual byte FilterWindow(const byte** window, uint w, uint h, uint cx, uint cy)
	{
		static const uint ws = windowSize / 2;
		const uint ws_cx = ws - cx, ws_cy = ws - cy;

		uint v = 0, t = 0;
		for (uint y = 0; y < h; ++y)
			for (uint x = 0; x < w; ++x)
			{
				uint m = Kernel::Matrix[x + ws_cx][y + ws_cy];
				v += m * window[y][x];
				t += m;
			}
		return (byte)(v / t);
	}
	virtual ~ConvFilter() { }
};

template<uint windowSize, typename kernel>
struct SepConvFilter : public ConvFilter<windowSize, kernel>
{
	// TODO: add code to do the optimized separate 1D convolutions
	// Will eventually use the Vert and Horz vectors of the kernel
	virtual ~SepConvFilter() { }
};

#endif
