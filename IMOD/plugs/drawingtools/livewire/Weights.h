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

#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "general.h"

#include <QMutex>
#include <QWaitCondition>

#include <QPoint>
#include <QVector>

#include "Threaded.h"
#include "PointPriorityQueue.h"

namespace Livewire
{
	class Weights : public Threaded
	{
		Q_OBJECT

	public:
		enum DataFormat
		{
			GrayscaleByte = 1,
			GrayscaleUShort = 2,
			RGB = 4,
		};

		enum CoalescingMethod
		{
			RedChannel, GreenChannel, BlueChannel,
			AvgRGB,
			Luma, Luma601, LumaSMPTE, // Rec. 709, Rec. 601, or SMPTE 240M
			WeightedHSV, WeightedHSL, WeightedHSI, // 0.6 * [V or L or I] + 0.3 * H + 0.1 * S
		};
		enum PixelReductionMethod
		{
			NoPixelReduction  = 0x01,
			Median2pxWindow   = 0x12, Median3pxWindow   = 0x13, Median4pxWindow   = 0x14, Median5pxWindow   = 0x15,
			Mean2pxWindow     = 0x22, Mean3pxWindow     = 0x23, Mean4pxWindow     = 0x24, Mean5pxWindow     = 0x25,
			/* same as mean 2x2 */    Gaussian3pxWindow = 0x33, Gaussian4pxWindow = 0x34, Gaussian5pxWindow = 0x35,
			// TODO: any noise reduction methods added
		};
		enum NoiseReductionMethod
		{
			NoNoiseReduction        = 0x01,
			MedianFilter3pxWindow   = 0x13, MedianFilter5pxWindow   = 0x15,
			MeanFilter3pxWindow     = 0x23, MeanFilter5pxWindow     = 0x25,
			GaussianFilter3pxWindow = 0x33, GaussianFilter5pxWindow = 0x35,
			// TODO: AnisotropicDiffusion, http://en.wikipedia.org/wiki/Anisotropic_diffusion
			// TODO: Other RCRS filters: max, min, LUM, http://en.wikipedia.org/wiki/Noise_reduction#Nonlinear_filters
			// TODO: k-Nearest Neighbor Filtering, http://www.anirudh.net/courses/cse585/project1/
		};
		enum AccentuationMethod
		{
			NoAccentuation = 0x01,
			Sigmoid = 0x11,
			// TODO: histogram equalization [+exact and/or CLAHE]
			// TODO: linear
		};
		enum EdgeDetectionMethod
		{
			NoEdgeDetection = 0x01,
			Sobel3  = 0x13, Sobel5  = 0x15,
			Scharr  = 0x23,
			Canny   = 0x33,
		};

		struct Settings
		{
			CoalescingMethod Method;
			bool Invert;
			PixelReductionMethod PixelReduction;
			NoiseReductionMethod NoiseReduction;
			AccentuationMethod Accentuation;
			EdgeDetectionMethod EdgeDetection;

			Settings(
				CoalescingMethod Method = BlueChannel,
				bool Invert = false,
				PixelReductionMethod PixelReduction = Mean2pxWindow,
				NoiseReductionMethod NoiseReduction = NoNoiseReduction,
				AccentuationMethod Accentuation = Sigmoid,
				EdgeDetectionMethod EdgeDetection = NoEdgeDetection
				);

			bool operator ==(const Settings& rhs) const;
			inline bool operator !=(const Settings& rhs) const { return !this->operator==(rhs); }
		};

		static const Settings GrayscaleSettings; // the default settings
		static const Settings ColorSettings; // HSV + Sobel5 (in addition to the grayscale settings)

	private:
		Settings _settings;
		uint _filter_overflow;

		uint _width_raw, _height_raw;
		uint _stride;
		DataFormat _format;
		const byte *_data_raw;

		uint _width, _height, _scale;
		byte **_data;

		uint _width_status, _height_status, _last_col_width;
		byte *_status;
		PointPriorityQueue _block_queue;

		mutable QMutex _queue_lock, _status_lock;
		mutable QWaitCondition _blocks_queued, _block_finished;

		// prevent copying
		Weights(const Weights&);
		Weights& operator=(const Weights&);

		void UpdatedScaleOrSize();

	public:
		Weights();
		virtual ~Weights();
		
		void SetImage(const byte* image, uint W, uint H, DataFormat format, uint stride);
		void SetSettings(const Settings& settings);
		const Settings& GetSettings() const;

		// points in the original size, not reduced
		// points are of the upper-left corner of each calculated block
		// the return value is the size of the blocks
		uint GetBlocksCalculated(QVector<QPoint>& done, QVector<QPoint>& doing, QVector<QPoint>& not_done) const;

		// Gets a weight
		// All indices are of the reduced size
		// x = x coord
		// y = y coord
		// w = (out) the weight
		// returns true if the point is obtainable, false if the point is out of bounds of the calculation region
		// if the requested point is not yet calculated but has been requested this function will block until it is ready
		bool Get(uint x, uint y, byte* w) const;

#ifndef IMOD_PLUGIN
		// will block while there are pending pixels to be calculated then save the weights as an image
		void SaveImage(const char* name) const;
#endif
		
		inline uint GetScale()          const { return this->_scale;      }

		inline uint GetOriginalWidth()  const { return this->_width_raw;  }
		inline uint GetOriginalHeight() const { return this->_height_raw; }

		inline uint GetReducedWidth()   const { return this->_width;      }
		inline uint GetReducedHeight()  const { return this->_height;     }

		void CalculateRegion(uint x, uint y, uint min_room); // indices of the reduced size
		virtual void Stop();
	protected:
		void Run();

	signals:
		void Stopping();
		void SettingsChanged(bool size, bool image);

	private:
		void CalcBlock(uint x, uint y, uint I);

	};
}

#endif
