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

#ifndef LIVEWIRE_H
#define LIVEWIRE_H

#include "general.h"

#include "PointPriorityQueue.h"
#include "Threaded.h"
#include "Weights.h"

#include <QPoint>
#include <QVector>

#ifndef IMOD_PLUGIN
#include <QPainter>
#endif

namespace Livewire
{
	class LivewireCalculator : public Threaded
	{
		Q_OBJECT

	private:
		/// <summary>The weights backing this livewire calculator</summary>
		Weights *_weights;

		/// <summary>If this livewire calculator owns its weights or not</summary>
		bool _owns_weights;

		/// <summary>The X coordinate of the point at which the livewire starts</summary>
		uint _x;
		
		/// <summary>The Y coordinate of the point at which the livewire starts</summary>
		uint _y;

		/// <summary>The minimum amount of room to have weights calculated for</summary>
		uint _min_room;

		/// <summary>True if the livewire algorithm has reached a given pixel</summary>
		SparseMatrix<bool> _visited;

		/// <summary>The priority queue representing the points in the edge while running the algorithm</summary>
		PointPriorityQueue _edge;

		/// <summary>The livewire algorithm trace: the point which led to the pixel's best score (in the form x+y*w)</summary>
		SparseMatrix<uint> _trace;

	public:
		/// <summary>Create a new livewire calculator with a private set of weights</summary>
		LivewireCalculator();
		/// <summary>Create a new livewire calculator with an external set of weights - you must delete the weights object when apporpiate</summary>
		LivewireCalculator(Weights* weights);
		virtual ~LivewireCalculator();

		void SetImage(const byte* image, uint w, uint h, Weights::DataFormat format, uint stride);
		void SetSettings(const Weights::Settings& settings);
		const Weights::Settings& GetSettings() const;

		/// <summary>Starts the livewire-calculating thread for the given point, see prepareLivewire for the function that actually does the work</summary>
		/// <param name="x">The X coordinate at which to calculate the livewire data for, it is given a score of 0</param>
		/// <param name="y">The Y coordinate at which to calculate the livewire data for, it is given a score of 0</param>
		/// <param name="y">The minimum amount of room around the point for which weights should be calculating</param>
		void Start(uint x, uint y, uint min_room);

	protected:
		/// <summary>Run the livewire LivewireCalculator</summary>
		void Run();

	private:
		inline void CalcPoint(const uint x, const uint y, const uint i, const bool diagonal, const uint I, const uint S);

	private slots:
		void OnWeightsStopping();
		void OnWeightsSettingsChanged(bool size, bool image);

	public:
		/// <summary>
		/// Gets the livewire trace using all of the precomputed data.
		/// If the livewire data is not completely computed yet the available livewire data is used if the end point is computed.
		/// </summary>
		/// <param name="x">The X coordinate of the end of the trace</param>
		/// <param name="y">The Y coordinate of the end of the trace</param>
		QVector<QPoint> GetTrace(uint x, uint y) const;

#ifndef IMOD_PLUGIN
		/// <summary>
		/// Draws the livewire trace using all of the precomputed data.
		/// If the livewire data is not completely computed yet the available livewire data is used if the end point is computed.
		/// </summary>
		/// <param name="x">The X coordinate of the end of the trace</param>
		/// <param name="y">The Y coordinate of the end of the trace</param>
		/// <param name="painter">The painter to use</param>
		void DrawTrace(uint x, uint y, QPainter &painter) const;
#endif
	};
}

#endif
