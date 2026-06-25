# Makefile for livewire modules

include ../../../configure

O     = $(OBJEXT)
CXXFLAGS = $(LIBCFLAGS)  $(QTFLAGS) -DIMOD_PLUGIN

all : $(OBJS)

$(OBJS) : LivewireCalculator.$(O) moc_LivewireCalculator.$(O) PointPriorityQueue.$(O) Threaded.$(O) \
	 moc_Threaded.$(O) Weights.$(O) moc_Weights.$(O) general.$(O)

moc_Threaded.cpp : Threaded.h
	$(MOC) Threaded.h -o moc_Threaded.cpp

moc_Weights.cpp : Weights.h
	$(MOC) Weights.h -o moc_Weights.cpp

moc_LivewireCalculator.cpp : LivewireCalculator.h
	$(MOC) LivewireCalculator.h -o moc_LivewireCalculator.cpp
