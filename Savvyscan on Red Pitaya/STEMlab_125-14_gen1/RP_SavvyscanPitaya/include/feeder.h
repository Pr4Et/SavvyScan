// feeder.h
#pragma once
#include <cstdint>
#include <string>

bool monitor_init(); // e.g., "~/serialem_in.txt"
void monitor_shutdown();
void monitor_pump_once();
void checkworking();

void prepare_AWG(int binning, int scXsize, int scYsize, double AspectRatio, double rotation, int oversampling, double pixeltime_us,int ScanMode,double flyback_us, int width, int height,double sampling_pixeltime_us, int fsize);


// If you keep using this name in feeder.cpp:
bool acquire_image(const int width, const int height, short* pData);

// When you later connect hardware, implement this in feeder.cpp.
// For software-only builds, provide a stub that returns 0 (success).
int activate_scan(double integration_time_us, int number_of_triggers);
