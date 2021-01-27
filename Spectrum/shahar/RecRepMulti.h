int Initialize_cards();
void Close_cards();
int Acquire_scan();
int set_cards_samplerate_andmore(double pixeltime_us, int fsize, int & oversampling);
void prepare_AWG(int binning, int scXsize, int scYsize, double AspectRatio, double rotation, int oversampling, double pixeltime_us,int ScanMode,double flyback_us, int width, int height); 
bool whoiswho (int & cardnoREC, int & cardnoREP, long & lStarHubCarrierIdx);
void retrieve_images(int width, int height, short * pData);
double correctedLOCxy(double LOCxy, double delay, int thesamplerate, int counter, double * prev_correctedLOCxy);
