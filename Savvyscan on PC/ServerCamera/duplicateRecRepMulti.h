extern int Initialize_cards();
extern void Close_cards();
extern int Acquire_scan();
extern int set_cards_samplerate_andmore(double pixeltime_us, int fsize, int & oversampling);
extern void prepare_AWG(int binning, int scXsize, int scYsize, double AspectRatio, double rotation, int oversampling, double pixeltime_us,int ScanMode,double flyback_us, int width, int height); 
extern bool whoiswho (int & cardnoREC, int & cardnoREP, long & lStarHubCarrierIdx);
extern void retrieve_images(int width, int height, short * pData);
extern short correctedLOCxy(double LOCxy, double delay, int thesamplerate, int counter, double * prev_correctedLOCxy);