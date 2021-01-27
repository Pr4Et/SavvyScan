/*
**************************************************************************
Patern generation, key preparation, hardware settings and acqusition.  
Written by Shahar Seifer 2020
Weizmann institute of Science
Using Spectrum cards. Currently installed:
8ch ADC: M2p.5923-x4  called /dev/spcm0 (card 0)  , 1024MB, 20MS/s
2ch DAC: M2p.6541-x4  called /dev/spcm1 (card 1) , 1024MB, 40MS/s, 16 bits
Star-Hub: Master card1, idx1, connected to hub connector 0
		  Slave  card0, idx0, connected to hub connector 2
Spectrum card commands are based on Spectrum source code
**************************************************************************
*/
 

// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"

#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

#include "RecRepMulti.h"
// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
//added by shahar
#include "..\..\IMOD\include\mrcfiles.h" //add capability to save mrc files
#include <vector>
#include "..\..\MAT\include\mat.h"

/*
**************************************************************************
Global Variables
**************************************************************************
*/
ST_SPCM_CARDINFO    stCardREC;             // info structure of my RECording card (ADC)
ST_SPCM_CARDINFO    stCardREP;             // info structure of my REPlaying card (generator, DAC)
drv_handle			hSync = 0; //handle to STAR-HUB  
drv_handle			handlemaster;
drv_handle			handleslave;
bool                bTimestampInstalled=false;
bool			    bBaseXIOInstalled;
bool                bTrigSrcAvailable = false;
void*               pvBufferREC; //input buffer for ADC
void*               pvBufferREP; //output buffer for DAC
void*               pvBufferKey; //output buffer for DAC
uint64              qwMemInBytesREC;
uint64              qwMemInBytesREP;
//At 5MS/s we read up to 10sec -> 800MByte in 16bits resolution, 8 ch, up to 50M samples per channel, enough for X11 oversampling in 2048X2048)
uint64				MaxNumSamplesPerCh=50*1024*1024; //Num of samples default value, so 50[M]*2[bytes/sample]*8[ch] will not exceed memory on Spectrum AWG card 1024MByte 
uint64				NumSamplesPerCh=MaxNumSamplesPerCh; 
int					MaxSampleRate=MEGA(20); //max capability of ADC card m2p-59
int					samplerate;
int64               llMemsize;
int32               lSegmentsize, lPosttrigger;
uint64              qwTSBufferLen_bytes = 0;
uint64*             pqwTimestamps = NULL;
const int			successful=1;
const double PI = 3.141592653589793;
const int BytesPerSample = 2;
const int REPnumchannels = 2;
const int RECnumchannels = 8;
typedef unsigned short ushort;
//for mat -files
const char* MAT_field_names[] = { "fileversion","tomography_type", "tiltXangle", "tiltYangle", "tiltZangle", "scan_type", "time_resolution_us", "pixelXnm","pixelYnm", "fullXsize", "fullYsize", "netXsize", "netYsize", "netXstart", "netYstart", "samples_per_pixel"  };
const int MAT_numfields = 16;
const int MAT_dim = 1;
mwSize MAT_dims[MAT_dim] = { 1 };
double MAT_fullXsize = 0;
double MAT_fullYsize = 0;
double MAT_netXsize = 0;
double MAT_netYsize = 0;
double MAT_netXstart = 0;
double MAT_netYstart = 0;
double MAT_samples_per_pixel = 0;



//shahar, global vairbales from Shadow GUI to RecRepMulti.cpp and JeolCamera.cpp
 extern char GUI_foldername[239];
 extern int GUI_scanmode;
 extern double GUI_scanargument;
 extern double GUI_aspectratio;
 extern int GUI_InputAmpmV;
 extern int GUI_OutputAmpmV;
 extern int GUI_ch1_avg;
 extern int GUI_ch2_avg;
 extern int GUI_ch3_avg;
 extern int GUI_ch4_avg;
 extern  bool GUI_ch_update_ready;
 extern bool GUI_tomography;
 extern int GUI_tomoIndex;
 extern int GUI_numberOfSlices;
 extern int GUI_chosenCH;//channel in card to be seen in SeiralEM (others will just be saved in files on Shadowscan)
 extern int GUI_cella;
 extern int PYTHON_ask_defocus;
 extern int GUI_Lothar;
 extern double GUI_tiltangle;
 extern bool GUI_saveMAT;
 extern bool GUI_saveMATKEY;
 MrcHeader hdr;
MrcHeader hdr_python;



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int32 lSegmentsize, int32 lPosttrigger, bool ismaster)
{
    int     i;
    int64   llChannelMask;

    // Multiple Recording setup
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;
    if (pstCard->eCardFunction==AnalogIn) //Analog Input
	{	// set memory size, program all input channels to +/-1 V and 50 ohm termination
		bSpcMSetupModeRecStdSingle (pstCard, llChannelMask, llMemsize, lPosttrigger);
		
		for (i=0; i < pstCard->lMaxChannels; i++)
		{	if (pstCard->bM2i || pstCard->bM2p)
				bSpcMSetupInputChannel (pstCard, i, GUI_InputAmpmV, true);
			else if (pstCard->bM3i || pstCard->bM4i)
				bSpcMSetupPathInputCh (pstCard, i, 0, GUI_InputAmpmV, true, true, false,false);
		}
	}
	else if (pstCard->eCardFunction==AnalogOut) //Analog Output
	{	// define single replay mode , set memory size
		bSpcMSetupModeRepStdSingle (pstCard, llChannelMask, llMemsize);

			// set up the timestamp mode to standard if timestamp is installed
		if (bTimestampInstalled && (pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0)
		{
			int32 lAvailTSModes = 0;
			spcm_dwGetParam_i32 (pstCard->hDrv, SPC_TIMESTAMP_AVAILMODES, &lAvailTSModes);

			int32 lTSCmd = SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL;
			if (lAvailTSModes & SPC_TSFEAT_TRGSRC)
				lTSCmd |= SPC_TSFEAT_TRGSRC; // acquire trigger sources
		   bSpcMSetupTimestamp (pstCard, lTSCmd, 0);
		}
	}


}




/*
**************************************************************************
Initialize Spectrum cards 
**************************************************************************
*/

int Initialize_cards ()
{
    char                szBuffer[1024];     // a character buffer for any messages
    int32               lOversampling = 1;
	int					cardnoREC=1; //RECording card (ADC)
	int					cardnoREP=0; //REPlaying card (DAC, generator card)
	long				cardnoHUBmaster=0; //master card connected directly to the STAR-HUB
    bool				bError = false;



	
	// ------------------------------------------------------------------------
    //Fetch who is who in cards on the computer in case of change in default
	bError=whoiswho (cardnoREC, cardnoREP,  cardnoHUBmaster);   //decides the numbers of cardnoREC, cardnoREP, not necessarily as initiated
     	
	
    // init ADC card, get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    if (bSpcMInitCardByIdx (&stCardREC, cardnoREC))    //Opens and initialize the card
    //if (bSpcMInitCardByIdx (&stCardREC, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCardREC, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCardREC, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCardREC, "Error: Could not open card\n", true);
    // check whether we support this card type in the example
    if ((stCardREC.eCardFunction != AnalogIn) )
        return nSpcMErrorMessageStdOut (&stCardREC, "Error: Card function not supported \n", false);

    // if timestamp is installed we set a flag to support this mode in the example also
    bTimestampInstalled =false;// ((stCardREC.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);
    bBaseXIOInstalled = false;//((stCardREC.lFeatureMap & SPCM_FEAT_BASEXIO) != 0);
    if (stCardREC.bM4i || stCardREC.bM2p)
        bTrigSrcAvailable = true;

	// init DAC card, get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    if (bSpcMInitCardByIdx (&stCardREP, cardnoREP))    //Opens and initialize the card
    //if (bSpcMInitCardByIdx (&stCardREC, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCardREP, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCardREP, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCardREP, "Error: Could not open card\n", true);
    // check whether we support this card type in the example
    if ((stCardREP.eCardFunction != AnalogOut) )
        return nSpcMErrorMessageStdOut (&stCardREP, "Error: Card function not supported \n", false);
 
	bSpcMSetupTrigSoftware (&stCardREP, true); //make sure he knows to use software trigger

	//PARK THE BEAM at +3V,+3V
	printf("Praking beam ... ");
	// set master internal clock
	spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_INTPLL);
	//set sample rate in cards 
	spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_SAMPLERATE, 1000000); //DAC sampling rate 
    //the amplitudes are true only with external 50 Ohm termination, otherwise the amplitude doubles
    bSpcMSetupAnalogOutputChannel (&stCardREP, 0, GUI_OutputAmpmV, 0, 0); //3000mV amplitude
    bSpcMSetupAnalogOutputChannel (&stCardREP, 1, GUI_OutputAmpmV, 0, 0);
    spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_TIMEOUT, 5000);
    spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_TIMEOUT, 5000);
	llMemsize=16;
	bSpcMSetupModeRepStdSingle (&stCardREP, CHANNEL0 | CHANNEL1, llMemsize);
	spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE); //trigger when master card starts
	//park signal at 32,767, 32,767 after replay
    spcm_dwSetParam_i32 (stCardREP.hDrv,SPC_CH0_STOPLEVEL,SPCM_STOPLVL_HIGH); //HIGH:+3V , LOW: -3V
    spcm_dwSetParam_i32 (stCardREP.hDrv,SPC_CH1_STOPLEVEL,SPCM_STOPLVL_HIGH);
	short * pData=new short[llMemsize*2]; //dummy buffer to test replay
	for (int i=0; i<llMemsize*2; i++)
	{	
		pData[i]=(short)0;
	}
	qwMemInBytesREP=llMemsize*2*BytesPerSample;
    spcm_dwDefTransfer_i64 (stCardREP.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pData, 0, qwMemInBytesREP); //CHECK parameter !!!!
    spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
	if (spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT) //for single card we write stCardREC.hDrv instead of hSync
    {
		return nSpcMErrorMessageStdOut (&stCardREC, "... Timeout for parking beam\n", false);
    }
	printf("Successful\n");
	
	//STAR-HUB
	hSync = spcm_hOpen ("sync0");  //open STAR-HUB sycnhronized gang of cards
    if (!hSync)
        {
        printf ("no star-hub found.\n");
        bError = true;
        }
    else
        printf ("Found Star-Hub\n");

	//Use STAR-HUB to distribute CLOCK and TRIGGER
	spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << cardnoREP) | (1<< cardnoREC)); //set synchronization to the ADC and DAC cards
	if (cardnoREP==cardnoHUBmaster)
	{	handlemaster=stCardREP.hDrv;
		handleslave=stCardREC.hDrv;
	}
	else if (cardnoREC==cardnoHUBmaster)
	{	handlemaster=stCardREC.hDrv;
		handleslave=stCardREP.hDrv;
	}
	else
        return nSpcMErrorMessageStdOut (&stCardREP, "Error: STAR-HUB not connected to one of the main cards \n", false);


    spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_TIMEOUT, 70000);
    spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_TIMEOUT, 70000);
	//TRIGGER setting: Triggering is needed only in master card, the rest follow the shared line via STAR-HUB (all trigger modes are shared)
	spcm_dwSetParam_i32 (handleslave, SPC_TRIG_ORMASK, SPC_TMASK_NONE);
	spcm_dwSetParam_i32 (handlemaster, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE); //trigger when master card starts
	printf("STEM system is ready for imaging\n");


	return successful;
	
}


/*
**************************************************************************
SET Spectrum cards with samplerate and other parameters 
**************************************************************************
*/

int set_cards_samplerate_andmore(double pixeltime_us, int fsize, int & oversampling)
{
	oversampling=(int)floor((double)MaxNumSamplesPerCh/fsize);
	int oversampling2=(int)floor((double)MaxSampleRate/(1000000.0/pixeltime_us));
	oversampling=oversampling2<oversampling? oversampling2:oversampling; //take miniumum, since two limitations must be met: max samplerate and max number of samples 

	if (oversampling==0 || samplerate>MaxSampleRate)
	{
		printf("ADC Card cannot fit requirements \n");
		return 0; //error
	}
	if (oversampling>10)
	{
		NumSamplesPerCh=10*fsize;
		oversampling=10;
	}
	else
		NumSamplesPerCh=fsize*oversampling;//It must be that fsize*oversampling=NumSamplesPerCh, but here is also the MaxSampleRate
	//size must be in stpes of 8 (otherwise error)
	NumSamplesPerCh=((int)floor((double)NumSamplesPerCh/8.0)+1.0)*8;

	samplerate=(int)(oversampling*1000000.0/pixeltime_us);

	printf("Sample Rate = %d,  oversampling = %d \n",samplerate,oversampling);

	// set master internal clock
	spcm_dwSetParam_i32 (handlemaster, SPC_CLOCKMODE, SPC_CM_INTPLL);
	//set sample rate in cards 
	spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_SAMPLERATE, samplerate); //ADC sampling rate 
	spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_SAMPLERATE, samplerate); //DAC sampling rate 

    char                szBuffer[1024];     // a character buffer for any messages
	
	// do the card setup
    llMemsize =     NumSamplesPerCh; //the number of samples per channel: 50 Millions samples 
    lSegmentsize =  NumSamplesPerCh; //only one timestamp for full block
    lPosttrigger =  NumSamplesPerCh; //all samples after the trigger
	if (!stCardREC.bSetError)
		vDoCardSetup (&stCardREC,  lSegmentsize, lPosttrigger, stCardREC.hDrv==handlemaster);
    if (!stCardREP.bSetError)
 		vDoCardSetup (&stCardREP,  lSegmentsize, lPosttrigger, stCardREP.hDrv==handlemaster);
 
    // REC: calculate the amount of data we need and allocate memory buffer
        qwMemInBytesREC = NumSamplesPerCh * BytesPerSample * RECnumchannels;
        pvBufferREC =(void *) new char[qwMemInBytesREC];// pvAllocMemPageAligned (qwMemInBytesREC);

        if (!pvBufferREC )
            return nSpcMErrorMessageStdOut (&stCardREC, "Memory allocation error\n", false);
    
	// REPcalculate the amount of data we need and allocate memory buffer
        qwMemInBytesREP = NumSamplesPerCh * BytesPerSample * REPnumchannels;
        pvBufferREP = (void *) new char[qwMemInBytesREP];//pvAllocMemPageAligned (qwMemInBytesREP);
        pvBufferKey = (void *) new char[NumSamplesPerCh * sizeof(double) * REPnumchannels];//pvAllocMemPageAligned (qwMemInBytesREP);
        if (!pvBufferREP || !pvBufferKey )
            return nSpcMErrorMessageStdOut (&stCardREP, "Memory allocation error\n", false);
 
 
	if (stCardREC.bSetError)
		return nSpcMErrorMessageStdOut (&stCardREC, "Fault in card settings \n", false);
	else if (stCardREP.bSetError)
		return nSpcMErrorMessageStdOut (&stCardREP, "Fault in preparation \n", false);



		if (bTimestampInstalled) //store timestamp only for DAC card (REP)
		{
			// M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
			if (stCardREP.bM2i || stCardREP.bM3i)
				qwTSBufferLen_bytes = sizeof (uint64) * llMemsize / lSegmentsize;
			else if (stCardREP.bM4i || stCardREP.bM2p)
				qwTSBufferLen_bytes = 2ULL*sizeof (uint64) * llMemsize / lSegmentsize;
            pqwTimestamps =new uint64[2ULL* llMemsize / lSegmentsize]; //(uint64*)pvAllocMemPageAligned (qwTSBufferLen_bytes);
			// if using timestamps we need to start the transfer before the card start to avoid an overrun of the timestamp memory
            if (pqwTimestamps)
			{
				printf ("Starting the timestamp DMA transfer\n");
				// M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
				if (stCardREC.bM2i || stCardREC.bM3i)
					spcm_dwDefTransfer_i64 (stCardREC.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, llMemsize / lSegmentsize * sizeof (uint64));
				else if (stCardREC.bM4i || stCardREC.bM2p)
					spcm_dwDefTransfer_i64 (stCardREC.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, 2 * llMemsize / lSegmentsize * sizeof (uint64));
				spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_M2CMD, M2CMD_EXTRA_STARTDMA);
			}
		}

		if (stCardREC.bSetError || stCardREP.bSetError)
		{
			return nSpcMErrorMessageStdOut (&stCardREP, "Fault in timestamp engagement\n", false);
		}
 
	return successful;
}

/*
**************************************************************************
Prepare AWG (Arbitraty wave generator) data for scanx/y 
**************************************************************************
*/
void prepare_AWG(int binning, int scXsize, int scYsize, double AspectRatio, double rotation, int oversampling, double pixeltime_us,int ScanMode,double flyback_us, int width, int height) 
{
	//pvBufferREP will contain the x-scan (line) and y-scan (frame) data as two channels
	// scXsize, scYsize are the total number of pixels in a scan, scXsize includes also the margins that will be ignored in the image transfered to SerialEM
	//AspectRatio is delY/delX of the pixel size. In Hoppe's scheme we should use AspectRatio=cos(tiltangle)
	//The pixel size in x direction should be universaly constant, with same voltage swing
	//Both Gathan and FEI scan from positive voltage to negative both in x and y, so you see the lines travel from left to right and the bright point of acceleation is on the front side
	//when the beam is parking the voltage is minus,minus and the point is at the rare right corner
	double standard_amplitude=GUI_OutputAmpmV;
	int xscan_millivolts_amp=int(standard_amplitude*(scXsize*binning)/2048.0);
	int yscan_millivolts_amp=int(standard_amplitude* AspectRatio * (scYsize*binning)/2048.0);
	double costheta=cos(rotation*PI/180);
	double sintheta=sin(rotation*PI/180);
	int xscan_millivolts_amp_tag=xscan_millivolts_amp*costheta+yscan_millivolts_amp*sintheta;
	int yscan_millivolts_amp_tag=xscan_millivolts_amp*sintheta+yscan_millivolts_amp*costheta;
	double scanx,scany;//would be scan in millivolt if not for rotation
	int DACx_tag=0;
	int DACy_tag=0;
	double prevLOCx=0;
	double prevLOCy=0;
	double locx;
	double locy;
	double cor_locx=0.0; //was short
	double cor_locy=0.0; //was short
	int marginx,marginy;
	int leftover=0;
	int counter=0;
	int rnbs,phnmbs;
	int scanr;
	double scanph,scanphcor;
	double xdrift;
	//int pixeldelay=(int)(33*((30.0/16.0)/(pixeltime_us/(binning*binning))));
	double delay=0.000220; //FEI deflection coils L/R, in seconds 
	marginx=(short)(0.5*(scXsize-width));
	marginy=(short)(0.5*(scYsize-height));

	MAT_fullXsize = scXsize;
	MAT_fullYsize = scYsize;
	MAT_netXsize = width;
	MAT_netYsize = height;
	MAT_netXstart = marginx;
	MAT_netYstart = marginy;
	MAT_samples_per_pixel = oversampling;

	//the scan output is set in milivolts, 1mV resolution for 50 Ohm load. the value x means a range +/-x[mV]
	//the amplitudes are true only with external 50 Ohm termination, otherwise the amplitude doubles
     bool low_impedance_load=true;//if true then up to 6V alowed, otherwise up to 12V
	 if (xscan_millivolts_amp_tag>6000)
	 {
		 low_impedance_load=false;
		 printf("$$$ WARNING, AMPLITUDE EXCEEDS RECOMEDNED RANGE $$$");
	 }
	 bSpcMSetupAnalogOutputChannel (&stCardREP, 0, xscan_millivolts_amp_tag, 0, 0,16L, !low_impedance_load);
     bSpcMSetupAnalogOutputChannel (&stCardREP, 1, yscan_millivolts_amp_tag, 0, 0,16L, !low_impedance_load);

	//The numbers are stored in short format (signed 16bits)
	short * AWGtable=(short *)pvBufferREP;
	double * Keytable=(double *)pvBufferKey;
	
	//here is scan mode 1. Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
	switch (ScanMode) 
	{
	case 1: //### standard scan: Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
		for (int scanyvar=0;scanyvar<scYsize; scanyvar++)
		{
			//scan from left to right
			scany=yscan_millivolts_amp*(scanyvar-scYsize/2.0)/(scYsize/2.0);
			for (int scanxvar=0;scanxvar<scXsize*oversampling; scanxvar++)
			{

				scanx=xscan_millivolts_amp*(scanxvar-scXsize*oversampling/2.0)/(scXsize*oversampling/2.0);
				DACx_tag=-((int)((scanx*costheta+scany*sintheta)*32768.0/xscan_millivolts_amp_tag));
				DACy_tag=(int)((-scanx*sintheta+scany*costheta)*32768.0/yscan_millivolts_amp_tag);
				if (DACx_tag>32767) DACx_tag=32767;
				if (DACy_tag>32767) DACy_tag=32767;
				if (DACx_tag<-32768) DACx_tag=-32768;
				if (DACy_tag<-32768) DACy_tag=-32768;
				* AWGtable++=(short)DACx_tag;//store scanx_tag value
				* AWGtable++=(short)DACy_tag;//store scany_tag value
				locx=floor((double)scanxvar/oversampling)-marginx;
				locy=scanyvar-marginy;
				cor_locx=correctedLOCxy(locx, delay, samplerate, counter, & prevLOCx);
				cor_locy=correctedLOCxy(locy, delay, samplerate, counter, & prevLOCy);
				counter++;
				if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
				{* Keytable++=cor_locx; //location in the final image (+margin) for the data to fit into
				* Keytable++=cor_locy;
				}
				/*else
				{* Keytable++= -1.0;
				* Keytable++= -1.0;
				}*/
			}
			//flyback
			for (int scanxvar=(int)(flyback_us/pixeltime_us)*oversampling-1; scanxvar>=0 ; scanxvar--)
			{
				scanx=xscan_millivolts_amp*(scanxvar-(int)(flyback_us/pixeltime_us)*oversampling/2)/((int)(flyback_us/pixeltime_us)*oversampling/2);
				DACx_tag=-((int)((scanx*costheta+scany*sintheta)*32768.0/xscan_millivolts_amp_tag));
				DACy_tag=(int)((-scanx*sintheta+scany*costheta)*32768.0/yscan_millivolts_amp_tag);
				if (DACx_tag>32767) DACx_tag=32767;
				if (DACy_tag>32767) DACy_tag=32767;
				if (DACx_tag<-32768) DACx_tag=-32768;
				if (DACy_tag<-32768) DACy_tag=-32768;
				* AWGtable++=(short)DACx_tag;//store scanx_tag value
				* AWGtable++=(short)DACy_tag;//store scany_tag value
				counter++;
				* Keytable++= -1.0;
				* Keytable++= -1.0;
			}
		}
		leftover=llMemsize-scYsize*(scXsize*oversampling+(int)(flyback_us/pixeltime_us)*oversampling);
	break;

	case 2: //### SPIRAL SCAN (from outside in, round, with sure cover)
		rnbs=(int)(scXsize/2);//number of radius points
		counter=0;
		for (int scanrvar=rnbs;scanrvar>=0; scanrvar--) //radius of circle, enclosing a square widthXwidth
		{
			if (scanrvar!=0)
				phnmbs=((int)(2.0*PI*scanrvar))*oversampling; ////number of phase points: so the arc steps are 1pixel/oversampling long
			else
				phnmbs=1*oversampling;;
			for (int scanphvar=0; scanphvar<phnmbs; scanphvar++)
			{
				scanph=scanphvar*2.0*PI/phnmbs;
				//scanphcor=pixeldelay*2*PI/phnmbs;
				DACx_tag=-((int)(sin(scanph+rotation*PI/180.0)*scanrvar*32767.0/(rnbs-1))); //added (-) on 6.8.2020 since stage navigation in x contradicted the image
				DACy_tag=(int)(cos(scanph+rotation*PI/180.0)*scanrvar*32767.0/(rnbs-1));
				locx=((double)scanrvar*sin(scanph+rotation*PI/180.0)+rnbs-marginx);
				locy=((double)scanrvar*cos(scanph+rotation*PI/180.0)+rnbs-marginy);
				cor_locx=correctedLOCxy(locx, delay, samplerate, counter, & prevLOCx);
				cor_locy=correctedLOCxy(locy, delay, samplerate, counter, & prevLOCy);
				counter++;
				if (DACx_tag>32767) DACx_tag=32767;
				if (DACy_tag>32767) DACy_tag=32767;
				if (DACx_tag<-32768) DACx_tag=-32768;
				if (DACy_tag<-32768) DACy_tag=-32768;
				* AWGtable++=(short)DACx_tag;//store scanx_tag value
				* AWGtable++=(short)DACy_tag;//store scany_tag value
				if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
				{
					* Keytable++=cor_locx;
					* Keytable++=cor_locy;
				}
				/*else
				{
					* Keytable++= -1.0;
					* Keytable++= -1.0;
				}*/
			}
		}
		leftover=llMemsize-counter;
		break;


	case 3: //### Linear Mandala (circle motion smoothly drifting along x direction)
		rnbs=(int)(width/2);//number of circling times, so the right hand of the circle interlaces with left hand, jumping 2 pixels per cycle, over 2 width distance
		//rnbs=(int)(width);//number of circling times, case with 2/3 aspect ratio
		scanr=int(height/2); //radius of the circle
		phnmbs=((int)(2.0*PI*scanr))*oversampling; ////number of phase points in one round: so the arc steps are 1pixel/oversampling long
		counter=0;
		for (double scanphvar=0; scanphvar<phnmbs*rnbs; scanphvar++) //repeat not only phases in one round, but the number of round as well, exceeding 2PI many times 
		{
				scanph=scanphvar*2.0*PI/phnmbs;
				scanphcor=scanph - int(scanph/(2.0*PI))*2.0*PI; //modulus 2PI
				xdrift=-width/2+scanph/PI; //drift smoothly one pixel per half a round
				//xdrift=-width+scanph/PI; //drift : case with 2/3 AspectRatio
				DACx_tag=(int)((sin(scanphcor)*scanr+xdrift)*32767.0/(scXsize/2.0));
				DACy_tag=(int)(cos(scanphcor)*scanr*32767.0/(scYsize/2.0));
				locx=((double)scanr*sin(scanphcor)+xdrift+scXsize/2.0-marginx);
				locy=((double)scanr*cos(scanphcor)+scYsize/2.0-marginy);
				cor_locx=correctedLOCxy(locx, delay, samplerate, counter, & prevLOCx);
				cor_locy=correctedLOCxy(locy, delay, samplerate, counter, & prevLOCy);
				counter++;
				if (DACx_tag>32767) DACx_tag=32767;
				if (DACy_tag>32767) DACy_tag=32767;
				if (DACx_tag<-32768) DACx_tag=-32768;
				if (DACy_tag<-32768) DACy_tag=-32768;
				* AWGtable++=(short)DACx_tag;//store scanx_tag value
				* AWGtable++=(short)DACy_tag;//store scany_tag value
				if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
				{
					* Keytable++=cor_locx;
					* Keytable++=cor_locy;
				}
				/*else
				{
					* Keytable++=-1.0;
					* Keytable++=-1.0;
				}*/
		}
		leftover=llMemsize-counter;
		break;


		default:
			printf("Error:  No such ScanMode");

	}
	
	for (int i=0; i<leftover; i++)
	{
			* AWGtable++=(short)DACx_tag;//store last scanx_tag value
			* Keytable++= cor_locx;
			* AWGtable++=(short)DACy_tag;//store last scany_tag value
			* Keytable++= cor_locy;
	}

}

/*
**************************************************************************
Calculate DAC of scanx/y with correction of LR circuit 
**************************************************************************
*/
double correctedLOCxy(double LOCxy, double delay, int thesamplerate, int counter, double * prev_correctedLOCxy)
{
	double result;
	double previous= *prev_correctedLOCxy;
	double dt=1.0/thesamplerate;
	if (counter==0) previous=LOCxy;
	result=(LOCxy+(delay/dt)*previous)/(1+(delay/dt));
	if (result>32766) result=32766;
	if (result<-32767) result=-32767;
	*prev_correctedLOCxy=result;
	return (double)result;
	
}


/*
**************************************************************************
Retrieve and save images 
**************************************************************************
*/
void retrieve_images(const int width, const int height, short * pData)
{
	//split pvBufferREC to 8 channels, deceipher the signals according to pvBufferKey
	int ch,counter;
	short* fillcounter = new short[height * width];
	//int32* fillresults = new int32[height * width];
	double* fillresults = new double[height * width];
	short* oneimage = new short[height * width];
	
	//fillcounter=pvAllocMemPageAligned (height*width*2);
	
	
	//MAT_field_names[] = {"fileversion", "tomography_type", "tiltXangle", "tiltYangle", "tiltZangle", "scan_type", "time_resolution_us", "pixelXnm","pixelYnm", "fullXsize", "fullYsize", "netXsize", "netYsize", "netXstart", "netYstart", "samples_per_pixel" };
	MATFile* pmat;
	mxArray* series;
	mxArray* header;
	mxArray* mx;
	char varname[12];
	header=mxCreateStructArray(MAT_dim, MAT_dims, MAT_numfields, MAT_field_names);
	//maybe field_value = mxCreateDoubleMatrix(1,1,mxREAL);
	if (GUI_saveMAT)
	{
		mx = mxCreateString("2020Sep5");
		mxSetField(header, 0, "fileversion", mx);
		mx = mxCreateString("Single axis tilt tomography");
		mxSetField(header, 0, "tomography_type", mx);
		mx = mxCreateDoubleScalar(GUI_tiltangle);
		mxSetField(header, 0, "tiltXangle", mx);
		mx = mxCreateDoubleScalar(0);
		mxSetField(header, 0, "tiltYangle", mx);
		mx = mxCreateDoubleScalar(0);
		mxSetField(header, 0, "tiltZangle", mx);
		mx = mxCreateString((GUI_scanmode==2)? "savvy-spiral1":((GUI_scanmode == 1) ? "savvy-raster1": "savvy-mandala2"));
		mxSetField(header, 0, "scan_type", mx);
		mx = mxCreateDoubleScalar(1000000.0 / samplerate);
		mxSetField(header, 0, "time_resolution_us", mx);
		mx = mxCreateDoubleScalar((double)GUI_cella * 0.1 / width);
		mxSetField(header, 0, "pixelXnm", mx);
		mx = mxCreateDoubleScalar((double)GUI_cella * 0.1 * GUI_aspectratio / height);//cella is the size in angstroms of net image x size
		mxSetField(header, 0, "pixelYnm", mx);
		mx = mxCreateDoubleScalar(MAT_fullXsize);
		mxSetField(header, 0, "fullXsize", mx);
		mx = mxCreateDoubleScalar(MAT_fullYsize);
		mxSetField(header, 0, "fullYsize", mx);
		mx = mxCreateDoubleScalar(MAT_netXsize);
		mxSetField(header, 0, "netXsize", mx);
		mx = mxCreateDoubleScalar(MAT_netYsize);
		mxSetField(header, 0, "netYsize", mx);
		mx = mxCreateDoubleScalar(MAT_netXstart);
		mxSetField(header, 0, "netXstart", mx);
		mx = mxCreateDoubleScalar(MAT_netYstart);
		mxSetField(header, 0, "netYstart", mx);
		mx = mxCreateDoubleScalar(MAT_samples_per_pixel);
		mxSetField(header, 0, "samples_per_pixel", mx);
	}


	short * ppTChannelData;
	ppTChannelData=new short[llMemsize];
	double* ChannelData;
	ChannelData = new double[llMemsize];
	int32 meanvalue[8];
	double sumval;
	double factorx, factory;
	int intlocx, intlocy, intloc;
	double valsave;
	for (ch=0; ch<8; ch++)
	{
	}
	double locx,locy; //was short
	char shiftno;
	int32 minval,maxval;
	if (GUI_tomoIndex==0)
	{
		mrc_head_new(&hdr, width, height, GUI_numberOfSlices, MRC_MODE_USHORT);	//hdr.mode=MRC_MODE_USHORT  : we say it is ushort but serialEM expect short values
		//Lothar's program extract the x and y pixel sizes from the spans cella.x and y (according to MRC 2014 standard), which is xlen and ylen in MRCHeader
		mrc_set_scale(&hdr, (double)GUI_cella/width, (double)GUI_cella*GUI_aspectratio/height, 0);//double x, double y, double z=0 so will not update, since this is a tilt series. Remember the rotation axis is x so hoppe scheme requires updating the y scale
	}
	//save the files for python program
	mrc_head_new(&hdr_python, width, height, 1, MRC_MODE_USHORT);	//hdr.mode=MRC_MODE_USHORT  : we say it is ushort but serialEM expect short values
	mrc_set_scale(&hdr_python, (double)GUI_cella/width, (double)GUI_cella*GUI_aspectratio/height, 0);//double x, double y, double z=0 so will not update, since this is a tilt series. Remember the rotation axis is x so hoppe scheme requires updating the y scale
	char fname_python[320];
	FILE *fp_python;

	char fname[320];
	FILE *fp;
	int sign_factor=1;
    // split data function (segment 0) to 8 channels
	for (ch=7; ch>=0; ch--)
	{
		//if (ch>=1 && ch<=6) //was until 22Dec23
		if (ch >= 5 && ch <= 7) //turn all to bright field
				sign_factor=-1;
		else
			sign_factor=1;
		bSpcMDemuxAnalogDataOneCH<short> (ch, &stCardREC, pvGetSegmentDataPointer (&stCardREC, pvBufferREC, lSegmentsize, 0, BytesPerSample), lSegmentsize, ppTChannelData);
		short * CHdata=(short *)ppTChannelData; //pointer to search over channel data
		double * LOCdata=(double *)pvBufferKey; //pointer to search over key (location) data
		if (GUI_saveMAT)
		{
			for (int j = 0; j < lSegmentsize; j++)
			{
				ChannelData[j] = (double)ppTChannelData[j];
			}
		}
		for (int j = 0; j <height*width ; j++)
		{
			fillresults[j]=0.0;
			fillcounter[j]=0;
		}
		for (int j = 0; j <lSegmentsize ; j++, CHdata++)
		{
			locx=* LOCdata++;
			locy=* LOCdata++;
			intlocx = (int)(locx+0.5);
			intlocy = (int)(locy+0.5);
			if (intlocx >= 0 && intlocy >= 0 && intlocx < width && intlocy < height)
			{
				intloc = intlocx + intlocy * width;
				fillresults[intloc] += (double)(*CHdata);
				fillcounter[intloc] += 1;
				//valsave = (double)(*CHdata);
				//fillresults[intloc] += valsave;
				//factorx = 1.0 - (locx - (double)intlocx);
				//factory = 1.0 - (locy - (double)intlocy);
			    //fillcounter[intloc] += (factorx < factory) ? factorx : factory;
				/*if (intlocx+1 < width && intlocy+1 <height && 1==0)
				{
					intlocx++;
					intloc++;
					factorx = 1.0 - factorx;
					fillresults[intloc] += valsave;
					fillcounter[intloc] += (factorx < factory) ? factorx : factory;
					intlocy++;
					intloc+=width;
					factory = 1.0 - factory;
					fillresults[intloc] += valsave;
					fillcounter[intloc] += (factorx < factory) ? factorx : factory;
					intlocx--;
					intloc--;
					factorx = 1.0 - factorx;
					fillresults[intloc] += valsave;
					fillcounter[intloc] += (factorx < factory) ? factorx : factory;
				}*/
			}

		}
		counter=0;
		sumval=0;
		 for (int j = 0; j <height*width ; j++)
		{
			if (fillcounter[j]>0)
			   fillresults[j]=sign_factor*(fillresults[j]/fillcounter[j]); //use average (all points marked with same x,y are oversampling points
			counter++;	
			sumval+=fillresults[j];
		}
		meanvalue[ch]=(long)(sumval/counter);
		double ptN = 0, ptS = 0, ptE = 0, ptW = 0;
		for (int j = 0; j < height * width; j++)
		{
			if (fillcounter[j] == 0)
			{
				fillresults[j] = meanvalue[ch];
				if (j > width + 1 && j < (height - 1) * width - 1)
				{
					ptN = fillresults[j - width];
					ptS = fillresults[j + width];
					ptE = fillresults[j + 1];
					ptW = fillresults[j - 1];
					if (ptN > 0 && ptS > 0 && ptE > 0 && ptW > 0)
						fillresults[j] = (ptN + ptS + ptE + ptW) / 4.0;
				}

			}
		}

		for (int j = 0; j <height*width ; j++)
		{
			oneimage[j]=(ushort)((int)(fillresults[j])+(int)32767); //also can write (ushort)result+ushort(32767)
			if (ch==GUI_chosenCH)
			{	
				//image to SerialEM is short: choose the right 16bit settings in serialEM, or choose camera- devide 16 bits by 2 
				pData[j]=(short)fillresults[(j%width)+(height-(j/width)-1)*width]; //copy image and invert y axis since serialEM flips it even in saved mrc files compared to saved mrc file here
			}
		}


		//save fillresults as mrc file
		sprintf(fname, "d:/SavvyscanData/CH%d.mrc", ch);
		if (GUI_tomoIndex==0)
			fp = fopen(fname, "wb");
		else
			fp = fopen(fname, "ab");
		if (fp) {
			if (GUI_tomoIndex==0)
				mrc_head_write(fp, &hdr);
			mrc_write_slice((void*)oneimage, fp, &hdr, GUI_tomoIndex, 'Z'); //save mrc stack as slice=GUI_tomoIndex
			fclose(fp);
		}
		
		/*sprintf(fname, "d:/ShadowData/CH%d.dat", ch);
		fp = fopen(fname, "wb");
		if (fp) {
			fwrite((void*)ppTChannelData, BytesPerSample, llMemsize, fp);
			fclose(fp);
		}*/

		if (GUI_saveMAT)
		{
			sprintf(fname, "d:/SavvyscanData/CH%d.mat", ch);
			if (GUI_tomoIndex == 0)
				pmat = matOpen(fname, "wz"); //wz opens compressed MAT-file
			else
				pmat = matOpen(fname, "u"); //update compressed MAT-file (append slots)
			if (pmat != NULL) {
				sprintf(varname, "header%d", GUI_tomoIndex);
				matPutVariable(pmat, varname, header);
				series = mxCreateDoubleMatrix(1, llMemsize, mxREAL);
				memcpy((void*)(mxGetPr(series)), (void*)ChannelData, llMemsize*sizeof(double));
				sprintf(varname, "series%d", GUI_tomoIndex);
				matPutVariable(pmat, varname, series);
				matClose(pmat);
			}

		}


		//for python program
		if (ch>=1 && ch<=4)
		{
			sprintf(fname_python, "C:/Users/stem/Lothar/CH%d.mrc", ch);
			fp_python = fopen(fname_python, "wb");
			if (fp_python) {
				mrc_head_write(fp_python, &hdr_python);
				mrc_write_slice((void*)oneimage, fp_python, &hdr_python, 0, 'Z'); //save mrc stack as slice=GUI_tomoIndex
				fclose(fp_python);
			}

		}

		}
	printf("ch0=%d, ch1=%d, ch2=%d, ch3=%d, ch4=%d, ch5=%d, ch6=%d, ch7=%d\n",meanvalue[0],meanvalue[1],meanvalue[2],meanvalue[3],meanvalue[4],meanvalue[5],meanvalue[6],meanvalue[7]);
	GUI_ch1_avg=meanvalue[1];
	GUI_ch2_avg=meanvalue[2];
	GUI_ch3_avg=meanvalue[3];
	GUI_ch4_avg=meanvalue[4];
	GUI_ch_update_ready=true;

	if (GUI_saveMATKEY && GUI_saveMAT)
	{
		sprintf(fname, "d:/SavvyscanData/KEY-mode%d-%dx%d.mat", GUI_scanmode, MAT_netXsize, MAT_netYsize);
		pmat = matOpen(fname, "wz"); //wz opens compressed MAT-file
		if (pmat != NULL) {
			sprintf(varname, "header");
			matPutVariable(pmat, varname, header);
			series = mxCreateDoubleMatrix(1, llMemsize, mxREAL);
			double* LOCdata = (double*)pvBufferKey;
			double* LOCtarget = mxGetPr(series);
			for (int j = 0; j < llMemsize; j++)
			{
				*LOCtarget++ = *LOCdata++;
				LOCdata++;
			}
			sprintf(varname, "locationsX");
			matPutVariable(pmat, varname, series);
			series = mxCreateDoubleMatrix(1, llMemsize, mxREAL);
			LOCdata = (double*)pvBufferKey;
			LOCtarget = mxGetPr(series);
			for (int j = 0; j < llMemsize; j++)
			{
				LOCdata++;
				*LOCtarget++ = *LOCdata++;
			}
			sprintf(varname, "locationsY");
			matPutVariable(pmat, varname, series);
			matClose(pmat);
		}

	}




	delete [] ppTChannelData;
	delete [] oneimage;	
	delete [] fillcounter;
	delete [] fillresults;
	
	delete [] pvBufferREC;
    delete [] pvBufferREP;
	delete [] pvBufferKey;
	//RUN Lothar's program to analyse defocus
	if (GUI_Lothar==1 && PYTHON_ask_defocus==0)
		PYTHON_ask_defocus=1; //signal to send request to Python program to check the files
		//system("C:\\Users\\stem\\Lothar\\dist\\DPC-Shift-To-Aberration-wconfig.exe");
		//system("C:\\Users\\stem\\anaconda3\\Scripts\\conda.exe run -n base python C:\\Users\\stem\\Lothar\\DPC-Shift-To-Aberration-wconfig.py");

}


void Close_cards()
{
    // signal will park at zero since system closes off
	spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_RESET) ; 

	// clean up and close the driver
    vSpcMCloseCard (&stCardREC);
    //vFreeMemPageAligned (pvBufferREC, qwMemInBytesREC);
	vSpcMCloseCard (&stCardREP);
    //vFreeMemPageAligned (pvBufferREP, qwMemInBytesREP);
}

/*
**************************************************************************
Perform DAC and ADC 
**************************************************************************
*/

int Acquire_scan()
{
    
	//--------------------------------------------------------------------------
    // Copy the DAC AWG buffer to DAC card memory
        //spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_OUTBUFSIZE, llHWBufSize);
        //spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP);
		//spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, llSWBufSize);
		char                szBuffer[1024];     // a character buffer for any messages
	

		printf ("Starting the DMA buffer transfer to AWG, and waiting until data is in board memory\n");
        spcm_dwDefTransfer_i64 (stCardREP.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBufferREP, 0, qwMemInBytesREP); //CHECK parameter !!!!
        spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        // check for error code
        if (spcm_dwGetErrorInfo_i32 (stCardREP.hDrv, NULL, NULL, szBuffer))
            {
            //vFreeMemPageAligned (pvBufferREP, qwMemInBytesREP);
			delete [] pvBufferREC;
			delete [] pvBufferREP;
			delete [] pvBufferKey;
			return nSpcMErrorMessageStdOut (&stCardREP, szBuffer, false);
            }
        printf ("... data has been transferred to board memory\n");

		//bSpcMSetupTrigSoftware (&stCardREP, true); //make sure he knows to use software trigger  ############################

		//set rule to keep last replayed sample
		spcm_dwSetParam_i32 (stCardREP.hDrv,SPC_CH0_STOPLEVEL,SPCM_STOPLVL_HIGH);
		spcm_dwSetParam_i32 (stCardREP.hDrv,SPC_CH1_STOPLEVEL,SPCM_STOPLVL_HIGH);

		
    // ------------------------------------------------------------------------
    // make acquisition and get data
    // We'll start and wait untill the card has finished or until a timeout occurs

	//spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);

    printf ("Starting all cards and waiting for ready interrupt\n"); 
    if (spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT) //for single card we write stCardREC.hDrv instead of hSync
        {
		delete [] pvBufferREC;
		delete [] pvBufferREP;
		delete [] pvBufferKey;
		//vFreeMemPageAligned (pvBufferREC, qwMemInBytesREC);
        return nSpcMErrorMessageStdOut (&stCardREC, "... Timeout\n", false);
        }
    else
        {
		//spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_STOP ); //try to park the beam at -3V,-3V and not reset the voltage

        // we define the buffer for transfer and start the DMA transfer of the DATA
        printf ("Starting the DMA transfer and waiting until data is in PC memory\n");
        spcm_dwDefTransfer_i64 (stCardREC.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pvBufferREC, 0, qwMemInBytesREC);
        spcm_dwSetParam_i32 (stCardREC.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        // check for error code
        if (spcm_dwGetErrorInfo_i32 (stCardREC.hDrv, NULL, NULL, szBuffer))
            {
            //vFreeMemPageAligned (pvBufferREC, qwMemInBytesREC);
			delete [] pvBufferREC;
			delete [] pvBufferREP;
			delete [] pvBufferKey;
			printf ("... acquisition did not end well ...\n");
			return nSpcMErrorMessageStdOut (&stCardREC, szBuffer, false);
            }
		else
        printf ("... acquisition ended, data has been transferred to PC memory\n");

        // wait for the timestamps (should be already done as we started the transfer before the card start)
        if (bTimestampInstalled)
            {
            spcm_dwSetParam_i32 (stCardREP.hDrv, SPC_M2CMD, M2CMD_EXTRA_WAITDMA);
            printf ("... timestamps have been transferred to PC memory\n");
            }
        }



    // ------------------------------------------------------------------------
    // we go through the segments, split the data in separate channels and show some results
    if (!stCardREC.bSetError)
    {

        // some additional information on the acquisition
        printf ("Each segment is %.3f ms long\n", 1000.0 * lSegmentsize / samplerate);

 
        if (bTimestampInstalled)
		{	int32 lSegmentIdx=0;
            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (stCardREC.bM2i || stCardREC.bM3i)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[lSegmentIdx]) / stCardREC.llSetSamplerate / stCardREC.lOversampling);
            else if (stCardREC.bM4i || stCardREC.bM2p)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[2 * lSegmentIdx]) / stCardREC.llSetSamplerate);
            //vFreeMemPageAligned (pqwTimestamps, qwTSBufferLen_bytes);
			delete pqwTimestamps;
		}
    }
    


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCardREC.bSetError)
        return nSpcMErrorMessageStdOut (&stCardREC, "An error occured while programming the card:\n", true);



    return EXIT_SUCCESS;
 }

	/*
**************************************************************************
szTypeToName: doing name translation
**************************************************************************
*/

char* szTypeToName (int32 lCardType)
    {
    static char szName[50];
    switch (lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     sprintf (szName, "M2i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M2IEXPSERIES:  sprintf (szName, "M2i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M3ISERIES:     sprintf (szName, "M3i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M3IEXPSERIES:  sprintf (szName, "M3i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M4IEXPSERIES:  sprintf (szName, "M4i.%04x-x8", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M4XEXPSERIES:  sprintf (szName, "M4x.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M2PEXPSERIES:  sprintf (szName, "M2p.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
        default:                sprintf (szName, "unknown type");                               break;
        }
    return szName;
    }

/*
**************************************************************************
identify the cards 
**************************************************************************
*/
bool whoiswho (int & cardnoREC, int & cardnoREP, long & lStarHubCarrierIdx)
{
    //if error return false
	drv_handle  ahCard[16];
    int32       lCardCount;
    int32       lCardType, lSerialNumber, lFncType, lFeatures;
    int16*      apnData[16];
    char        szErrorTextBuffer[ERRORTEXTLEN], szName[50];
 
    // ------------------------------------------------------------------------
    // Open all cards and close them just to get who is who information
    for (lCardCount = 0; lCardCount < 8; lCardCount++)
    {
        // uncomment the second line and replace the IP address to use remote
        // cards like in a digitizerNETBOX
        sprintf (szName, "/dev/spcm%d", lCardCount);
        // sprintf (szName, "TCPIP::192.168.1.10::inst%d::INSTR", lCardCount);
        ahCard[lCardCount] = spcm_hOpen (szName);
        apnData[lCardCount] = NULL;

        // not one card found
        if (!lCardCount && !ahCard[lCardCount])
            {
            printf ("no card found...\n");
            return false; //error
            }

        // no more cards found in system
        if (!ahCard[lCardCount])
            break;

        // read out some info and print it
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCITYP,            &lCardType);
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCISERIALNO,       &lSerialNumber);
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_FNCTYPE,           &lFncType);


        // we check which card carries the StarHub
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCIFEATURES,       &lFeatures);
        if (lFeatures & (SPCM_FEAT_STARHUB4 | SPCM_FEAT_STARHUB16))
            lStarHubCarrierIdx = lCardCount;
        
        switch (lFncType)
            {
            case SPCM_TYPE_AI:  
                printf ("ADC Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
                cardnoREC=lCardCount;
                break;
            case SPCM_TYPE_AO:  
                printf ("DAC Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
                cardnoREP=lCardCount;
                break;
            default:
                printf ("Card: %s sn %05d not supported\n", szTypeToName (lCardType), lSerialNumber);            
                break;
            }
            
        spcm_vClose (ahCard[lCardCount]);

   }
	return true; //finished OK
}


