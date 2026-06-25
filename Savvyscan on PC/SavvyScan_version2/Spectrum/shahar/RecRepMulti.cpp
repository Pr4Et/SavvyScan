/*
**************************************************************************
Patern generation, key preparation, hardware settings and acqusition.  
Written by Shahar Seifer 2020
Weizmann institute of Science
Using Spectrum cards. Currently installed:
8ch ADC: M2p.5923-x4  called /dev/spcm0 (card 0)  , 1024MB, 20MS/s
2ch DAC: M2p.6541-x4  called /dev/spcm1 (card 1) , 1024MB, 40MS/s, 16 bits
Star-Hub: 24Apr22 Slave card1, idx1, connected to hub connector 0
		  24Apr22 Master card0, idx0, connected to hub connector 2
Based on Spectrum source code
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
#include <iostream>// calibration file
#include <fstream>//calibration file
#include <string>//calibration file
#include <sstream>//calibration file
#include <string.h>
#include <math.h> 
#include<windows.h>//for sleepand and cmd
//added by shahar
#include "..\..\IMOD\include\mrcfiles.h" //add capability to save mrc files
#include <vector>
#include "..\..\MAT\include\mat.h"
#include<thread> //added 25aug22

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
uint64				MaxNumSamplesPerCh=512*1024*1024; //Max Num of samples per channel that is reasonable to work with 
//uint64				LimitNumSamplesPerCh = 50 * 1024 * 1024; //Num of samples per channel,above which we start to decrease oversampling
//uint64				NumSamplesPerCh=MaxNumSamplesPerCh;
int					MaxSampleRate=MEGA(20); //max capability of ADC card m2p-59
int					samplerate,DACsamplerate;
uint64              NumSamplesPerCh, NumSamplesPerChREC,NumLocPerCh;//the number of samples per channel
int32               lPosttrigger;
uint64              qwTSBufferLen_bytes = 0;
uint64*             pqwTimestamps = NULL;
const int			successful=1;
const double PI = 3.141592653589793;
const int BytesPerSample = 2;
const int REPnumchannels = 2;
const int RECnumchannels = 8;
const int64 MaxllBlockToRec =496*1024*1024/ RECnumchannels; // size of card allows up to 512MB samples - I give a little less for chance of header data is also stored
const int MaxOversampling = 64;
int NotifyBlock = 32 * 1024 * 1024;
int NotifyBlock_slow = 32 * 1024 * 1024;
int NotifyBlock_fast = 8 * 1024 * 1024;
int ScanTimeMax_ms = 10000;
int suggest_pretrigger = 512;// //necessary 13apr2022,  minimum is 8 but will not work for short scans, so use 512. Note that AWG card will get trigger only after this count of pretrigger samples
bool use_star_hub = false; //if false, connect x0 of ADC to clk-in in AWG and x1 of ADC to trg-in in AWG
int REPtrigger_to_output_delay_base = 64;// AWG spec trigger to output delay=63 sample clocks+7ns, is only for external input trigger. Output trigger is immediatley at play start
int trigger_to_output_delay;
int suggest_more_samples = 2* 64 * MaxOversampling; //acquire more samples per channel to compensate for delay between trigger to AWG and actual output  
int saved_bias[8] = {0,0,0,0,0,0,0,0};
//static char *file_calibration= "c:\\Users\\stem\\Savvy\\calibration.csv";
static char file_calibration[36] = "c:\\Users\\stem\\Savvy\\calibration.csv";
static char file_runArina[34] = "c:\\Users\\stem\\Savvy\\run_Arina.bat";
static byte ArinaGate_pattern[5000*8];
static int max_record_count = 5000*8;
static int packet_id_last = 0;
int recent_oversampling = 0;
int pretrigger=0,more_samples=0;
bool alreadyInitalized = false;

typedef unsigned short ushort;
//for mat -files
const char* MAT_field_names[] = { "fileversion","tomography_type", "tiltXangle", "tiltYangle", "tiltZangle", "scan_type", "time_resolution_us", "pixelXnm","pixelYnm", "fullXsize", "fullYsize", "netXsize", "netYsize", "netXstart", "netYstart", "samples_per_pixel", "SerialEM_exposure_us" };
const int MAT_numfields = 17;
const int MAT_dim = 1;
mwSize MAT_dims[MAT_dim] = { 1 };
double MAT_fullXsize = 0;
double MAT_fullYsize = 0;
double MAT_netXsize = 0;
double MAT_netYsize = 0;
double MAT_netXstart = 0;
double MAT_netYstart = 0;
double MAT_samples_per_pixel = 0;
double MAT_exposure_us = 0;
bool   LowRes = false;
bool  repeat_parameters = false;
int prev_binning = 0, prev_scXsize = 0, prev_scYsize = 0, prev_oversampling = 0, prev_ScanMode = 0, prev_width = 0, prev_height = 0;
double prev_AspectRatio = 0, prev_rotation = 0, prev_pixeltime_us = 0, prev_GUI_OutputAmpmV = 0;


//shahar, global vairbales from Shadow GUI to RecRepMulti.cpp and JeolCamera.cpp
 extern char GUI_foldername[239];
 extern int GUI_scanmode;
 extern double GUI_scanargument;
 extern double GUI_aspectratio;
 extern int GUI_InputAmpmV;
 extern int GUI_OutputAmpmV;
 extern int GUI_BiasOutputP;
 extern int GUI_ch1_avg;
 extern int GUI_ch2_avg;
 extern int GUI_ch3_avg;
 extern int GUI_ch4_avg;
 extern int GUI_ch5_avg;
 extern  bool GUI_ch_update_ready;
 extern bool GUI_tomography;
 extern int GUI_tomoIndex;
 extern int GUI_numberOfSlices;
 extern int GUI_chosenCH;//channel in card to be seen in SeiralEM (others will just be saved in files on Shadowscan)
 extern int GUI_cella;
 extern int GUI_LRcella;
 extern int PYTHON_ask_defocus;
 extern int GUI_Lothar;
 extern double GUI_tiltangle;
 extern bool GUI_saveMAT;
 extern bool GUI_saveMATKEY;
 extern bool GUI_align;
 extern int GUI_CalibrateBias;
 extern  bool save4python_latch;
 extern int GUI_SerialEMAligned;
 extern int GUI_LowResTimeS;
 extern bool GUI_ArinaON;
 extern char GUI_ArinaFileName[200] ;
 FILE** fpv=new FILE*[8];
 MrcHeader* hdr=new MrcHeader[8];
 MrcHeader hdr_python;
 MrcHeader hdr_align;
 int LR_last_tomoindex = 0;

 //global valiables for multithreading tasks
 MATFile** pmat=new MATFile*[8];
 mxArray** series=new mxArray*[8];
 mxArray* header;
 mxArray* mx;
 char varname[12];
 FILE* fp_align;
 bool for_alignment = false;
 bool ref_for_alignment = false;
 char fname_python[320];
 FILE* fp_python;
 int32 meanvalue[8];
 bool chwrite_done[8];

/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, bool ismaster)
{
    int     i;
    int64   llChannelMask;
	bool dwError;

    // Multiple Recording setup.  Note: all numbers are in terms of sample count per channel
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;
    if (pstCard->eCardFunction==AnalogIn) //Analog Input
	{	// set memory size, program all input channels to +/-1 V and 50 ohm termination
		//bSpcMSetupModeRecStdSingle (pstCard, llChannelMask, NumSamplesPerCh, lPosttrigger);
		//int64 llBlockToRec = 10 * 1024 * 1024; //chunks of data passage to PC memory between lines in program, must be multiplies of 4k
		int64 llLoopToRec = 1;// NumSamplesPerCh / llBlockToRec;
		int64 llBlockToRec = NumSamplesPerChREC;//25apr2022
		if (llBlockToRec > MaxllBlockToRec)
		{
			llBlockToRec = MaxllBlockToRec;
			llLoopToRec = NumSamplesPerChREC / llBlockToRec;
			if (NumSamplesPerChREC - llLoopToRec * llBlockToRec > 0) llLoopToRec++;
		}
		//if (NumSamplesPerChREC>llLoopToRec*llBlockToRec) llLoopToRec++;
		dwError=bSpcMSetupModeRecFIFOSingle (pstCard, llChannelMask, pretrigger, llBlockToRec, llLoopToRec); //pretrigger of 8 pulses is minimum, but start command is followed seamlessly by a trigger 


		//spcm_dwSetParam_i32(pstCard, SPCM_X0_MODE, SPCM_XMODE_CLKOUT); // X0 set to clock output ADDED on 29Dec21 for implementation of extranal switching buffer/integrator

		for (i=0; i < pstCard->lMaxChannels; i++)
		{	if (pstCard->bM2i || pstCard->bM2p)
				bSpcMSetupInputChannel (pstCard, i, GUI_InputAmpmV, false);//bool bTerm(50 Ohm), int32 lInputOffset, bool bDiffInput
			else if (pstCard->bM3i || pstCard->bM4i)
				bSpcMSetupPathInputCh (pstCard, i, 0, GUI_InputAmpmV, false, false, true,false); //bool bTerm, bool bACCoupling, bool bBWLimit, bool bDiffInput
		}
	}
	else if (pstCard->eCardFunction==AnalogOut) //Analog Output
	{	// define single replay mode , set memory size
		bSpcMSetupModeRepStdSingle (pstCard, llChannelMask, NumLocPerCh);
		
	}

}




/*
**************************************************************************
Initialize Spectrum cards 
**************************************************************************
*/

int Initialize_cards ()
{
   
	if (alreadyInitalized)	return successful; //caused by mutiple savvyscan cameras in serialEM, skip

	
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
 
	//bSpcMSetupTrigSoftware (&stCardREP, true); //make sure he knows to use software trigger on the master

	//TOOL TEST ONLY
	printf("Parking beam ... ");
	// set master internal clock
	spcm_dwSetParam_i64 (stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_INTPLL);
	//set sample rate in cards 
	spcm_dwSetParam_i64 (stCardREP.hDrv, SPC_SAMPLERATE, 1000000); //DAC sampling rate 
    //the amplitudes are true only with external 50 Ohm termination, otherwise the amplitude doubles
    bSpcMSetupAnalogOutputChannel (&stCardREP, 0, GUI_OutputAmpmV, 0, 0); //3000mV amplitude
    bSpcMSetupAnalogOutputChannel (&stCardREP, 1, GUI_OutputAmpmV, 0, 0);
    spcm_dwSetParam_i64 (stCardREP.hDrv, SPC_TIMEOUT, 5000);
    spcm_dwSetParam_i64 (stCardREC.hDrv, SPC_TIMEOUT, 5000);
	NumSamplesPerCh=16;
	bSpcMSetupModeRepStdSingle (&stCardREP, CHANNEL0 | CHANNEL1, NumSamplesPerCh);
	spcm_dwSetParam_i64 (stCardREC.hDrv, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE); //trigger when master card starts
	//park signal at 32,767, 32,767 after replay
    spcm_dwSetParam_i64 (stCardREP.hDrv,SPC_CH0_STOPLEVEL,SPCM_STOPLVL_HIGH); //HIGH:+3V , LOW: -3V
    spcm_dwSetParam_i64 (stCardREP.hDrv,SPC_CH1_STOPLEVEL,SPCM_STOPLVL_HIGH);
	short * pData=new short[NumSamplesPerCh *2]; //dummy buffer to test replay
	for (int i=0; i< NumSamplesPerCh *2; i++)
	{	
		pData[i]=(short)0;
	}
	qwMemInBytesREP= NumSamplesPerCh *2*BytesPerSample;
    spcm_dwDefTransfer_i64 (stCardREP.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pData, 0, qwMemInBytesREP); //CHECK parameter !!!!
    spcm_dwSetParam_i64 (stCardREP.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
	if (spcm_dwSetParam_i64 (stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT) //for single card we write stCardREC.hDrv instead of hSync
    {
		return nSpcMErrorMessageStdOut (&stCardREC, "... Timeout for parking beam\n", false);
    }
	printf("Successful\n");
	
	if (use_star_hub)
	{
		//STAR-HUB
		hSync = spcm_hOpen("sync0");  //open STAR-HUB sycnhronized gang of cards
		if (!hSync)
		{
			printf("no star-hub found.\n");
			//bError = true;
		}
		else
			printf("Found Star-Hub\n");

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
		spcm_dwSetParam_i64 (handleslave, SPC_TRIG_ORMASK, SPC_TMASK_NONE); //was SPC_TMASK_NONE
		spcm_dwSetParam_i64(handlemaster, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE); //trigger when master card starts, should be the REP generator card 
	}
 	//TRIGGER setting: Triggering is needed only in master card, the rest follow the shared line via STAR-HUB (all trigger modes are shared)

	else //if (!use_star_hub)
	{
		//hSync = spcm_hOpen("sync0");
		//spcm_dwSetParam_i32(hSync, SPC_SYNC_ENABLEMASK, 0); //set no synchronization 
		spcm_dwSetParam_i64(stCardREC.hDrv, SPC_CLOCKMODE, SPC_CM_INTPLL);
		//printf("Clock out available\n");
		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_EXTERNAL); // set sample clock input to clock-in socket

		spcm_dwSetParam_i32(stCardREC.hDrv, SPCM_X1_MODE, SPCM_XMODE_TRIGOUT); // X1 set as trigger output, showing when trigger is decided after software request
		spcm_dwSetParam_i32(stCardREC.hDrv, SPCM_X1_MODE, SPCM_XMODE_TRIGOUT); // X1 set as trigger output, showing when trigger is decided after software request
		//ADC X0 always emit clock,  spcm_dwSetParam_i32(stCardREC.hDrv, SPCM_X0_MODE, SPCM_XMODE_CLKOUT); // X0 in ADC is clk output (similar to internal clock)
		printf("Connect trig out from ADC card [x1] to AWG [trig-in] instead of HUB syncing.\n");
		printf("Connect Clk out from ADC card [x0] to AWG [clk-in] for external sync.\n");
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT0); // set sample trigger input to trigger-in socket
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_TRIG_EXT0_MODE, SPC_TM_POS); // Set trig mode to ext. TTL mode (rising edge)
		spcm_dwSetParam_i32(stCardREP.hDrv, SPCM_X1_MODE, SPCM_XMODE_TRIGOUT); // X1 set as trigger output, showing when trigger is decided after software request
		//Instead I set SPC_CLOCKOUT=1 later.  spcm_dwSetParam_i32(stCardREP.hDrv, SPCM_X0_MODE, SPCM_XMODE_CLKOUT); // X0 in AWG clk output (similar to internal clock)
		//To make the ADC reference clock stable for output, according to bjoern.schormann@spec.de: 
		spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP); //should start the clock but not the ADC

		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_INTPLL); // Set to internal clock mode TENTATIVELY
		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_EXTREFCLOCK); // Set to reference clock mode
		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_REFERENCECLOCK, 3200000); // Reference clock that is fed in is 50KHz*64 MHz
		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCK_THRESHOLD, 1500); // set threshold to 1.5V, suitable for 3.3V LVCMOS clock 
		//spcm_dwSetParam_i64(stCardREP.hDrv, SPC_CLOCK50OHM, 1);//set 50 Ohm termination to Clk-in  (0 sets 5K termination)

		spcm_dwSetParam_i64(stCardREC.hDrv, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE);
	}

	spcm_dwSetParam_i64(stCardREP.hDrv, SPC_TIMEOUT, ScanTimeMax_ms);
	spcm_dwSetParam_i64(stCardREC.hDrv, SPC_TIMEOUT, ScanTimeMax_ms);



	std::ifstream myfile(file_calibration);
	std::string line;
	if (myfile.is_open())
	{
		getline(myfile, line);
		std::vector<std::string> split;
		std::stringstream recStream(line);
		std::string field;
		split.clear();
		while (std::getline(recStream, field, ','))
			split.push_back(field);
		for (int ii=0; ii < 8; ii++)
		{
			saved_bias[ii]=std::stoi(split[ii]);
		}
		myfile.close();
	}




	printf("STEM system is ready for imaging\n");
	alreadyInitalized = true;
	return successful;
	
}


/*
**************************************************************************
SET Spectrum cards with samplerate and other parameters 
**************************************************************************
*/

int set_cards_samplerate_andmore(double sampling_pixeltime_us, int fsize, int& oversampling)
{
	if (sampling_pixeltime_us * fsize <= (int)GUI_LowResTimeS * 1000000)
	{
		LowRes = true;
	}
	else
	{
		LowRes = false;
	}
	oversampling = (int)floor((double)MaxNumSamplesPerCh / fsize);
	if (oversampling < 1) oversampling = 1;
	int oversampling2 = (int)floor((double)MaxSampleRate / (1000000.0 / sampling_pixeltime_us));//For 20MHz sampling rate
	oversampling = oversampling2 < oversampling ? oversampling2 : oversampling; //take miniumum, since two limitations must be met: max samplerate and max number of samples 
	if (!LowRes && oversampling > MaxOversampling)
	{
		printf("Pixel time too long, acquisition of 0.1us pulses compromised \n");
		oversampling = MaxOversampling;
	}
	if (LowRes && oversampling > 10)
	{
		NumSamplesPerCh = 10 * fsize;//not for record, we may speed up
		oversampling = 10;
	}
	else
		NumSamplesPerCh = fsize * oversampling;//It must be that fsize*oversampling=NumSamplesPerCh, but here is also the MaxSampleRate
	//size must be in steps of 8 (otherwise error)
	//NumSamplesPerCh=((int)floor((double)NumSamplesPerCh/8.0)+1.0)*8;
	int incomplete8 = NumSamplesPerCh % 8;
	if (incomplete8 > 0) NumSamplesPerCh = NumSamplesPerCh + 8 - incomplete8;

	recent_oversampling = oversampling;
	pretrigger = suggest_pretrigger;// (int)(suggest_pretrigger / oversampling)* oversampling;;
	more_samples = suggest_more_samples;// (int)(suggest_more_samples / oversampling)* oversampling;


	NumLocPerCh = fsize;
	incomplete8 = NumLocPerCh % 8;
	if (incomplete8 > 0) NumLocPerCh = NumLocPerCh + 8 - incomplete8;

	trigger_to_output_delay = REPtrigger_to_output_delay_base * oversampling;

	//>>>>samplerate=(int)(oversampling*1000000.0/ sampling_pixeltime_us);//established value, before roundup
	DACsamplerate = (int)(1000000.0 / sampling_pixeltime_us); //For scanning pattern we only ask one point per pixel
	samplerate = DACsamplerate * oversampling;//ADC clock
	//note the value of samplerate is updates according to true value of DACsamplerate determined by the PLL


	long long lSamplerate;
	if (use_star_hub)
	{
		// set master internal clock
		spcm_dwSetParam_i64(handlemaster, SPC_CLOCKMODE, SPC_CM_INTPLL);
		//set sample rate in cards 
		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_SAMPLERATE, DACsamplerate); //DAC sampling rate 
		//spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKOUT, 1); // enable the clock out
		spcm_dwGetParam_i64(stCardREP.hDrv, SPC_SAMPLERATE, &lSamplerate); // Read back the programmed sample rate
		DACsamplerate = (int)lSamplerate;
		samplerate = (int)(lSamplerate * oversampling);//established value, must be pure integer ratio between clocks
		spcm_dwSetParam_i64(stCardREC.hDrv, SPC_SAMPLERATE, samplerate); //ADC sampling rate 
		spcm_dwGetParam_i64(stCardREC.hDrv, SPC_SAMPLERATE, &lSamplerate); // Read back the programmed sample rate and print
		samplerate = (int)lSamplerate;
	}
	else
	{
		spcm_dwSetParam_i64(stCardREC.hDrv, SPC_CLOCKMODE, SPC_CM_INTPLL);
		spcm_dwSetParam_i32(stCardREC.hDrv, SPC_CLOCKOUT, 1); // enable the clock out
		spcm_dwSetParam_i32(stCardREC.hDrv, SPC_SAMPLERATE, samplerate); //ADC sampling rate 
		spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP); //should start the clock but not the ADC
		Sleep(50);

		//spcm_dwSetParam_i32(stCardREC.hDrv, SPC_CLOCKOUT, 1); // enable the clock out
		spcm_dwGetParam_i64(stCardREC.hDrv, SPC_SAMPLERATE, &lSamplerate); // Read back the programmed ADC sample rate
		if (samplerate != (int)lSamplerate)
		{
			samplerate = (int)lSamplerate;
			printf("####  Sample rate not supported, using approximated !  ####");
		}
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKMODE, SPC_CM_EXTREFCLOCK); // Set AWG to use reference clock mode (will obtain reference from ADC card)
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_REFERENCECLOCK, samplerate); // Reference clock that is fed is samplerate
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_SAMPLERATE, DACsamplerate); // PLL will generate the lower clock frequency according to requested
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_CLOCKOUT, 1); // enable the clock out
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP); //start the clock but not the AWG
		spcm_dwGetParam_i64(stCardREP.hDrv, SPC_SAMPLERATE, &lSamplerate); // Read back the programmed AWG sample rate
		DACsamplerate = (int)lSamplerate;
	}

	printf("Pixel time=%.1f [us],  Sample Rate= %d,  AWG rate= %d,  oversampling = %d \n", 1000000.0/DACsamplerate, samplerate, DACsamplerate,oversampling);

    char                szBuffer[1024];     // a character buffer for any messages
	
	// do the card setup
    //cancelled duplicate variable: llMemsize =     NumSamplesPerCh; //the number of samples per channel: 50 Millions samples 
	NumSamplesPerChREC = NumSamplesPerCh + pretrigger + more_samples;
	if (NumSamplesPerChREC* BytesPerSample * RECnumchannels < NotifyBlock_fast)
	{ 
		NotifyBlock = NotifyBlock_fast;
	}
	else
	{
		NotifyBlock = NotifyBlock_slow;
	}
	
	int NumSampNotifyBlock = (int)(NotifyBlock / (BytesPerSample * RECnumchannels));
	int incomplete_Notify = NumSamplesPerChREC % NumSampNotifyBlock;
	if (incomplete_Notify > 0) NumSamplesPerChREC = NumSamplesPerChREC + NumSampNotifyBlock - incomplete_Notify;
	
	
	if (!stCardREC.bSetError)
		vDoCardSetup (&stCardREC,   stCardREC.hDrv==handlemaster);
    if (!stCardREP.bSetError)
 		vDoCardSetup (&stCardREP,   stCardREP.hDrv==handlemaster);
	

    // REC: calculate the amount of data we need and allocate memory buffer
        qwMemInBytesREC = NumSamplesPerChREC * BytesPerSample * RECnumchannels;
        pvBufferREC =(void *) new char[qwMemInBytesREC ];// pvAllocMemPageAligned (qwMemInBytesREC);
		pvBufferKey = (void*) new char[NumSamplesPerCh * sizeof(double) * REPnumchannels];//will contain the expected locations of the beam in x,y

        if (!pvBufferREC || !pvBufferKey)
            return nSpcMErrorMessageStdOut (&stCardREC, "Memory allocation error\n", false);
    
	// REPcalculate the amount of data we need and allocate memory buffer
        qwMemInBytesREP = NumLocPerCh * BytesPerSample * REPnumchannels; //was NumSamplesPerCh * BytesPerSample * REPnumchannels;
        pvBufferREP = (void *) new char[qwMemInBytesREP];//pvAllocMemPageAligned (qwMemInBytesREP);
        if (!pvBufferREP )
            return nSpcMErrorMessageStdOut (&stCardREP, "Memory allocation error\n", false);
 
 
	if (stCardREC.bSetError)
		return nSpcMErrorMessageStdOut (&stCardREC, "Fault in card settings \n", false);
	else if (stCardREP.bSetError)
		return nSpcMErrorMessageStdOut (&stCardREP, "Fault in preparation \n", false);


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
	//pixeltime_us=exposure, the GUI entered overall time divided by number of pixels in the camera
	if (prev_binning == binning && prev_scXsize == scXsize && prev_scYsize == scYsize && prev_rotation == rotation && prev_oversampling == oversampling && prev_pixeltime_us == pixeltime_us && prev_ScanMode == ScanMode && prev_height == height && prev_width == width && prev_GUI_OutputAmpmV== GUI_OutputAmpmV)
	{
		repeat_parameters = true;
		//continue anyhow since the key is determined here and erased each time.  
	}
	else
	{
		repeat_parameters = false;
		prev_binning = binning;
		prev_scXsize = scXsize;
		prev_scYsize = scYsize;
		prev_rotation = rotation;
		prev_oversampling = oversampling;
		prev_pixeltime_us = pixeltime_us;
		prev_ScanMode = ScanMode;
		prev_height = height;
		prev_width = width;
		prev_GUI_OutputAmpmV = GUI_OutputAmpmV;
	}
	
	double standard_amplitude=GUI_OutputAmpmV;
	int xscan_millivolts_amp=int(standard_amplitude*(scXsize*binning)/8192.0);
	int yscan_millivolts_amp=int(standard_amplitude* AspectRatio * (scYsize*binning)/8192.0);
	double costheta = 1.0;// cos(rotation * PI / 180); Rotation in properties file caused problem
	double sintheta = 0;// sin(rotation * PI / 180);
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
	int counterLoc = 0;
	int rnbs,phnmbs;
	int scanr;
	double scanph,scanphcor;
	double xdrift;
	int epoch = 0;
	double Tcx = 0;
	double Tcy = 0;
	//int pixeldelay=(int)(33*((30.0/16.0)/(pixeltime_us/(binning*binning))));
	double delay=0.000220; //FEI deflection coils L/R, in seconds 
	marginx=(int)(0.5*(scXsize-width));
	marginy=(int)(0.5*(scYsize-height));
	
	MAT_fullXsize = scXsize;
	MAT_fullYsize = scYsize;
	MAT_netXsize = width;
	MAT_netYsize = height;
	MAT_netXstart = marginx;
	MAT_netYstart = marginy;
	MAT_samples_per_pixel = oversampling;
	MAT_exposure_us = pixeltime_us;


	//the scan output is set in milivolts, 1mV resolution for 50 Ohm load. the value x means a range +/-x[mV]
	//the amplitudes are true only with external 50 Ohm termination, otherwise the amplitude doubles
     bool low_impedance_load=true;//if true then up to 6V alowed, otherwise up to 12V
	 if (xscan_millivolts_amp_tag>6000)
	 {
		 low_impedance_load=false;
		 printf("$$$ WARNING, AMPLITUDE EXCEEDS RECOMEDNED RANGE $$$");
	 }

	 int BiasOutput_mv = 0; //ignore GUI_BiasOutputP
	 if (ScanMode == 1) BiasOutput_mv = (int)((double)standard_amplitude * 0.414*(1.0+ (double)GUI_BiasOutputP/100.0) );
	 if (ScanMode == 5) BiasOutput_mv = (int)((double)standard_amplitude * 0.125 * (1.0 + (double)GUI_BiasOutputP / 100.0));
	 bSpcMSetupAnalogOutputChannel (&stCardREP, 0, xscan_millivolts_amp_tag, BiasOutput_mv, 0,16L, !low_impedance_load);
     bSpcMSetupAnalogOutputChannel (&stCardREP, 1, yscan_millivolts_amp_tag, 0, 0,16L, !low_impedance_load);

	//The numbers are stored in short format (signed 16bits)
	short * AWGtable=(short *)pvBufferREP;
	double * Keytable=(double *)pvBufferKey;
	
	//here is scan mode 1. Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
	switch (ScanMode)
	{
	case 1: //### standard scan: Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
		for (int scanyvar = scYsize-1; scanyvar >=0; scanyvar--) //order reversed 8May2022
		{
			//scan from left to right
			scany = (double)yscan_millivolts_amp * ((double)scanyvar - (double)scYsize / 2.0) / ((double)scYsize / 2.0);
			for (int scanxvar = 0; scanxvar < scXsize * oversampling; scanxvar++)
			{

				//scanx = (double)xscan_millivolts_amp * ((double)scanxvar - (double)(scXsize * oversampling) / 2.0) / ((double)(scXsize * oversampling) / 2.0);
				scanx = (double)xscan_millivolts_amp * ((double)scanxvar - 0.5*(double)scXsize * (double)oversampling) / ((double)(scXsize * oversampling) / 2.0);
				DACx_tag = (int)(-(scanx * costheta + scany * sintheta) * 32768.0 / (double)xscan_millivolts_amp_tag);
				DACy_tag = (int)((-scanx * sintheta + scany * costheta) * 32768.0 / (double)yscan_millivolts_amp_tag);
				if (DACx_tag > 32767) DACx_tag = 32767;
				if (DACy_tag > 32767) DACy_tag = 32767;
				if (DACx_tag < -32768) DACx_tag = -32768;
				if (DACy_tag < -32768) DACy_tag = -32768;
				if (counter % oversampling == 0)
				{
					*AWGtable++ = (short)DACx_tag;//store scanx_tag value
					*AWGtable++ = (short)DACy_tag;//store scany_tag value
					counterLoc++;
				}
				locx = (double)(scanxvar / oversampling - 2 * marginx);
				locy = (double)(scanyvar - marginy);
				//cor_locx=correctedLOCxy(locx, delay, samplerate, counter, & prevLOCx);  before 24Feb22
				//cor_locy=correctedLOCxy(locy, delay, samplerate, counter, & prevLOCy);
				cor_locx = checkLOCxy(locx);
				cor_locy = checkLOCxy(locy);
				counter++;
				*Keytable++ = cor_locx; //location in the final image (+margin) for the data to fit into
				*Keytable++ = cor_locy;
			}
			//flyback
			int flypix =(int)(flyback_us / pixeltime_us);
			for (int scanxvar = flypix * oversampling - 1; scanxvar >= 0; scanxvar--)
			{
				//scanx = xscan_millivolts_amp * (scanxvar - (int)(flyback_us / pixeltime_us) * oversampling / 2) / ((int)(flyback_us / pixeltime_us) * oversampling / 2);
				scanx = xscan_millivolts_amp* (double(width + 2 * marginx) * (double)scanxvar / (double)flypix + ((-0.5 * (double)width - (double)marginx) * (double)oversampling)) / ((double)(scXsize * oversampling) / 2.0);
				DACx_tag = (int)(-(scanx * costheta + scany * sintheta) * 32768.0 / (double)xscan_millivolts_amp_tag);
				DACy_tag = (int)((-scanx * sintheta + scany * costheta) * 32768.0 / (double)yscan_millivolts_amp_tag);
				if (DACx_tag > 32767) DACx_tag = 32767;
				if (DACy_tag > 32767) DACy_tag = 32767;
				if (DACx_tag < -32768) DACx_tag = -32768;
				if (DACy_tag < -32768) DACy_tag = -32768;
				if (counter % oversampling == 0)
				{
					*AWGtable++ = (short)DACx_tag;//store scanx_tag value
					*AWGtable++ = (short)DACy_tag;//store scany_tag value
					counterLoc++;
				}
				counter++;
				*Keytable++ = -1.0;
				*Keytable++ = -1.0;
			}
		}
		// = NumSamplesPerCh - scYsize * (scXsize * oversampling + (int)(flyback_us / pixeltime_us) * oversampling);
		break;

	case 2: //### SPIRAL SCAN (from outside in, round, with sure cover)
		rnbs = (int)(scXsize / 2);//number of radius points
		counter = 0;
		for (int scanrvar = rnbs; scanrvar >= 0; scanrvar--) //radius of circle, enclosing a square widthXwidth
		{
			if (scanrvar != 0)
				phnmbs = ((int)(2.0 * PI * scanrvar)) * oversampling; ////number of phase points: so the arc steps are 1pixel/oversampling long
			else
				phnmbs = 1 * oversampling;;
			for (int scanphvar = 0; scanphvar < phnmbs; scanphvar++)
			{
				scanph = scanphvar * 2.0 * PI / phnmbs;
				//scanphcor=pixeldelay*2*PI/phnmbs;
				DACx_tag = -((int)(sin(scanph + rotation * PI / 180.0) * scanrvar * 32767.0 / (rnbs - 1))); //added (-) on 6.8.2020 since stage navigation in x contradicted the image
				DACy_tag = (int)(cos(scanph + rotation * PI / 180.0) * scanrvar * 32767.0 / (rnbs - 1));
				locx = ((double)scanrvar * sin(scanph + rotation * PI / 180.0) + rnbs - marginx);
				locy = ((double)scanrvar * cos(scanph + rotation * PI / 180.0) + rnbs - marginy);
				cor_locx = correctedLOCxy(locx, delay, samplerate, counter, &prevLOCx);
				cor_locy = correctedLOCxy(locy, delay, samplerate, counter, &prevLOCy);
				if (DACx_tag > 32767) DACx_tag = 32767;
				if (DACy_tag > 32767) DACy_tag = 32767;
				if (DACx_tag < -32768) DACx_tag = -32768;
				if (DACy_tag < -32768) DACy_tag = -32768;
				if (counter % oversampling == 0)
				{
					*AWGtable++ = (short)DACx_tag;//store scanx_tag value
					*AWGtable++ = (short)DACy_tag;//store scany_tag value
					counterLoc++;
				}
				counter++;
				if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
				{
					*Keytable++ = cor_locx;
					*Keytable++ = cor_locy;
				}
				/*else
				{
					* Keytable++= -1.0;
					* Keytable++= -1.0;
				}*/
			}
		}
		break;


	case 3: //### Linear Mandala (circle motion smoothly drifting along x direction)
		rnbs = (int)(width);//number of circling times, so the right hand of the circle interlaces with left hand, jumping 2 pixels per cycle, over 2 width distance
		scanr = (int)(height / 2); //radius of the circle
		phnmbs = ((int)(2.0 * PI * (double)scanr)) * oversampling; ////number of phase points in one round: so the arc steps are 1pixel/oversampling long
		counter = 0;
		for (double scanphvar = 0; scanphvar < (double)phnmbs * (double)rnbs; scanphvar++) //repeat not only phases in one round, but the number of round as well, exceeding 2PI many times 
		{
			scanph = scanphvar * 2.0 * PI / (double)phnmbs;
			scanphcor = scanph - (double)((int)(scanph / (2.0 * PI))) * 2.0 * PI; //modulus 2PI
			xdrift = -(double)width + scanph / PI; //drift smoothly one pixel per half a round
			DACx_tag = (int)((sin(scanphcor) * scanr + xdrift) * 32767.0 / ((double)scXsize / 2.0));
			DACy_tag = (int)(cos(scanphcor) * scanr * 32767.0 / ((double)scYsize / 2.0));
			locx = ((double)scanr * sin(scanphcor) + xdrift + (double)scXsize / 2.0 - (double)marginx);
			locy = ((double)scanr * cos(scanphcor) + (double)scYsize / 2.0 - (double)marginy);
			cor_locx = correctedLOCxy(locx, delay, samplerate, counter, &prevLOCx);
			cor_locy = correctedLOCxy(locy, delay, samplerate, counter, &prevLOCy);
			if (DACx_tag > 32767) DACx_tag = 32767;
			if (DACy_tag > 32767) DACy_tag = 32767;
			if (DACx_tag < -32768) DACx_tag = -32768;
			if (DACy_tag < -32768) DACy_tag = -32768;
			if (counter % oversampling == 0)
			{
				*AWGtable++ = (short)DACx_tag;//store scanx_tag value
				*AWGtable++ = (short)DACy_tag;//store scany_tag value
				counterLoc++;
			}
			counter++;
			if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
			{
				*Keytable++ = cor_locx;
				*Keytable++ = cor_locy;
			}
			/*else
			{
				* Keytable++=-1.0;
				* Keytable++=-1.0;
			}*/
		}
		break;

	case 4: //### Lissajous
		// x,y near frequencies with detuning
		rnbs = (int)(width);
		scanr = (int)(width / 2);
		Tcx = 3.0 * rnbs * oversampling; ////number of points in one round: so the arc steps are about 1pixel/oversampling long
		Tcy = (3.0 * rnbs * oversampling) / (1 + 0.5 / rnbs);
		epoch = 3 * rnbs * rnbs * oversampling;
		counter = 0;
		for (int np = 0; np < epoch; np++)
		{
			DACx_tag = (int)(sin(2 * PI * np / Tcx) * scanr * 32767.0 / ((double)scXsize / 2.0));
			DACy_tag = (int)(sin(2 * PI * np / Tcy) * scanr * 32767.0 / ((double)scYsize / 2.0));
			locx = ((double)scanr * sin(2 * PI * np / Tcx) + (double)scXsize / 2.0 - (double)marginx);
			locy = ((double)scanr * sin(2 * PI * np / Tcy) + (double)scYsize / 2.0 - (double)marginy);
			cor_locx = correctedLOCxy(locx, delay, samplerate, counter, &prevLOCx);
			cor_locy = correctedLOCxy(locy, delay, samplerate, counter, &prevLOCy);
			if (DACx_tag > 32767) DACx_tag = 32767;
			if (DACy_tag > 32767) DACy_tag = 32767;
			if (DACx_tag < -32768) DACx_tag = -32768;
			if (DACy_tag < -32768) DACy_tag = -32768;
			if (counter % oversampling == 0)
			{
				*AWGtable++ = (short)DACx_tag;//store scanx_tag value
				*AWGtable++ = (short)DACy_tag;//store scany_tag value
				counterLoc++;
			}
			counter++;
			if (true)//(cor_locx>=0 && cor_locx<width && cor_locy>=0 && cor_locy<height)
			{
				*Keytable++ = cor_locx;
				*Keytable++ = cor_locy;
			}
			/*else
			{
				* Keytable++=-1.0;
				* Keytable++=-1.0;
			}*/
		}
		break;

	case 5: //### 4X4 raster scans:     (in each raster doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left)
		for (double NY = 0; NY < 4; NY++)
		{
			for (double NX = 0; NX < 4; NX++)
			{
				for (int scanyvar = 0; scanyvar < height/4; scanyvar++)
				{
					//scan from left to right
					scany = (double)yscan_millivolts_amp * ((double)scanyvar + ((-0.5+0.25* (double)NY)*(double)scYsize) ) / ((double)scYsize / 2.0);
					for (int scanxvar = 0; scanxvar < (width/4 +2*marginx) * oversampling; scanxvar++)
					{

						scanx = (double)xscan_millivolts_amp * ((double)scanxvar + (((-0.5+0.25*(double)NX)*(double)width -(double)marginx ) * (double)oversampling) ) / ((double)(scXsize * oversampling) / 2.0);
						DACx_tag = (int)(-(scanx * costheta + scany * sintheta) * 32768.0 / (double)xscan_millivolts_amp_tag);
						DACy_tag = (int)((-scanx * sintheta + scany * costheta) * 32768.0 / (double)yscan_millivolts_amp_tag);
						if (DACx_tag > 32767) DACx_tag = 32767;
						if (DACy_tag > 32767) DACy_tag = 32767;
						if (DACx_tag < -32768) DACx_tag = -32768;
						if (DACy_tag < -32768) DACy_tag = -32768;
						if (counter % oversampling == 0)
						{
							*AWGtable++ = (short)DACx_tag;//store scanx_tag value
							*AWGtable++ = (short)DACy_tag;//store scany_tag value
							counterLoc++;
						}
						if (scanxvar >= (2 * marginx) * oversampling)
						{
							locx = (double)(scanxvar / oversampling - 2 * marginx + NX * (width / 4));
							locy = (double)(scanyvar - marginy + NY * (height / 4));
							cor_locx = checkLOCxy(locx);
							cor_locy = checkLOCxy(locy);
						}
						else
						{
							cor_locx = -1;
							cor_locy = -1;
						}
						counter++;
						*Keytable++ = cor_locx; //location in the final image (+margin) for the data to fit into
						*Keytable++ = cor_locy;
					}
					//flyback
					int flypix = (int)(flyback_us / pixeltime_us);
					for (int scanxvar = flypix * oversampling - 1; scanxvar >= 0; scanxvar--)
					{
						//scanx = xscan_millivolts_amp * (scanxvar - (int)(flyback_us / pixeltime_us) * oversampling / 2) / ((int)(flyback_us / pixeltime_us) * oversampling / 2);
						//scanx = xscan_millivolts_amp * (double(width/4+2*marginx)*(double)scanxvar / (double)flypix + (((-0.5 + 0.25 * (double)NX) * (double)width ) * (double)oversampling) ) / ((double)(scXsize * oversampling) / 2.0);
						scanx = xscan_millivolts_amp * (double(width / 4 + 2 * marginx) * (double)scanxvar / (double)flypix + (((-0.5 + 0.25 * (double)NX) * (double)width - (double)marginx) * (double)oversampling)) / ((double)(scXsize * oversampling) / 2.0);
						DACx_tag = -((int)((scanx * costheta + scany * sintheta) * 32768.0 / xscan_millivolts_amp_tag));
						DACy_tag = (int)((-scanx * sintheta + scany * costheta) * 32768.0 / yscan_millivolts_amp_tag);
						if (DACx_tag > 32767) DACx_tag = 32767;
						if (DACy_tag > 32767) DACy_tag = 32767;
						if (DACx_tag < -32768) DACx_tag = -32768;
						if (DACy_tag < -32768) DACy_tag = -32768;
						if (counter % oversampling == 0)
						{
							*AWGtable++ = (short)DACx_tag;//store scanx_tag value
							*AWGtable++ = (short)DACy_tag;//store scany_tag value
							counterLoc++;
						}
						counter++;
						*Keytable++ = -1.0;
						*Keytable++ = -1.0;
					}
				}
			}
		}
		break;

	case 6: //### SPIRAL SCAN FOR DIFFRACTION DISK SHIFT ALIGNMENT
		rnbs = (int)(scXsize / 2);//number of radius points
		counter = 0;
		for (int scanrvar = rnbs; scanrvar >= 0; scanrvar--) //radius of circle, enclosing a square widthXwidth
		{
			if (scanrvar != 0)
				phnmbs = ((int)(2.0 * PI * scanrvar)) * oversampling; ////number of phase points: so the arc steps are 1pixel/oversampling long
			else
				phnmbs = 1 * oversampling;;
			for (int scanphvar = 0; scanphvar < phnmbs; scanphvar++)
			{
				scanph = scanphvar * 2.0 * PI / phnmbs;
				DACx_tag = -((int)(sin(scanph + rotation * PI / 180.0) * scanrvar * 32767.0 / (rnbs - 1))); //added (-) on 6.8.2020 since stage navigation in x contradicted the image
				DACy_tag = (int)(cos(scanph + rotation * PI / 180.0) * scanrvar * 32767.0 / (rnbs - 1));
				cor_locx = ((double)scanrvar * sin(scanph + rotation * PI / 180.0) + rnbs - marginx);
				cor_locy = ((double)scanrvar * cos(scanph + rotation * PI / 180.0) + rnbs - marginy);
				if (DACx_tag > 32767) DACx_tag = 32767;
				if (DACy_tag > 32767) DACy_tag = 32767;
				if (DACx_tag < -32768) DACx_tag = -32768;
				if (DACy_tag < -32768) DACy_tag = -32768;
				if (counter % oversampling == 0)
				{
					*AWGtable++ = (short)DACx_tag;//store scanx_tag value
					*AWGtable++ = (short)DACy_tag;//store scany_tag value
					counterLoc++;
				}
				counter++;
				*Keytable++ = cor_locx;
				*Keytable++ = cor_locy;
				
			}
		}
		break;


		default:
			printf("Error:  No such ScanMode");

	}
	
	leftover = NumSamplesPerCh - counter;
	for (int i=0; i<leftover; i++)
	{
			* Keytable++= cor_locx;
			* Keytable++= cor_locy;
	}
	leftover = NumLocPerCh - counterLoc;
	for (int i = 0; i < leftover; i++)
	{
		*AWGtable++ = (short)32766;//store last scanx_tag value
		*AWGtable++ = (short)32766;//store last scany_tag value
	}


	//determine ArinaGate pattern table
	if ( !LowRes) //leaving out !repeat_parameters &&, to send each time a single pattern so ArinaGate will only transmit pulses when I specifically ask for
	{
		int intlocx, intlocy;
		double* LOCdata = (double*)pvBufferKey; //pointer to search over key (location) data
		int record_position = 0;
		int state_counter = 0;//ArinaGate connected to AWG trig-out, otherwise= REPtrigger_to_output_delay_base* recent_oversampling; // after the trigger-in signal the output is still delayed by REPtrigger_to_output_delay_base AWG clocks
		bool record_state = false;
		bool previous_record_state = false;
		int packet_id = 0;
		int half_oversampling = (int)(recent_oversampling/2);
		for (int j = 0; j < NumSamplesPerCh; j++)
		{
			locx = *LOCdata++;
			locy = *LOCdata++;
			if (j % recent_oversampling == half_oversampling)
			{
				intlocx = (int)(locx);
				intlocy = (int)(locy);
				record_state = (intlocx >= 0 && intlocy >= 0 && intlocx < width && intlocy < height);
				if (record_state == previous_record_state && state_counter < 32765 )
				{
					state_counter++;
				}
				else
				{
					if (record_position <= max_record_count - 4)
					{
						if (record_position % 5000 == 4998) //control stamp at the end of each 5000 byte packet
						{
							if (record_position < max_record_count - 5000)
							{
								ArinaGate_pattern[record_position] = packet_id; //stamp current packet id
								packet_id++;
								record_position++;
								ArinaGate_pattern[record_position] = packet_id; //stamp end packet id that is above current id to signal that another packet should be loaded 
								record_position++;
							}
						}
						//cancelled:state_counter = (int)(state_counter / recent_oversampling); //waiting in terms of AWG clock ticks
						ArinaGate_pattern[record_position] = (byte)(state_counter & 0b11111111);
						record_position++;
						if (previous_record_state)
						{
							ArinaGate_pattern[record_position] = (byte)((state_counter >> 8 & 0b01111111) | 0b10000000);
						}
						else
						{
							ArinaGate_pattern[record_position] = (byte)((state_counter >> 8 & 0b01111111));
						}
						record_position++;
					}
					state_counter = 1;
				}
				previous_record_state = record_state;
			}
		}
		for (int pos = record_position; pos < 4998 + 5000 * packet_id; pos++)
		{
			ArinaGate_pattern[pos] = 0;
		}
		ArinaGate_pattern[4998 + 5000 * packet_id] = packet_id; //Stamp the end of packet series
		ArinaGate_pattern[4999 + 5000 * packet_id] = packet_id;
		packet_id_last = packet_id; //say how many packets to transmit
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
double checkLOCxy(double result)
{
	if (result > 32766.0) result = 32766.0;
	if (result < -32767.0) result = -32767.0;
	return result;

}


/*
**************************************************************************
Retrieve and save images 
**************************************************************************
*/
void retrieve_images(const int width, const int height, short * pData)
{
	//split pvBufferREC to 8 channels, deceipher the signals according to pvBufferKey
	int ch;
	for_alignment=false;
	ref_for_alignment = false;
	if (LowRes && GUI_tomography && GUI_SerialEMAligned==0 && (GUI_tomoIndex == LR_last_tomoindex || GUI_tomoIndex==0))
	{
		ref_for_alignment = true;
	}
	if (LowRes && GUI_tomoIndex >= 0 && GUI_tomography)
	{
		for_alignment = true;
		LR_last_tomoindex = GUI_tomoIndex;
	}

	for (ch = 0; ch < 8; ch++ )
	{
		chwrite_done[ch] = false;
	}
	char fname[320];

	//MAT_field_names[] = {"fileversion", "tomography_type", "tiltXangle", "tiltYangle", "tiltZangle", "scan_type", "time_resolution_us", "pixelXnm","pixelYnm", "fullXsize", "fullYsize", "netXsize", "netYsize", "netXstart", "netYstart", "samples_per_pixel" };
	header=mxCreateStructArray(MAT_dim, MAT_dims, MAT_numfields, MAT_field_names);
	//maybe field_value = mxCreateDoubleMatrix(1,1,mxREAL);
	if (GUI_saveMAT && !for_alignment && GUI_tomoIndex >= 0)
	{
		mx = mxCreateString("2022apr26");
		mxSetField(header, 0, "fileversion", mx);
		mx = mxCreateString("Single axis tilt tomography");
		mxSetField(header, 0, "tomography_type", mx);
		mx = mxCreateDoubleScalar(GUI_tiltangle);
		mxSetField(header, 0, "tiltXangle", mx);
		mx = mxCreateDoubleScalar(0);
		mxSetField(header, 0, "tiltYangle", mx);
		mx = mxCreateDoubleScalar(0);
		mxSetField(header, 0, "tiltZangle", mx);
		mx = mxCreateString((GUI_scanmode == 5) ? "savvy-4x4" : ((GUI_scanmode == 4) ? "savvy-lisajous1" :((GUI_scanmode==2)? "savvy-spiral1":((GUI_scanmode == 1) ? "savvy-raster4": "savvy-mandala2"))));
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
		mx = mxCreateDoubleScalar(MAT_exposure_us);
		mxSetField(header, 0, "SerialEM_exposure_us", mx);
	}


	//save the files for python program
	if (save4python_latch)
	{
		mrc_head_new(&hdr_python, width, height, 1, MRC_MODE_USHORT);	//hdr.mode=MRC_MODE_USHORT  : we say it is ushort but serialEM expect short values
		mrc_set_scale(&hdr_python, (double)GUI_cella / width, (double)GUI_cella * GUI_aspectratio / height, 0);//double x, double y, double z=0 so will not update, since this is a tilt series. Remember the rotation axis is x so hoppe scheme requires updating the y scale
	}
	
	mrc_head_new(&hdr_align, width, height, 1, MRC_MODE_USHORT);	//hdr.mode=MRC_MODE_USHORT  : we say it is ushort but serialEM expect short values

    // ####   Work on 8 channels in parallel  #####
	std::thread th7(process_channel_parallel, 7, width, height, pData);
	std::thread th6(process_channel_parallel, 6, width, height, pData);
	std::thread th5(process_channel_parallel, 5, width, height, pData);
	std::thread th4(process_channel_parallel, 4, width, height, pData);
	std::thread th3(process_channel_parallel, 3, width, height, pData);
	std::thread th2(process_channel_parallel, 2, width, height, pData);
	std::thread th1(process_channel_parallel, 1, width, height, pData);
	std::thread th0(process_channel_parallel, 0, width, height, pData);
	//wait to finish
	int countsum = 0;
	do
	{
		countsum = 0;
		for (ch = 7; ch >= 0; ch--)
		{
			if (chwrite_done[ch]) countsum++;
		}
	} while (countsum < 8);
	//finish formally so error will not occur on sub exit
	th7.join();
	th6.join();
	th5.join();
	th4.join();
	th3.join();
	th2.join();
	th1.join();
	th0.join();

	printf("ch0=%d, ch1=%d, ch2=%d, ch3=%d, ch4=%d, ch5=%d, ch6=%d, ch7=%d (including bias)\n",meanvalue[0],meanvalue[1],meanvalue[2],meanvalue[3],meanvalue[4],meanvalue[5],meanvalue[6],meanvalue[7]);
	//For calibration of ADC reading with zero illumination (with no beam)
	if (GUI_CalibrateBias==1)
	{
		GUI_CalibrateBias=0;
		std::ofstream savefile(file_calibration);
		if (savefile.is_open())
		{
			for (int ii = 0; ii < 8; ii++)
			{
				if (ii < 7)
				{
					saved_bias[ii] = (int)meanvalue[ii];
					savefile << saved_bias[ii] << ",";
				}
				else //ch7 is HAADF so make it always positive
				{
					savefile << "-32768\n";
				}
			}
			savefile.close();
		}

	}
	
	if (for_alignment)
		GUI_align = true; //mark to communicate with shadow GUI
	GUI_ch1_avg=(int)meanvalue[1]- saved_bias[1];
	GUI_ch2_avg= (int)meanvalue[2] - saved_bias[2];
	GUI_ch3_avg= (int)meanvalue[3] - saved_bias[3];
	GUI_ch4_avg= (int)meanvalue[4] - saved_bias[4];
	GUI_ch5_avg = (int)meanvalue[5] - saved_bias[5];
	GUI_ch_update_ready=true;

	if (GUI_saveMATKEY && GUI_saveMAT)
	{
		sprintf(fname, "A:/SavvyscanData/KEY-mode%d-%dx%d.mat", GUI_scanmode, MAT_netXsize, MAT_netYsize);
		pmat[0] = matOpen(fname, "wz"); //wz opens compressed MAT-file
		if (pmat[0] != NULL) {
			sprintf(varname, "header");
			matPutVariable(pmat[0], varname, header);
			series[0] = mxCreateDoubleMatrix(1, NumSamplesPerCh, mxREAL);
			double* LOCdata = (double*)pvBufferKey;
			double* LOCtarget = mxGetPr(series[0]);
			for (int j = 0; j < NumSamplesPerCh; j++)
			{
				*LOCtarget++ = *LOCdata++;
				LOCdata++;
			}
			sprintf(varname, "locationsX");
			matPutVariable(pmat[0], varname, series[0]);
			series[0] = mxCreateDoubleMatrix(1, NumSamplesPerCh, mxREAL);
			LOCdata = (double*)pvBufferKey;
			LOCtarget = mxGetPr(series[0]);
			for (int j = 0; j < NumSamplesPerCh; j++)
			{
				LOCdata++;
				*LOCtarget++ = *LOCdata++;
			}
			sprintf(varname, "locationsY");
			matPutVariable(pmat[0], varname, series[0]);
			matClose(pmat[0]);
		}

	}

	delete [] pvBufferREC;
    delete [] pvBufferREP;
	delete [] pvBufferKey;
	//RUN Lothar's program to analyse defocus
	if (GUI_Lothar==1 && PYTHON_ask_defocus==0 && !for_alignment)
		PYTHON_ask_defocus=1; //signal to send request to Python program to check the files
		//system("C:\\Users\\stem\\Lothar\\dist\\DPC-Shift-To-Aberration-wconfig.exe");
		//system("C:\\Users\\stem\\anaconda3\\Scripts\\conda.exe run -n base python C:\\Users\\stem\\Lothar\\DPC-Shift-To-Aberration-wconfig.py");
}


//sub module to run in multithreading: process and write channels
void process_channel_parallel(int ch, const int width, const int height, short* pData)
{
	short* fillcounter = new short[height * width];
	//int32* fillresults = new int32[height * width];
	double* fillresults = new double[height * width];
	short* oneimage = new short[height * width];
	short* ppTChannelData;
	ppTChannelData = new short[NumSamplesPerChREC];
	double* ChannelData;
	ChannelData = new double[NumSamplesPerCh];
	double sumval;
	double factorx, factory;
	int intlocx, intlocy, intloc;
	double valsave;
	double locx, locy; //was short
	char shiftno;
	int32 minval, maxval;
	char multi_fname[320];

	//if (ch>=1 && ch<=6) //was until 2022Dec23
	//if (false && ch >= 5 && ch <= 7) //turn all to bright field
	//		sign_factor=-1;
	//else
	bSpcMDemuxAnalogDataOneCH<short>(ch, &stCardREC, pvGetSegmentDataPointer(&stCardREC, pvBufferREC, NumSamplesPerChREC, 0, BytesPerSample), NumSamplesPerChREC, ppTChannelData);
	short* CHdata = (short*)(ppTChannelData + pretrigger + trigger_to_output_delay); //careful: shift (26apr22). pointer to search over channel data
	double* LOCdata = (double*)pvBufferKey; //pointer to search over key (location) data
	if (GUI_saveMAT && !for_alignment) //to save time, done only if mat file is saved
	{
		for (int j = 0; j < NumSamplesPerCh-(pretrigger + trigger_to_output_delay); j++)
		{
			ChannelData[j] = (double)(ppTChannelData[j + pretrigger + trigger_to_output_delay]); //careful: shift (26apr22). 
		}
	}
	//26apr22: note that the output of replay is delayed by (pretrigger) number of sample clocks +63 cycles compared to start of ADC acquisition

	for (int j = 0; j < height * width; j++)
	{
		fillresults[j] = 0.0;
		fillcounter[j] = 0;
	}
	//determine filling points, both for pixel calculation and for building gate pattern
	int record_position = 0;
	int state_counter = REPtrigger_to_output_delay_base*recent_oversampling; // after the active trigger signal the output is still delayed by REPtrigger_to_output_delay_base AWG clocks
	bool record_state = false;
	for (int j = 0; j < NumSamplesPerCh; j++, CHdata++)
	{
		locx = *LOCdata++;
		locy = *LOCdata++;
		intlocx = (int)(locx );
		intlocy = (int)(locy);
		if (intlocx >= 0 && intlocy >= 0 && intlocx < width && intlocy < height)
		{
			intloc = intlocx + intlocy * width;
			fillresults[intloc] += (double)(*CHdata);
			fillcounter[intloc] += 1;
		}

	}
	int counter = 0;
	sumval = 0;
	double result_in_pixel, truncated_result_in_pixel;
	for (int j = 0; j < height * width; j++)
	{
		if (fillcounter[j] > 0)
		{
			result_in_pixel = fillresults[j] / (double)fillcounter[j]; //use average (all points marked with same x,y are oversampling points
			truncated_result_in_pixel = result_in_pixel - (double)saved_bias[ch]; //remove the measured bias in zero illumination
			if (truncated_result_in_pixel < 0) truncated_result_in_pixel = 0;
			fillresults[j] = truncated_result_in_pixel;
			counter++;
			sumval += result_in_pixel;
		}
	}
	meanvalue[ch] = (long)(sumval / (double)counter);
	double ptN = 0, ptS = 0, ptE = 0, ptW = 0;
	for (int j = 0; j < height * width; j++)
	{
		if (fillcounter[j] == 0)
		{
			truncated_result_in_pixel = (double)meanvalue[ch] - (double)saved_bias[ch];
			if (truncated_result_in_pixel < 0) truncated_result_in_pixel = 0;
			fillresults[j] = truncated_result_in_pixel;
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

	for (int j = 0; j < height * width; j++)
	{
		oneimage[j] = (ushort)((int)(fillresults[j])); //was fillresults[j])+(int)32767 before I removed the bias measured in zero illumination
		if (ch == GUI_chosenCH)
		{
			//image to SerialEM is short: choose the right 16bit settings in serialEM, or choose camera- devide 16 bits by 2 
			pData[j] = (short)fillresults[(j % width) + (height - (j / width) - 1) * width]; //copy image and invert y axis since serialEM flips it even in saved mrc files compared to saved mrc file here
		}
	}

	//save fillresults as mrc file
	if (GUI_tomoIndex == 0 && !for_alignment)
	{
		mrc_head_new(&hdr[ch], width, height, GUI_numberOfSlices, MRC_MODE_USHORT);	//hdr.mode=MRC_MODE_USHORT  : we say it is ushort but serialEM expect short values
		mrc_set_scale(&hdr[ch], (double)GUI_cella / width, (double)GUI_cella * GUI_aspectratio / height, 0);//double x, double y, double z=0 so will not update, since this is a tilt series. Remember the rotation axis is x so hoppe scheme requires updating the y scale
	}
	if (GUI_tomoIndex == 0) fpv[ch] = NULL;
	if (!for_alignment)
	{
		sprintf(multi_fname, "A:/SavvyscanData/CH%d.mrc", ch);
		if (GUI_tomoIndex == 0 )
			fpv[ch] = fopen(multi_fname, "wb");//write replace
		else //if (GUI_tomoIndex > 0 && !GUI_tomography)
			fpv[ch] = fopen(multi_fname, "ab"); // write append. //XXXXsould not happen, I remove since closeing and opening consumes time
		if (true || fpv[ch] != NULL) {
			if (GUI_tomoIndex == 0)
				mrc_head_write(fpv[ch], &hdr[ch]);
			mrc_write_slice((void*)oneimage, fpv[ch], &hdr[ch], GUI_tomoIndex, 'Z'); //save mrc stack as slice=GUI_tomoIndex
			//if (!GUI_tomography)  fclose(fpv[ch]);
			fclose(fpv[ch]);
		}
	}

	if (ch == GUI_chosenCH && ref_for_alignment && GUI_tomography && !GUI_SerialEMAligned)
	{
		if (GUI_tiltangle >= -3)
		{
			sprintf(multi_fname, "A:/SavvyscanData/lastHR_positive.mrc", ch);
			fp_align = fopen(multi_fname, "wb");
			mrc_head_write(fp_align, &hdr_align);
			mrc_write_slice((void*)oneimage, fp_align, &hdr_align, 0, 'Z');
			fclose(fp_align);
		}
		if (GUI_tiltangle <= 3)
		{
			sprintf(multi_fname, "A:/SavvyscanData/lastHR_negative.mrc", ch);
			fp_align = fopen(multi_fname, "wb");
			mrc_head_write(fp_align, &hdr_align);
			mrc_write_slice((void*)oneimage, fp_align, &hdr_align, 0, 'Z');
			fclose(fp_align);
		}
	}

	if (ch == GUI_chosenCH && for_alignment && !GUI_SerialEMAligned)
	{
		sprintf(multi_fname, "A:/SavvyscanData/newLR4alignment.mrc", ch);
		fp_align = fopen(multi_fname, "wb");
		mrc_head_write(fp_align, &hdr_align);
		mrc_write_slice((void*)oneimage, fp_align, &hdr_align, 0, 'Z');
		fclose(fp_align);
	}

	/*sprintf(multi_fname, "d:/ShadowData/CH%d.dat", ch);
	fp = fopen(multi_fname, "wb");
	if (fp) {
		fwrite((void*)ppTChannelData, BytesPerSample, NumSamplesPerCh, fp);
		fclose(fp);
	}*/

	if (GUI_saveMAT && !for_alignment && GUI_tomoIndex >= 0)
	{
		sprintf(multi_fname, "A:/SavvyscanData/CH%d.mat", ch);
		if (GUI_tomoIndex == 0)
			pmat[ch] = matOpen(multi_fname, "wz"); //wz opens compressed MAT-file
		else
			pmat[ch] = matOpen(multi_fname, "u"); //update compressed MAT-file (append slots)
		if (pmat[ch] != NULL) {
			sprintf(varname, "header%d", GUI_tomoIndex);
			matPutVariable(pmat[ch], varname, header);
			series[ch] = mxCreateDoubleMatrix(1, NumSamplesPerCh, mxREAL);
			memcpy((void*)(mxGetPr(series[ch])), (void*)ChannelData, NumSamplesPerCh * sizeof(double));
			sprintf(varname, "series%d", GUI_tomoIndex);
			matPutVariable(pmat[ch], varname, series[ch]);
			matClose(pmat[ch]);
		}

	}

	//for python program
	if (ch >= 1 && ch <= 4 && save4python_latch)
	{
		sprintf(fname_python, "C:/Users/stem/Lothar/CH%d.mrc", ch);
		fp_python = fopen(fname_python, "wb");
		if (fp_python) {
			mrc_head_write(fp_python, &hdr_python);
			mrc_write_slice((void*)oneimage, fp_python, &hdr_python, 0, 'Z'); //save mrc stack as slice=GUI_tomoIndex
			fclose(fp_python);
		}

	}
	delete[] oneimage;
	delete[] fillcounter;
	delete[] fillresults;
	delete[] ppTChannelData;
	delete[] ChannelData;
	chwrite_done[ch] = true;

}



/* close cards */
void Close_cards()
{
	// signal will park at zero since system closes off
	spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_RESET);

	// clean up and close the driver
	vSpcMCloseCard(&stCardREC);
	//vFreeMemPageAligned (pvBufferREC, qwMemInBytesREC);
	vSpcMCloseCard(&stCardREP);
	//vFreeMemPageAligned (pvBufferREP, qwMemInBytesREP);
}


/*
**************************************************************************
Perform DAC and ADC 
**************************************************************************
*/

int Acquire_scan()
{
	// Send skip patterns to MCU for Arina trigger
	if (!LowRes)
	{
		if (GUI_ArinaON) //Arm Arina detector before each scan
		{
			// use MAT_netXsize = width , MAT_netYsize = height, true pixel time: 1/DACsamplerate
			int integration_time_us = (int)(1000000.0 / DACsamplerate) - 4;
			int number_of_triggers = MAT_netXsize * MAT_netYsize;
			char* namestr;
			namestr = GUI_ArinaFileName;
			if (integration_time_us >= 10)
			{
				
				char line2[250];
				sprintf(line2,"python.exe -m acquireimages -i 192.168.100.70 -t %.6f -n %d -e 200 -o %s_s%d -d c://ArinaData -x",0.000001*(double)integration_time_us, number_of_triggers, GUI_ArinaFileName,GUI_tomoIndex);


				std::ofstream savefile(file_runArina);
				if (savefile.is_open())
				{
					savefile << "c:" << "\n";
					savefile << "cd C:\\Users\\stem\\Desktop\\Velan\\python" <<"\n";
					savefile << line2 <<"\n";
					savefile.close();
					WinExec(file_runArina, SW_SHOWNORMAL);
				}

			}
			else
			{
				printf(">>>   Could not arm Arina detector, use longer pixel time  <<<");
			}

		}

		Send_ArinaGate_pattern();
		Sleep(1000);
	}


	//--------------------------------------------------------------------------
	// Copy the DAC AWG buffer to DAC card memory

	char                szBuffer[1024];     // a character buffer for any messages
	bool dwError=false;

	//if (use_star_hub)
	dwError=spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_INVALIDATEDATA); //necessary 13apr2022
	//spcm_dwSetParam_i64(stCardREC.hDrv, SPC_M2CMD, M2CMD_ALL_STOP);

	if (!repeat_parameters)
	{
		printf("Starting the DMA buffer transfer to AWG, and waiting until data is in board memory\n");
		spcm_dwDefTransfer_i64(stCardREP.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBufferREP, 0, qwMemInBytesREP); //CHECK parameter !!!!
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);///?

		// check for error code
		if (spcm_dwGetErrorInfo_i32(stCardREP.hDrv, NULL, NULL, szBuffer))
		{
			//vFreeMemPageAligned (pvBufferREP, qwMemInBytesREP);
			printf("error in preparation");
			delete[] pvBufferREC;
			delete[] pvBufferREP;
			delete[] pvBufferKey;
			return nSpcMErrorMessageStdOut(&stCardREP, szBuffer, false);
		}
		printf("... data has been transferred to board memory\n");

		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_M2CMD, M2CMD_DATA_STOPDMA); ///? //added to make sure interrupts not comming from AWG card anymore


		//set rule to keep last replayed sample
		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_CH0_STOPLEVEL, SPCM_STOPLVL_HIGH);///?
		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_CH1_STOPLEVEL, SPCM_STOPLVL_HIGH);
	}
	else
	{
		printf("Repeating the AWG pattern\n");
		//dwError=spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_START );//ready to use (since I stoped in previous loop)
		//if (dwError == true)
		//{
		//	printf("Error in AWG card");
		//	return nSpcMErrorMessageStdOut(&stCardREP, szBuffer, false);
		//}
		//set rule to keep last replayed sample
		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_CH0_STOPLEVEL, SPCM_STOPLVL_HIGH);///?
		spcm_dwSetParam_i64(stCardREP.hDrv, SPC_CH1_STOPLEVEL, SPCM_STOPLVL_HIGH);

	}

	//spcm_dwSetParam_i32(stCardREP.hDrv, SPCM_X0_MODE, SPCM_XMODE_CLKOUT); // X0 set to clock output, added 27 Nov 2022
	//spcm_dwSetParam_i32(stCardREP.hDrv, SPCM_X1_MODE, SPCM_XMODE_TRIGOUT); // X1 set to trigout output, True during acquisition/operation, added 29 Nov 2022
	//spcm_dwSetParam_i64(stCardREP.hDrv, SPC_SAMPLERATE, DACsamplerate); //DAC sampling rate //>>
	//spcm_dwSetParam_i64(stCardREC.hDrv, SPC_SAMPLERATE, samplerate); //ADC sampling rate  //>>

	// ------------------------------------------------------------------------
	// make acquisition and get data
	// We'll start and wait untill the card has finished or until a timeout occurs
	int64 llAvailBytes, llBytesPos, lStatus, memcount = 0;
	//NotifyBlock = 32 * 1024 * 1024;
	if (qwMemInBytesREC < NotifyBlock) NotifyBlock = qwMemInBytesREC;
	int64 llLoopToRec = qwMemInBytesREC / NotifyBlock;
	if (qwMemInBytesREC > llLoopToRec * NotifyBlock) llLoopToRec++;

	//try some reset
	//dwError = spcm_dwSetParam_i64(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_INTERNALRESET );

		printf ("Setting DMA transfer to PC memory\n");
    	spcm_dwDefTransfer_i64(stCardREC.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, NotifyBlock, pvBufferREC, 0, qwMemInBytesREC);
		//spcm_dwGetParam_i64(stCardREC.hDrv, SPC_DATA_AVAIL_USER_LEN, &llAvailBytes);
		//spcm_dwGetParam_i64(stCardREC.hDrv, SPC_DATA_AVAIL_USER_POS, &llBytesPos);
		//spcm_dwGetParam_i64(stCardREC.hDrv, SPC_M2STATUS, &lStatus);

		//printf ("Starting all ADC cards and waiting for ready interrupt\n"); 
		printf("Starting all cards: one trigger for all, signal recording in FIFO mode\n");
		printf("Target record size [bytes]: %lld >>> ", qwMemInBytesREC);

		//dwError = spcm_dwSetParam_i32(hSync, SPC_M2CMD, M2CMD_CARD_DISABLETRIGGER);
		spcm_dwSetParam_i64(stCardREC.hDrv, SPC_DATA_AVAIL_CARD_LEN, 0); //added 21June22

		// start DMA
		dwError = spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA);//M2CMD_CARD_START
		Sleep(1);
		//dwError=spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) ; //for single card we write stCardREC.hDrv instead of hSync
		if (use_star_hub)
		{
			dwError = spcm_dwSetParam_i32(hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);
			if (dwError)
			{
				printf("Sync problem (change exposure time).\n");
				return 1;
			}
			//Sleep(1);
			//dwError = spcm_dwSetParam_i32(hSync, SPC_M2CMD, M2CMD_CARD_START |M2CMD_CARD_FORCETRIGGER);
		}
		else
		{
			dwError = spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);//will start only after REC card finished all pretrigger acquisition
			//Sleep(1);
			dwError = spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_FORCETRIGGER);//start pretrigger acquisition followed by trigger pulse sent to AWG card
		}
		

		//starts to count () pretrigger clocks and all post-trigger clocks so to collect total recorded length
		//ForceTrigger: like software trigger
		//Sleep(20);//2 ms delay
		//dwError = spcm_dwSetParam_i32(hSync, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER);

		// we wait for the next data to be available. After this call we get at least Blocksize of data to proceed
		//if (!dwError) dwError = spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_DATA_WAITDMA);

		//FIFO loop according to page 77 in m2p_59xx_manual from Spectra
		for (int loopn=0; loopn<llLoopToRec; loopn++)
		{
			if (!dwError)
			{
				//spcm_dwGetParam_i64 (stCardREC.hDrv, SPC_DATA_AVAIL_USER_POS, &llBytesPos);
				//printf ("We now have %lld new bytes available\n", llAvailBytes);
				//printf ("The available data starts at position %lld\n", llBytesPos);
				// we take care not to go across the end of the buffer
				//if ((llBytesPos + llAvailBytes) >= qwMemInBytesREC)
				//llAvailBytes = qwMemInBytesREC - llBytesPos;
				// our do function gets a pointer to the start of the available data section and the length
				//vDoSomething (&pcData[llBytesPos], llAvailBytes);
				// the buffer section is now immediately set available for the card
				// we wait for the next data to be available / remaining data
				dwError = spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_DATA_WAITDMA);
				// if there was no error we can proceed and read out the available bytes that are free again
				spcm_dwGetParam_i64(stCardREC.hDrv, SPC_DATA_AVAIL_USER_LEN, &llAvailBytes);
				//spcm_dwGetParam_i64(stCardREC.hDrv, SPC_M2STATUS, &lStatus);
				if (llAvailBytes != NotifyBlock)
				{
					printf("?");
				}
				if (llAvailBytes > NotifyBlock)
				{
					llAvailBytes = NotifyBlock;
				}
				if (qwMemInBytesREC - memcount < NotifyBlock)
				{
					llAvailBytes = qwMemInBytesREC - memcount;
				}
				spcm_dwSetParam_i64 (stCardREC.hDrv, SPC_DATA_AVAIL_CARD_LEN, llAvailBytes);
				memcount = memcount + llAvailBytes;
				//printf("mem=%lld out-of %lld,  acquired=%lld \n", memcount, qwMemInBytesREC, llAvailBytes);
				printf(" %lld.", memcount);
			}
			else
			{
				printf ("\n... acquisition did not end well ...\n");
				break;
			}
		}

		//spcm_dwSetParam_i64(stCardREC.hDrv, M2CMD_DATA_SGFLUSH, 0); //(added on 13apr2022)


																			  //Stop to avoid filling the PC memory "by intertia",  M2CMD_CARD_DISABLETRIGGER M2CMD_CARD_INVALIDATEDATA M2CMD_CARD_INTERNALRESET
		if (use_star_hub)
		{
			dwError = spcm_dwSetParam_i32(hSync, SPC_M2CMD, M2CMD_CARD_STOP);///? M2CMD_ALL_STOP M2CMD_CARD_DISABLETRIGGER  22June22: NECESSARY, OTHERWISE MEMORY CONTINUE TO UPDATE AND ERASE
		}
		//dwError = spcm_dwSetParam_i64(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_FLUSHFIFO);///?
		//spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_WAITREADY);//25apr22
		spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);//25apr22
		spcm_dwSetParam_i32(stCardREP.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);//25apr22

		dwError = spcm_dwSetParam_i32(stCardREC.hDrv, SPC_M2CMD, M2CMD_DATA_STOPDMA);
		//test spcm_dwInvalidateBuf(stCardREC.hDrv, SPCM_BUF_DATA);

        // check for error code
        if (spcm_dwGetErrorInfo_i32 (stCardREC.hDrv, NULL, NULL, szBuffer))
            {
            //vFreeMemPageAligned (pvBufferREC, qwMemInBytesREC);
			printf("\n error indication \n");
			delete [] pvBufferREC;
			delete [] pvBufferREP;
			delete [] pvBufferKey;
			return nSpcMErrorMessageStdOut (&stCardREC, szBuffer, false);
            }
		else
        printf ("\n... acquisition ended, data has been transferred to PC memory\n");

    


    // ------------------------------------------------------------------------
    // we go through the segments, split the data in separate channels and show some results
    if (!stCardREC.bSetError)
    {

        // some additional information on the acquisition
        printf ("Each segment is %.3f ms long\n", 1000.0 * NumSamplesPerCh / samplerate);

 
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


/*
**************************************************************************
Perform DAC and ADC
**************************************************************************
*/

bool Send_ArinaGate_pattern() 
{
	HANDLE hComm;
	DCB dcb;
	bool fSuccess;
	COMMTIMEOUTS timeouts;

	hComm = CreateFileA("COM3",                //port name
		GENERIC_READ | GENERIC_WRITE, //Read/Write  
		0,                            // 1=Yes Sharing
		NULL,                         // No Security
		OPEN_EXISTING,// Open existing port only
		0,            // Non Overlapped I/O
		NULL);        // Null for Comm Devices

	if (hComm == INVALID_HANDLE_VALUE)
	{
		printf("Error in opening serial port");
		return false;
	}
	else
	printf("Connection with ArinaGate COM Port \n");



	fSuccess = GetCommState(hComm, &dcb);
	dcb.BaudRate = CBR_38400;     //  baud rate
	dcb.ByteSize = DATABITS_8;             //  data size, xmit and rcv
	dcb.Parity = NOPARITY;      //  parity bit
	dcb.StopBits = ONESTOPBIT;    //  stop bit
	dcb.fDtrControl = 0;
	dcb.fRtsControl = 0;


	fSuccess=SetCommState(hComm, &dcb);
	//PurgeComm(hComm, PURGE_RXABORT |PURGE_RXCLEAR);

	DWORD dNoOfBytesWritten = 0;
	char message[6];

	fSuccess = GetCommTimeouts(hComm, &timeouts);
	/* Set timeout to 0 to force that:
	   If a character is in the buffer, the character is read,
	   If no character is in the buffer, the function do not wait and returns immediatly
	*/
	//timeouts.ReadIntervalTimeout = MAXDWORD;
	//timeouts.ReadTotalTimeoutMultiplier = 0;
	timeouts.ReadTotalTimeoutConstant = 300;
	timeouts.WriteTotalTimeoutConstant = 0;
	timeouts.ReadIntervalTimeout = MAXDWORD;
	timeouts.ReadTotalTimeoutMultiplier = MAXDWORD;
	timeouts.WriteTotalTimeoutMultiplier = MAXDWORD;
	fSuccess = SetCommTimeouts(hComm, &timeouts);

	Sleep(200);
	fSuccess = PurgeComm(hComm, PURGE_TXABORT |
		PURGE_RXABORT |
		PURGE_TXCLEAR |
		PURGE_RXCLEAR);

	bool Status=true;
	int ccount = 0;
	int packet_counter = 0;

	while (Status && packet_counter <= packet_id_last)
	{
		Status = WriteFile(hComm,        // Handle to the Serial port
			(byte*) (ArinaGate_pattern+5000* packet_counter),     // Data to be written to the port   LPCVOID
			5000,  //No of bytes to write
			&dNoOfBytesWritten, //Bytes written
			NULL);

		if (Status)
		{
			Sleep(200);
			Status = ReadFile(hComm,        // Handle to the Serial port
				message,     // Data to be written to the port
				6,  //No of bytes to read
				&dNoOfBytesWritten, //Bytes read
				NULL);
			ccount = strncmp(message, "_OK_", 4);
			if (ccount == 0)
			{
				printf("Pattern transmitted\n");
				Status = true;
			}
			else
			{
				printf("No acknowledgement from ArinaGate \n");
				Status = false;
			}
		}
		else
		{
			printf(">>> Failed transmission to ArinaGate <<< \n");
		}
		packet_counter++;
	}
	

	CloseHandle(hComm);//Closing the Serial Port

	return Status;
}


