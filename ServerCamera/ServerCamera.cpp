// Based on JeolCamera.cpp from another project- class shared between plugin and server

#ifdef CAM_SERVER
#include "stdafxServer.h"
#else
#include "stdafxPlugin.h"
#endif

#include "ServerCamera.h"
#include "BaseServer.h"
#include "../Spectrum/shahar/RecRepMulti.h"
#include <math.h>
#define MSG_STR_SIZE 320
#define IMAGE_MUTEX_WAIT  1000
#define CONTINUOUS_TIMEOUT 2000
#define SETTING_INDEX 9

// Even when running in a DLL, we need our own copy of the module, so here it is for both
CComModule _Module;

static bool sComInitialized = false;
static bool sObjInitialized = false;
static bool sCamObjExists = false;
static bool sScanObjExists = false;
static bool sCamInitialized = false;
static int sScanInitialized = 0;
static char sMessageStr[MSG_STR_SIZE];
static char sErrorStr[MSG_STR_SIZE];
static int sCurrentCam = 0;
static short *sLiveBuffer = NULL;
static int sLiveWidth, sLiveHeight, sArrSize;
static int sDivideBy2;
static int sLiveError = 0;
static int sFrameReady = 0;
static int sWaitingForFrame = 0;
static HANDLE sImageMutexHandle;
static int sDebug = 0;

static BOOL SleepMsg(DWORD dwTime_ms);
static int ProcessFrameImage(_FrameImageCom *frameImage, short *array, int arrSize, 
                             int divideBy2, int &width, int &height);
static  double flyback_us=40.0; //flyback time imposed in ScanMode=1 (shahar)
const double PI = 3.141592653589793;

 //shahar, global vairbales from Shadow GUI to RecRepMulti.cpp and ServerCamera.cpp
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
 extern bool GUI_ch_update_ready;
 extern bool GUI_tomography;
 extern int GUI_tomoIndex;
 extern int GUI_numberOfSlices;
 extern int GUI_chosenCH;
 extern int GUI_cella;

// Convenience functions for error and debug output
void ErrorFmt(char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vsprintf_s(sErrorStr, MSG_STR_SIZE, fmt, args);
  va_end(args);
  ErrorToLog(sErrorStr);
}

void DebugFmt(char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vsprintf_s(sMessageStr, MSG_STR_SIZE, fmt, args);
  va_end(args);
  DebugToLog(sMessageStr);
}

void ReportCOMError(_com_error E, const char *context)
{
  char messBuf[MSG_STR_SIZE];

  // These were %ls, when project was unicode
  if (E.Description().length() > 0) {
    sprintf_s(messBuf, MSG_STR_SIZE, "COM ERROR in call to %s: %s\n", 
      context, (char *)E.Description());
  } else if (E.ErrorMessage() != NULL) {
    sprintf_s(messBuf, MSG_STR_SIZE, "COM ERROR in call to %s: %s\n", 
      context, (char *)E.ErrorMessage());
  } else {
    sprintf_s(messBuf, MSG_STR_SIZE, 
      "COM ERROR in call to %s: No description\n", context);
  }
  EitherToLog("", messBuf, true);
}


CJediCamera::CJediCamera(void)
: m_pCameraComInterfaceClass( NULL)
, m_pCameraSettingCom( NULL)
, m_pFrameImageCom( NULL)
, m_pEventReceiver( NULL)
, mNumCameras(0)
, mUserSettingIndex(-1)
{
  sImageMutexHandle = CreateMutex(0, 0, 0);
}

CJediCamera::~CJediCamera(void)
{
}

int CJediCamera::InitializeObjects()
{ 
  HRESULT hr;
  int retVal = 0, retVal2 = 0;
  if (!sComInitialized)
    _Module.Init(NULL, GetModuleHandle(NULL));
  sComInitialized = true;
  if (sObjInitialized)
    return 0;


  if (retVal && retVal2)
    UninitializeObjects();
  else {
    mNumCameras = 2;
    sObjInitialized = true;
    sCamObjExists = !retVal;
    sScanObjExists = !retVal2;
  }
  return retVal + retVal2;
}

int CJediCamera::UninitializeObjects()
{
  sObjInitialized = false;
  if (sComInitialized)
    _Module.Term();
  sComInitialized = false;
  return 0;
}

/////////////////////////////////////
//   THE Function Calls

int CJediCamera::GetNumberOfCameras()
{
  return mNumCameras;
}

int CJediCamera::InitializeCamera(int camNum)
{
  int error = 0;
  long ind, uBound;
  int settingIndex = 9;
  HRESULT hr;
  VARTYPE pvt;
  SAFEARRAY *pSa;
  //FRAMESIZE size; //shahar moved to 
  long scanSpeeds[] = {0, 1, 2, 3, 4, 4, 4};
  long speedLevels[] = {8, 8, 8, 8, 8, 9, 10};
  if (!sObjInitialized || camNum < 0 || camNum > MAX_SCAN_CHANNELS)
    return 1;
  if (!camNum) {

    // Index 0 is CCD camera
    if (sCamInitialized)
      return 0;
    if (!sCamObjExists)
      return 1;
	//index 1 is out STEM camera
    sScanInitialized = camNum;
  }
  else
	sScanInitialized=Initialize_cards ();//1 is successful
  return 0;
}
  
int CJediCamera::GetSTEMProperties(int useSizeX, int useSizeY, double *minPixel, 
                                   double *maxPixel, double *pixelInc, 
                                   double *rotationInc, double *ddum, int *maxIntegration,
                                   int *idum)
{
  if (!sScanInitialized)
    return 1;
  if (useSizeX > 0 && useSizeY > 0) {
    mScanSizeX = useSizeX;
    mScanSizeY = useSizeY;
  }
  *minPixel = mMinPixelTime; 
  *maxPixel = mMaxPixelTime; 
  *pixelInc = mPixelTimeInc;
  *rotationInc = mRotationStep;
  *ddum = 0.;
  *maxIntegration = mMaxFrameIntegration;
  *idum = 0;
  return 0;
}

int CJediCamera::UninitializeCameras(void)
{
  if (sCamInitialized) {
    sCamInitialized = false;
    m_pCameraComInterfaceClass->raw_TemUnInitialize(); 
  }
  if (sScanInitialized) {
      Close_cards();
	  sScanInitialized = 0;
    //m_pScanningComInterface[0]->raw_UnInitialize();//removed shahar
  }
  return 0;
}

int CJediCamera::SelectCamera(int camNum)
{
  if (camNum / 2) 
    return 1;
  sCurrentCam = camNum;
  return 0;
}

int CJediCamera::SetDebugMode(int value)
{
  sDebug = value;
  return 0;
}


int CJediCamera::AcquireCCDImage(short int *array, long *arrSize, long *sizeX, 
                                 long *sizeY, int top, int left, int binning, 
                                 int processing, double exposure, int flags)
{
    return 1;
}

// Acquire STEM image given binned coordinates and size, and binning
//SavvyScan support: Shahar Seifer
int CJediCamera::AcquireSTEMImage(short **arrays, long *arrSize, long *sizeX, long *sizeY,
                                  int top, int left, int binning, double exposure,
                                  double rotation, int integration, int numChan, 
                                  int *channels, int flags, int *numAcquired)
{
  //exposure is the pixel time in us seen in the SerialEM GUI, binning and sizeX/Y are exactly like in the SerialEM GUI
  DWORD startTime = GetTickCount();
  HRESULT hr;
  int error = 0, width, height, scTop, scLeft, scXsize, scYsize, ind, wait;
  int camLeft, camTop, camXsize, camYsize;
  int fsize=0; //full sampling point size initially without oversampling, then overwriten with net size
  int ScanMode=GUI_scanmode; //shahar
  sDivideBy2 = flags & 1;
  DWORD selectSleep = 3000;
  bool doContinuous = (flags & 2) != 0;
  bool setScanSetArea;
  bool running, anyChanged = false;
  int imageMode = 0;
  int settingIndex = -1;
  if (!integration)
    integration = 1;
  if (!sScanInitialized) {
    ErrorToLog("You need to restart SerialEM to reinitialize the SavvyScan access");
    return 1;
  }

  if (!doContinuous && sLiveBuffer)
    StopContinuous();
  if (doContinuous && sLiveBuffer) {
    *numAcquired = 1;
    return ReturnLiveFrame(arrays[0], false, sizeX, sizeY);
  }
  if (flags & 8)
    settingIndex = flags >> 24;

   width = *sizeX;    //net horizontal(x) size
  height = *sizeY;  //net vertical (y) size
  scTop = top;
  scLeft = left;
  if (ScanMode==1)
  {
	  scXsize = (int)((*sizeX)*(1+0.414));//was 1+0.414 before sep17 12:00
	  scYsize = *sizeY;
 	  fsize =  height*scXsize+height*(int)(flyback_us/exposure); //must be in coordination with prepare_AWG(..)

  }
  else if (ScanMode==2)
  {
	  scXsize = (int)(width*sqrt(2.0)); //bruto size
	  scYsize = (int)(height*sqrt(2.0));
	  //simulate spiral with oversampling=1 to calculate fsize
	  int rnbs=(int)(scXsize/2);//number of radius points
   	  int phnmbs;//number of phase points
	  int counter=0;
	  for (int scanrvar=rnbs;scanrvar>=0; scanrvar--) //radius of circle, enclosing a square widthXwidth
	  {
			if (scanrvar!=0)
				phnmbs=(int)(2*PI*scanrvar); //so the arc steps are 1pixel long
			else
				phnmbs=1;
			counter=counter+phnmbs;
	  }
	  fsize=counter;
  }
  else if (ScanMode==3)
  {
	  //scXsize = (int)(3*width); //overall size in the case of 2/3 AspectRatio: width + two circles of same size on the left and right
	  scXsize = (int)(2*width); //overall size in normal case: width + two half circles of same diameter on the left and right
	  scYsize = (int)height;
	  //simulate with oversampling=1 to calculate fsize
		//int rnbs=(int)(width);//number of circling times in the case of 2/3 AspectRatio
		int rnbs=(int)(width/2);//number of circling times, so the right hand of the circle interlaces with left hand, jumping 2 pixels per cycle, drifting width=1 diameters
		int scanr=(int)(height/2); //radius of the circle
		int phnmbs=((int)(2.0*PI*scanr)); ////number of phase points in one round: so the arc steps are 1pixel, oversampling here is 1
		int counter=phnmbs*rnbs;
		fsize=counter;
  }
  setScanSetArea =true;


  if (*sizeX < mMinScanSizeX || *sizeY < mMinScanSizeY) {
    ErrorFmt("Scan image size of %d x %d is too small, must be at least %d x %d",
      scXsize, scYsize, mMinScanSizeX, mMinScanSizeY);
    return 1;
  }
  if (top < 0 || left < 0 || top + *sizeY > mScanSizeY / binning || 
    left + *sizeX > mScanSizeX / binning) {
      ErrorFmt("The binned scan coordinates (X %d to %d, Y %d to %d) are outside the "
        "range for the device (%d x %d)", left, left + *sizeX, top, top + *sizeY,
        mScanSizeX / binning, mScanSizeY / binning);
      return 1;
  }

  try {

    if (anyChanged){
      Sleep(selectSleep);

      DebugFmt("Time for setting camera = %.0f  at %d", TickInterval(startTime), 
        (GetTickCount() % 3600000));
      startTime = GetTickCount();
    }
  
    // Get the frames
    *numAcquired = 0;
    bool simulate=false;
	if (!simulate) //Here we activate STEM by ADC-DAC cards
	{		
		int oversampling,pass,fail;
		pass=set_cards_samplerate_andmore(exposure, fsize, oversampling); //calculate samplerate,oversampling, datasize and set cards accordingly, start timestamp DMA
		if (oversampling<=0 || !pass) {
			RestoreUserSettings();
			return false; //error
		}
		prepare_AWG(binning, scXsize, scYsize, GUI_aspectratio , rotation, oversampling, exposure, ScanMode, flyback_us, width, height); //prepare scan data, mode 1: like Gatan, flyback=40us
		fail=Acquire_scan(); //copy DAC buffer to card, run the cards, copy ADC data to buffer
		if (fail) {
			RestoreUserSettings();
			return error;
		}
		//Rearrnage data to different channels (all channels already acquired in parallel 
		fsize =  (width * height); //overwrite previous value
		//Note that in SERVER-SEMCamserver.cpp it is aleady defined: sChanArrays[ind] = new short[sArrSize], and its pointer is passed to arrays pointer here.
		retrieve_images(width,height,arrays[0]); //prepare images, save all in mrc files, pass one image to arrays[0] pointer of image
		(*numAcquired)+=1;
		error=0;
	
	} else
	{		
		//simulation
		for (ind = 0; ind < numChan; ind++) {
			//shahar: acquire images to arrays[0][], of length width X height, short (2 bytes) per pixel 
		    char filename[]="D:\\Savvyscan-source\\test.bmp";
			FILE* f = fopen(filename, "rb");
			unsigned char info[1024+54];
			// read the 54-byte header of BMP file
			fread(info, sizeof(unsigned char), 1024+54, f); 
			// extract image height and width from header
			width = *(int*)&info[18]; 
			height = *(int*)&info[22]; 
			// allocate 1 bytes per pixel (grayscale image)
			fsize =  (width * height); 
			unsigned char* data = new unsigned char[fsize];
			short *pdata=new short[fsize]; 
			// read the rest of the data at once
			fread(data, sizeof(unsigned char), fsize, f); 
			fclose(f);
			//copy pointers for each pixel since we trandform from uchar to short with different sizes in memory
			for (int i=0; i<width*height; i++)
			{
				pdata[i]=(short)data[i];
			}
			//memcpy(arrays[ind],pdata, width * height*sizeof(short)); //copy physical bytes
			arrays[ind]=pdata; //copy pointer to array of pointers
			(*numAcquired)+=1;
			error=0;
		}
	}
    //should be: *numAcquired = numChan;
    RestoreUserSettings();
  }
  catch (_com_error E) {
    ReportCOMError(E, "Error setting or accessing a variable in STEM interface");
    RestoreUserSettings();
    return 1;
  }
  *sizeX = width;
  *sizeY = height;
  *arrSize=(long)width*height;
  return 0;
}

void CJediCamera::RestoreUserSettings()
{
}

int CJediCamera::ReturnLiveFrame(short *array, bool firstTime, long *sizeX, long *sizeY)
{
  int error = 0;
  double startTime = GetTickCount();
  int gotmut, start = GetTickCount() % 3600000;
  sWaitingForFrame = 1;
  for (; ;) {
    if (sFrameReady || (!firstTime && TickInterval(startTime) > CONTINUOUS_TIMEOUT)) {

      WaitForSingleObject(sImageMutexHandle, IMAGE_MUTEX_WAIT);
      gotmut = GetTickCount() % 3600000;
      if (sLiveError) {
        StopContinuous();
        error = sLiveError;
      }
      if (!error && !sLiveBuffer)
        error = 1;
      if (!error) {
        memcpy(array, sLiveBuffer, sLiveWidth * sLiveHeight * sizeof(short));
        *sizeX = sLiveWidth;
        *sizeY = sLiveHeight;
      }
      sFrameReady = 0;
      sWaitingForFrame = 0;
      ReleaseMutex(sImageMutexHandle);
      return error;
    }
    SleepMsg(10);
  }
}


int CJediCamera::StopContinuous()
{
  int retVal = 0;
  DebugFmt("Called LiveStopped at %d", GetTickCount() % 3600000);
  delete [] sLiveBuffer;
  sLiveBuffer = NULL;
  if (sCurrentCam > 0)
    RestoreUserSettings();
  return retVal;
}


double TickInterval(double start)
{
  double interval = GetTickCount() - start;
  if (interval < 0)
    interval += 4294967296.;
  return interval;
}

/*
 * sleeps for the given amount of time while pumping messages
 * returns TRUE if successful, FALSE otherwise
 */
static BOOL SleepMsg(DWORD dwTime_ms)
{
  DWORD dwStart = GetTickCount();
  DWORD dwElapsed;
  while ((dwElapsed = GetTickCount() - dwStart) < dwTime_ms) {
    DWORD dwStatus = MsgWaitForMultipleObjects(0, NULL, FALSE,
                                               dwTime_ms - dwElapsed, QS_ALLINPUT);
    
    if (dwStatus == WAIT_OBJECT_0) {
      MSG msg;
      while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
        if (msg.message == WM_QUIT) {
          PostQuitMessage((int)msg.wParam);
          return FALSE; // abandoned due to WM_QUIT
        }

        TranslateMessage(&msg);
        DispatchMessage(&msg);
      }
    }
  }

  return TRUE;
}
