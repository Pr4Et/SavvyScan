/*
**************************************************************************

spcm_drv_def.h                                 (c) Spectrum GmbH , 04/2006

**************************************************************************

Contains a definetion of all driver function for external linking
  
Feel free to use this source for own projects and modify it in any kind

**************************************************************************
*/

// driver handling
typedef drv_handle (_stdcall SPCM_HOPEN)              (char* szDeviceName);
typedef void       (_stdcall SPCM_VCLOSE)             (drv_handle hDevice);

// software register handling
typedef uint32     (_stdcall SPCM_DWSETPARAM_I32)     (drv_handle hDevice, int32 lRegister, int32 lValue);           
typedef uint32     (_stdcall SPCM_DWSETPARAM_I64M)    (drv_handle hDevice, int32 lRegister, int32 lValueHigh, uint32 dwValueLow);
typedef uint32     (_stdcall SPCM_DWSETPARAM_I64)     (drv_handle hDevice, int32 lRegister, int64 llValue);

typedef uint32     (_stdcall SPCM_DWGETPARAM_I32)     (drv_handle hDevice, int32 lRegister, int32* plValue);
typedef uint32     (_stdcall SPCM_DWGETPARAM_I64M)    (drv_handle hDevice, int32 lRegister, int32* plValueHigh, uint32* pdwValueLow);
typedef uint32     (_stdcall SPCM_DWGETPARAM_I64)     (drv_handle hDevice, int32 lRegister, int64* pllValue);
  
// data transfer functions
typedef uint32     (_stdcall SPCM_DWDEFTRANSFER_I64M) (drv_handle hDevice, uint32 dwBufType, uint32 dwDirection, uint32 dwNotifySize, void* pvDataBuffer, uint32 dwBrdOffsH, uint32 dwBrdOffsL, uint32 dwTransferLenH, uint32 dwTransferLenL);
typedef uint32     (_stdcall SPCM_DWDEFTRANSFER_I64)  (drv_handle hDevice, uint32 dwBufType, uint32 dwDirection, uint32 dwNotifySize, void* pvDataBuffer, uint64 qwBrdOffs,  uint64 qwTransferLen);
typedef uint32     (_stdcall SPCM_DWINVALIDATEBUF)    (drv_handle hDevice, uint32 dwBufType);

// error handling
typedef uint32     (_stdcall SPCM_DWGETERRORINFO_I32) (drv_handle hDevice, uint32* pdwErrorReg, int32* plErrorValue, char pszErrorTextBuffer[ERRORTEXTLEN]);

// continuous memory
typedef uint32     (_stdcall SPCM_DWGETCONTBUF_I64)   (drv_handle hDevice, uint32 dwBufType, void** ppvDataBuffer, uint64* pqwContBufLen);
typedef uint32     (_stdcall SPCM_DWGETCONTBUF_I64M)  (drv_handle hDevice, uint32 dwBufType, void** ppvDataBuffer, uint32* pdwContBufLenH, uint32* pdwContBufLenL);

