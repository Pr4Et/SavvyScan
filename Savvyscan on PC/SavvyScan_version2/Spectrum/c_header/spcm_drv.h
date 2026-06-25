/*
**************************************************************************

spcm_interface.h                               (c) Spectrum GmbH , 08/2005

**************************************************************************

Interface of the spectrum driver for all M2I cards. 

**************************************************************************
*/



/*
**************************************************************************

open/close: 
opens and closes the specified board. The returned handle is used by all
function calls. The device name is operating system specific. Under
linux it's normally "/dev/spcm0" for card#0 and under windows it's "spc0"
In a synchronized system the starhub has it's own handle, received under
the device name "starhub"

***************************************************************************
*/

// tries to open the device and returns handle or error code
SPCM_IMPORT drv_handle _stdcall spcm_hOpen (
    const char* szDeviceName);          // name of the device to be opened

//*************************************************************************

// closes the device
SPCM_IMPORT void _stdcall spcm_vClose (             
    drv_handle  hDevice);               // handle to an already opened device





/*
**************************************************************************

SetParam and GetParam: 
handles the register based access to the driver. Each functionality is
programmed by accessing one of the software registers of the driver

Functions are declared as i32 for single 32 bit integer access, i64 for 
single 64 bit integer access or 64m for multiplexed 64 bit integer access 
consisting one 32 bit integer high-part and one 32 bit unsigned integer 
low part. 

Most software registers are only handled by the i32 or i64 function and 
are then not allowed to be accessed by i64m functions. Some registers 
can be more than 32 bit wide. They can be accessed by any of the 
functions. If accessed by the i32 function the value is limited to the 
32 bit signed integer area.

***************************************************************************
*/

// Sets a software register using 1 x 32 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwSetParam_i32 (   
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be modified
    int32       lValue);                // the value to be set

//*************************************************************************

// Sets a software register using 2 x 32 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwSetParam_i64m (  
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be modified
    int32       lValueHigh,             // upper 32 bit of the value. Containing the sign bit !
    uint32      dwValueLow);            // lower 32 bit of the value.

//*************************************************************************

// Sets a software register using 1 x 64 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwSetParam_i64 (   
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be modified
    int64       llValue);               // the value to be set



//*************************************************************************
//*************************************************************************

// Reads out a software register using 1 x 32 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwGetParam_i32 (   
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be read out
    int32*      plValue);               // pointer for the return value

//*************************************************************************

// Reads out a software register using 2 x 32 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwGetParam_i64m (  
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be read out
    int32*      plValueHigh,            // pointer for the upper part of the return value
    uint32*     pdwValueLow);           // pointer for the lower part of the return value

//*************************************************************************

// Reads out a software register using 1 x 64 bit integer value. Return value is an error code
SPCM_IMPORT uint32 _stdcall spcm_dwGetParam_i64 (   
    drv_handle  hDevice,                // handle to an already opened device
    int32       lRegister,              // software register to be read out
    int64*      pllValue);              // pointer for the return value




/*
**************************************************************************

DefTransfer:
sets up all needed information for the next data transfer. Data transfer
itself is started by an extra register command.

The function needs 64 bit unsigned integer values. Therefore it is 
available as an i64m type, consisting of one upper 32 bit uint and one
lower 32 bit uint value. And it is availabl as a true 64 bit version.

Offset and length are both given in samples. As data is multiplexed the
transfer buffer in PC memory must be large enough to handle 
[length x channels] entries

***************************************************************************
*/

// defintions of the transfer direction
#define SPCM_DIR_PCTOCARD   0           // transfer from PC memory to card memory
#define SPCM_DIR_CARDTOPC   1           // transfer from card memory to PC memory
#define SPCM_DIR_CARDTOGPU  2           // RDMA transfer from card memory to GPU memory
#define SPCM_DIR_GPUTOCARD  3           // RDMA transfer from GPU memory to card memory

// defintions of the different data buffers
#define SPCM_BUF_DATA       1000        // main data buffer for acquired or generated samples
#define SPCM_BUF_ABA        2000        // buffer for ABA data, holds the A-DATA (slow samples)
#define SPCM_BUF_TIMESTAMP  3000        // buffer for timestamps
#define SPCM_BUF_LOG        4000        // write content of buffer to log file

//*************************************************************************

// Defines the transer buffer by using 2 x 32 bit unsigned integer values for each 64 bit value
SPCM_IMPORT uint32 _stdcall spcm_dwDefTransfer_i64m(
    drv_handle  hDevice,                // handle to an already opened device
    uint32      dwBufType,              // type of the buffer to define as listed above under SPCM_BUF_XXXX
    uint32      dwDirection,            // the transfer direction as defined above
    uint32      dwNotifySize,           // amount of bytes after which i want do receive an event (0=end of transfer) 
    void*       pvDataBuffer,           // pointer to the data buffer
    uint32      dwBrdOffsH,             // high part of offset in board memory
    uint32      dwBrdOffsL,             // low part of offset in board memory
    uint32      dwTransferLenH,         // high part of transfer buffer length
    uint32      dwTransferLenL);        // low part of transfer buffer length

//*************************************************************************

// Defines the transer buffer by using 64 bit unsigned integer values
SPCM_IMPORT uint32 _stdcall spcm_dwDefTransfer_i64 (
    drv_handle  hDevice,                // handle to an already opened device
    uint32      dwBufType,              // type of the buffer to define as listed above under SPCM_BUF_XXXX
    uint32      dwDirection,            // the transfer direction as defined above
    uint32      dwNotifySize,           // amount of bytes after which i want do receive an event (0=end of transfer) 
    void*       pvDataBuffer,           // pointer to the data buffer
    uint64      qwBrdOffs,              // offset for transfer in board memory
    uint64      qwTransferLen);         // buffer length

//*************************************************************************

// invalidate the transfer buffer (is automatically performed if new transfer buffer is defined with DefTransfer)
SPCM_IMPORT uint32 _stdcall spcm_dwInvalidateBuf (  
    drv_handle  hDevice,                // handle to an already opened device
    uint32      dwBufType);             // type of the buffer to invalidate as listed above under SPCM_BUF_XXXX



/*
**************************************************************************

GetContBuf
reads out the internal continuous memory buffer if one has been allocated
this continuous buffer allows faster data transfer especially on Express cards

***************************************************************************
*/

SPCM_IMPORT uint32 _stdcall spcm_dwGetContBuf_i64 (
    drv_handle  hDevice,                // handle to an already opened device
    uint32      dwBufType,              // type of the buffer to read as listed above under SPCM_BUF_XXXX
    void**      ppvDataBuffer,          // address of available data buffer
    uint64*     pqwContBufLen);         // length of available continuous buffer

//*************************************************************************

SPCM_IMPORT uint32 _stdcall spcm_dwGetContBuf_i64m (
    drv_handle  hDevice,                // handle to an already opened device
    uint32      dwBufType,              // type of the buffer to read as listed above under SPCM_BUF_XXXX
    void**      ppvDataBuffer,          // address of available data buffer
    uint32*     pdwContBufLenH,         // high part of length of available continuous buffer
    uint32*     pdwContBufLenL);        // low part of length of available continuous buffer



/*
**************************************************************************

GetErrorInfo:
reads out the complete error information that is stored in the driver. 
internal error locking is afterwards reset.
If hDevice is zero the last open error is returned.

***************************************************************************
*/

SPCM_IMPORT uint32 _stdcall spcm_dwGetErrorInfo_i32 (
    drv_handle  hDevice,                // handle to an already opened device
    uint32*     pdwErrorReg,            // adress of the error register (can zero if not of interest)
    int32*      plErrorValue,           // adress of the error value (can zero if not of interest)
    char        pszErrorTextBuffer[ERRORTEXTLEN]); // text buffer for text error 



/*
**************************************************************************

StartEBox:
starts the ethernet box by instanciate the kernelhandle_lan_ebox.
The ethernet box is the server with a card with is waiting for 
request by lan-client (host pc).

***************************************************************************
*/

SPCM_IMPORT uint32 _stdcall spcm_dwStartEBox ();



/*
**************************************************************************

dwDiscovery:
the lan-client (host pc) starts a broadcast request and wait for a answer 
with the VISA string. 
A VISA string contains the ip address of the ethernet box.

***************************************************************************
*/

SPCM_IMPORT uint32 _stdcall spcm_dwDiscovery (
    char** pszVisaString,           // user-allocated array of C-strings to return the VISA strings
    uint32 dwMaxNoOfDevices,        // the maximum number of devices to be returned
    uint32 dwMaxVisaLen,            // maximum length of one entry in pszVisaString
    uint32 dwTimeout);              // time in milli seconds that the function will wait until each device has answered.



/*
**************************************************************************

dwWriteIDNRequest:
the lan-client (host pc) sends a "IDN" request to the ethernet box.
So the ethernet box will send a message which contains the manufacturer, 
the model, the serial number and the firmware version

***************************************************************************
*/

SPCM_IMPORT uint32 _stdcall spcm_dwSendIDNRequest (
    char** szIdnString,                 //the IDN string looks like that: <manufacturer>,<model>,<serial number>,<firmware version>
    uint32 dwNoOfDevices, 
    uint32 dwIdnStringLen);            //string which contains manufacturer, the model, the serial number and the firmware version
