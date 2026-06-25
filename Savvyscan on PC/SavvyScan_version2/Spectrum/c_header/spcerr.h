
// ***********************************************************************
//
// SpcErr.h                                        (c) Spectrum GmbH, 2006
//
// ***********************************************************************
//
// error codes of the Spectrum drivers. Until may 2004 this file was 
// errors.h. Name has been changed because errors.h has been already in 
// use by windows.
//
// ***********************************************************************

#define SPCM_ERROR_ORIGIN_MASK       0x80000000 // this bit marks the origin of the error
#define SPCM_ERROR_ORIGIN_LOCAL      0x00000000 // error occured on local system
#define SPCM_ERROR_ORIGIN_REMOTE     0x80000000 // error occured on remote system (netbox)

#define ERR_OK               0x0000      //   0 No Error
#define ERR_INIT             0x0001      //   1 Initialisation error
#define ERR_NR               0x0002      //   2 Board number out of range
#define ERR_TYP              0x0003      //   3 Unknown board Typ
#define ERR_FNCNOTSUPPORTED  0x0004      //   4 This function is not supported by the hardware
#define ERR_BRDREMAP         0x0005      //   5 The Board Index Remap table is wrong
#define ERR_KERNELVERSION    0x0006      //   6 The kernel version and the dll version are mismatching
#define ERR_HWDRVVERSION     0x0007      //   7 The driver version doesn't match the minimum requirements of the board
#define ERR_ADRRANGE         0x0008      //   8 The address range is disabled (fatal error)
#define ERR_INVALIDHANDLE    0x0009      //   9 Handle not valid
#define ERR_BOARDNOTFOUND    0x000A      //  10 Card with given name hasn't been found
#define ERR_BOARDINUSE       0x000B      //  11 Card with given name is already in use by another application
#define ERR_EXPHW64BITADR    0x000C      //  12 Express hardware version not able to handle 64 bit addressing -> update needed
#define ERR_FWVERSION        0x000D      //  13 Firmware versions of synchronized cards or for this driver do not match -> update needed
#define ERR_SYNCPROTOCOL     0x000E      //  14 Synchronization protocol of synchronized cards does not match -> update needed
#define ERR_KERNEL           0x000F      //  15 Some error occurred in the kernel driver
#define ERR_LASTERR          0x0010      //  16 Old Error waiting to be read
#define ERR_ABORT            0x0020      //  32 Abort of wait function
#define ERR_BOARDLOCKED      0x0030      //  48 Board acess already locked by another process. it's not possible to acess one board through multiple processes
#define ERR_DEVICE_MAPPING   0x0032      //  50 Device is mapped to an invalid device
#define ERR_NETWORKSETUP     0x0040      //  64 Network setup failed
#define ERR_NETWORKTRANSFER  0x0041      //  65 Network data transfer failed
#define ERR_FWPOWERCYCLE     0x0042      //  66 Power cycle needed to update card's firmware (simple PC reboot not sufficient !)
#define ERR_NETWORKTIMEOUT   0x0043      //  67 Network timeout
#define ERR_BUFFERSIZE       0x0044      //  68 Buffer too small
#define ERR_RESTRICTEDACCESS 0x0045      //  69 access to card has been restricted
#define ERR_INVALIDPARAM     0x0046      //  70 invalid parameter for function
#define ERR_TEMPERATURE      0x0047      //  71 card temperature too high

#define ERR_REG              0x0100      // 256 unknown Register for this Board
#define ERR_VALUE            0x0101      // 257 Not a possible value in this state
#define ERR_FEATURE          0x0102      // 258 Feature of the board not installed
#define ERR_SEQUENCE         0x0103      // 259 Channel sequence not allowed
#define ERR_READABORT        0x0104      // 260 Read not allowed after abort
#define ERR_NOACCESS         0x0105      // 261 Access to this register denied
#define ERR_POWERDOWN        0x0106      // 262 not allowed in Powerdown mode
#define ERR_TIMEOUT          0x0107      // 263 timeout occured while waiting for interrupt
#define ERR_CALLTYPE         0x0108      // 264 call type (int32 mux) is not allowed for this register
#define ERR_EXCEEDSINT32     0x0109      // 265 return value is int32 but software register exceeds the 32 bit integer range -> use 2x32 or 64
#define ERR_NOWRITEALLOWED   0x010A      // 266 register cannot be written, read only
#define ERR_SETUP            0x010B      // 267 the setup isn't valid
#define ERR_CLOCKNOTLOCKED   0x010C      // 268 clock section not locked: perhaps no external clock signal connected or not stable
#define ERR_MEMINIT          0x010D      // 269 on-board memory initialization error
#define ERR_POWERSUPPLY      0x010E      // 270 on-board power supply error
#define ERR_ADCCOMMUNICATION 0x010F      // 271 communication with ADC failed
#define ERR_CHANNEL          0x0110      // 272 Wrong number of Channel to be read out
#define ERR_NOTIFYSIZE       0x0111      // 273 Notify block size isn't valid
#define ERR_RUNNING          0x0120      // 288 Board is running, changes not allowed
#define ERR_ADJUST           0x0130      // 304 Auto Adjust has an error
#define ERR_PRETRIGGERLEN    0x0140      // 320 pretrigger length exceeds allowed values
#define ERR_DIRMISMATCH      0x0141      // 321 direction of card and memory transfer mismatch
#define ERR_POSTEXCDSEGMENT  0x0142      // 322 posttrigger exceeds segment size in multiple recording mode
#define ERR_SEGMENTINMEM     0x0143      // 323 memsize is not a multiple of segmentsize, last segment hasn't full length
#define ERR_MULTIPLEPW       0x0144      // 324 multiple pulsewidth counters used but card only supports one at the time
#define ERR_NOCHANNELPWOR    0x0145      // 325 channel pulsewidth can't be OR'd
#define ERR_ANDORMASKOVRLAP  0x0146      // 326 AND mask and OR mask overlap in at least one channel -> not possible
#define ERR_ANDMASKEDGE      0x0147      // 327 AND mask together with edge trigger mode is not allowed
#define ERR_ORMASKLEVEL      0x0148      // 328 OR mask together with level trigger mode is not allowed
#define ERR_EDGEPERMOD       0x0149      // 329 All trigger edges must be simular on one module
#define ERR_DOLEVELMINDIFF   0x014A      // 330 minimum difference between low output level and high output level not reached
#define ERR_STARHUBENABLE    0x014B      // 331 card holding the star-hub must be active for sync
#define ERR_PATPWSMALLEDGE   0x014C      // 332 Combination of pattern with pulsewidht smaller and edge is not allowed
#define ERR_XMODESETUP       0x014D      // 333 The chosen setup for (SPCM_X0_MODE .. SPCM_X19_MODE) is not valid. See hardware manual for details.

#define ERR_NOPCI            0x0200      // 512 No PCI bus found
#define ERR_PCIVERSION       0x0201      // 513 Wrong PCI bus version
#define ERR_PCINOBOARDS      0x0202      // 514 No Spectrum PCI boards found
#define ERR_PCICHECKSUM      0x0203      // 515 Checksum error on PCI board
#define ERR_DMALOCKED        0x0204      // 516 DMA buffer in use, try later
#define ERR_MEMALLOC         0x0205      // 517 Memory Allocation error
#define ERR_EEPROMLOAD       0x0206      // 518 EEProm load error
#define ERR_CARDNOSUPPORT    0x0207      // 519 no support for that card in the library
#define ERR_CONFIGACCESS     0x0208      // 520 error occured during config write or read

#define ERR_FIFOBUFOVERRUN   0x0300      // 768 Buffer overrun in FIFO mode
#define ERR_FIFOHWOVERRUN    0x0301      // 769 Hardware buffer overrun in FIFO mode
#define ERR_FIFOFINISHED     0x0302      // 770 FIFO transfer hs been finished. Number of buffers has been transferred
#define ERR_FIFOSETUP        0x0309      // 777 FIFO setup not possible, transfer rate to high (max 250 MB/s)

#define ERR_TIMESTAMP_SYNC   0x0310      // 784 Synchronisation to ref clock failed
#define ERR_STARHUB          0x0320      // 800 Autorouting of Starhub failed

#define ERR_INTERNAL_ERROR   0xFFFF      // 65535 Internal hardware error detected, please check for update


