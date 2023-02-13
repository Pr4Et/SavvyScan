#ifndef MD5_H
#define MD5_H
#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"


typedef struct 
    {
    uint32 dwState0;                                   
    uint32 dwState1;                                 
    uint32 dwState2;                                  
    uint32 dwState3;                                  
    } MD5_STATES;

typedef struct 
    {
    uint32 dwCount0;        
    uint32 dwCount1;
    } MD5_COUNTS;

typedef struct 
    {
    uint8 pcBuffer[64];                        
    } MD5_BUF;



void vMD5_Calculation (uint8 pbyCheckSum [16], MD5_STATES* pstStates, MD5_COUNTS* pstCounts, MD5_BUF* pstBuf);
void vStartCalc (MD5_STATES* pstStates, uint8 pcBuf[64]);

#endif
