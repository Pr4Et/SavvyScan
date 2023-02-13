#include "spcm_md5.h"

#include <math.h>


#define PUINT8_TO_UINT32(dwOut, pcIn) (dwOut = (pcIn[3] << 24) | (pcIn[2] << 16) | (pcIn[1] << 8) | pcIn[0])

#define UINT32_TO_PUINT8(pcOut, dwIn) { \
    pcOut[0] = (uint8) ((dwIn >>  0) & 0x000000ff); \
    pcOut[1] = (uint8) ((dwIn >>  8) & 0x000000ff); \
    pcOut[2] = (uint8) ((dwIn >> 16) & 0x000000ff); \
    pcOut[3] = (uint8) ((dwIn >> 24) & 0x000000ff); \
}

#define SIZEOF_UINT32               sizeof (uint32)
#define SIZEOF_ORIGIN_MSG_LEN       8                   //always 64 Bit->8Bytes
#define SHIFT3                      3
#define MASK_63BIT                  0x3F
#define SMALLEST_RESULT_MSG         64                  //after padding with the message len, the result ist a multiple of SMALLEST_RESULT_MSG (in Bytes)
#define MSG_WITHOUT_MSGLEN          SMALLEST_RESULT_MSG - SIZEOF_ORIGIN_MSG_LEN

//max 64 Bytes = 512 Bits, in any case there will be pad the 1
//the message is padding while the message len % 512 (bits) is 64 (bits). 
//The 64 bit will be used for the len of the origin message
static uint8 PADDING_BYTES_TO_MESSAGE [64] = 
    {
    0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

//this values are given by the algorithm
static uint32 r [64] = 
    {
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,         //round 1
    5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,             //round 2
    4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,         //round 3
    6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21          //round 4
    };

//hard coded sin values - this is for timing
static uint32 k [64] = 
    {
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501, //sin table entries for round 1
    0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
    0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x2441453, 0xd8a1e681, 0xe7d3fbc8,  //sin table entries for round 2
    0x21e1cde6, 0xc33707d6, 0xf4d50d87,    0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
    0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, //sin table entries for round 3
    0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x4881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1, //sin table entries for round 4
    0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391
    }; 



//the four functions for the MD5 algorithm
#define CALC_F1(dwState1, dwState2, dwState3) (((dwState1) & (dwState2)) | ((~dwState1) & (dwState3)))
#define CALC_F2(dwState1, dwState2, dwState3) (((dwState1) & (dwState3)) | ((dwState2) & (~dwState3)))
#define CALC_F3(dwState1, dwState2, dwState3) ((dwState1) ^ (dwState2) ^ (dwState3))
#define CALC_F4(dwState1, dwState2, dwState3) ((dwState2) ^ ((dwState1) | (~dwState3)))

//every bit will be shift. if the bit would fall out of the integer it will be set as lsb
#define ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32-(n))))

// functions for round 1-4
#define CALC_STATE_ROUND1(dwState0, dwState1, dwState2, dwState3, x, s, ac) { \
(dwState0) = ROTATE_LEFT (((dwState0) + CALC_F1 ((dwState1), (dwState2), (dwState3)) + (x) + (uint32)(ac)), (s)) + (dwState1); \
}
#define CALC_STATE_ROUND2(dwState0, dwState1, dwState2, dwState3, x, s, ac) { \
(dwState0) = ROTATE_LEFT (((dwState0) + CALC_F2 ((dwState1), (dwState2), (dwState3)) + (x) + (uint32)(ac)), (s)) + (dwState1); \
}
#define CALC_STATE_ROUND3(dwState0, dwState1, dwState2, dwState3, x, s, ac) { \
(dwState0) = ROTATE_LEFT (((dwState0) + CALC_F3 ((dwState1), (dwState2), (dwState3)) + (x) + (uint32)(ac)), (s)) + (dwState1); \
}
#define CALC_STATE_ROUND4(dwState0, dwState1, dwState2, dwState3, x, s, ac) { \
(dwState0) = ROTATE_LEFT (((dwState0) + CALC_F4 ((dwState1), (dwState2), (dwState3)) + (x) + (uint32)(ac)), (s)) + (dwState1); \
}



//makes the padding of the message. On a message will pad a '1' and so much '0' until the length of the message modulo 512 (bits) is 448 (512-64 bits). 
//The last 64 bit are for the origin len of the message.
//Here the calculation of the algorith is done too.
void vMD5_Calculation (uint8 pbyCheckSum [16], MD5_STATES* pstStates, MD5_COUNTS* pstCounts, MD5_BUF* pstBuf)
    {
    uint32 i = 0;
    uint8 pcBits[8];
    uint8 pcOut[4];
    uint32 dwIn = 0;
    uint32 dwIndex = 0;
    uint32 dwPaddingLen = 0;
    uint32 dwPartLen = 0;


    // Save number of bits 
    dwIn = pstCounts->dwCount0;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pcBits, pcOut, SIZEOF_UINT32);
    dwIn = pstCounts->dwCount1;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pcBits+4, pcOut, SIZEOF_UINT32);

  
    // Pad out to 56 mod 64.
    dwIndex = (uint32)((pstCounts->dwCount0 >> SHIFT3) &MASK_63BIT);      //runterrechnung auf bytes
    dwPaddingLen = (dwIndex < MSG_WITHOUT_MSGLEN) ? (MSG_WITHOUT_MSGLEN - dwIndex) : ((MSG_WITHOUT_MSGLEN + SMALLEST_RESULT_MSG) - dwIndex);           //länge die erweitert werden soll


    //Update number of bits
    pstCounts->dwCount0 += ((uint32)dwPaddingLen << SHIFT3);
    dwPartLen = SMALLEST_RESULT_MSG - dwIndex;

    //Transform as many times as possible.
    if (dwPaddingLen >= dwPartLen) 
        {
        memcpy ((uint8*)&pstBuf->pcBuffer[dwIndex], (uint8*)PADDING_BYTES_TO_MESSAGE, dwPartLen);
        vStartCalc (pstStates, pstBuf->pcBuffer);
        for (i = dwPartLen; i + (SMALLEST_RESULT_MSG - 1) < dwPaddingLen; i += SMALLEST_RESULT_MSG)
            vStartCalc (pstStates, &PADDING_BYTES_TO_MESSAGE[i]);
        dwIndex = 0;
        }
    else
        i = 0;

    //pad the bits (1 [0 ... 0])
    memcpy ((uint8*)&pstBuf->pcBuffer[dwIndex], (uint8*)&PADDING_BYTES_TO_MESSAGE[i], dwPaddingLen-i);

    //recalculate the Index and the count values
    dwIndex = (uint8)((pstCounts->dwCount0 >> SHIFT3) & MASK_63BIT);

    if ((pstCounts->dwCount0 += ((uint32)SIZEOF_ORIGIN_MSG_LEN << SHIFT3)) < ((uint32)8 << SHIFT3))
        pstCounts->dwCount1++;
    pstCounts->dwCount1 += ((uint32)SIZEOF_ORIGIN_MSG_LEN >> 29);
    dwPartLen = 64 - dwIndex;

    memcpy ((uint8*)&pstBuf->pcBuffer[dwIndex], (uint8*)pcBits, dwPartLen);
    vStartCalc (pstStates, pstBuf->pcBuffer);
    for (i = dwPartLen; i + (SMALLEST_RESULT_MSG - 1) < SIZEOF_ORIGIN_MSG_LEN; i += SMALLEST_RESULT_MSG)
        vStartCalc (pstStates, &pcBits[i]);
    dwIndex = 0;

    memcpy ((uint8*)&pstBuf->pcBuffer[dwIndex], (uint8*)&pcBits[i], SIZEOF_ORIGIN_MSG_LEN);

    dwIn = pstStates->dwState0;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pbyCheckSum, pcOut, SIZEOF_UINT32);
    dwIn = pstStates->dwState1;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pbyCheckSum+4, pcOut, SIZEOF_UINT32);
    dwIn = pstStates->dwState2;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pbyCheckSum+8, pcOut, SIZEOF_UINT32);
    dwIn = pstStates->dwState3;
    UINT32_TO_PUINT8 (pcOut, dwIn);
    memcpy (pbyCheckSum+12, pcOut, SIZEOF_UINT32);

    //reset structs
    dwPartLen = sizeof (*pstStates);
    memset ((uint8*)pstStates, 0, dwPartLen);
    dwPartLen = sizeof (*pstCounts);
    memset ((uint8*)pstCounts, 0, dwPartLen);
    dwPartLen = sizeof (*pstBuf);
    memset ((uint8*)pstBuf, 0, dwPartLen);
    }

//main calculation of this algorithm
void vStartCalc (MD5_STATES* pstStates, uint8 pcBuf[64])
    {
    uint32 dwState0 = pstStates->dwState0;
    uint32 dwState1 = pstStates->dwState1;
    uint32 dwState2 = pstStates->dwState2;
    uint32 dwState3 = pstStates->dwState3;
    uint32 x[16];
    uint32 i = 0;
    uint8 pcIn[4];
    uint32 dwOut = 0;


    for (i = 0; i < 16; i++)
        {
        memcpy (pcIn, pcBuf + i*SIZEOF_UINT32, SIZEOF_UINT32);
        PUINT8_TO_UINT32 (dwOut, pcIn);
        x[i] = dwOut;
        }

    CALC_STATE_ROUND1 (dwState0, dwState1, dwState2, dwState3, x[ 0], r[ 0], k[ 0]); 
    CALC_STATE_ROUND1 (dwState3, dwState0, dwState1, dwState2, x[ 1], r[ 1], k[ 1]); 
    CALC_STATE_ROUND1 (dwState2, dwState3, dwState0, dwState1, x[ 2], r[ 2], k[ 2]); 
    CALC_STATE_ROUND1 (dwState1, dwState2, dwState3, dwState0, x[ 3], r[ 3], k[ 3]); 
    CALC_STATE_ROUND1 (dwState0, dwState1, dwState2, dwState3, x[ 4], r[ 4], k[ 4]); 
    CALC_STATE_ROUND1 (dwState3, dwState0, dwState1, dwState2, x[ 5], r[ 5], k[ 5]); 
    CALC_STATE_ROUND1 (dwState2, dwState3, dwState0, dwState1, x[ 6], r[ 6], k[ 6]); 
    CALC_STATE_ROUND1 (dwState1, dwState2, dwState3, dwState0, x[ 7], r[ 7], k[ 7]); 
    CALC_STATE_ROUND1 (dwState0, dwState1, dwState2, dwState3, x[ 8], r[ 8], k[ 8]); 
    CALC_STATE_ROUND1 (dwState3, dwState0, dwState1, dwState2, x[ 9], r[ 9], k[ 9]); 
    CALC_STATE_ROUND1 (dwState2, dwState3, dwState0, dwState1, x[10], r[10], k[10]); 
    CALC_STATE_ROUND1 (dwState1, dwState2, dwState3, dwState0, x[11], r[11], k[11]); 
    CALC_STATE_ROUND1 (dwState0, dwState1, dwState2, dwState3, x[12], r[12], k[12]); 
    CALC_STATE_ROUND1 (dwState3, dwState0, dwState1, dwState2, x[13], r[13], k[13]); 
    CALC_STATE_ROUND1 (dwState2, dwState3, dwState0, dwState1, x[14], r[14], k[14]); 
    CALC_STATE_ROUND1 (dwState1, dwState2, dwState3, dwState0, x[15], r[15], k[15]); 
    
    //2nd round
    CALC_STATE_ROUND2 (dwState0, dwState1, dwState2, dwState3, x[ 1], r[16], k[16]); 
    CALC_STATE_ROUND2 (dwState3, dwState0, dwState1, dwState2, x[ 6], r[17], k[17]); 
    CALC_STATE_ROUND2 (dwState2, dwState3, dwState0, dwState1, x[11], r[18], k[18]); 
    CALC_STATE_ROUND2 (dwState1, dwState2, dwState3, dwState0, x[ 0], r[19], k[19]);
    CALC_STATE_ROUND2 (dwState0, dwState1, dwState2, dwState3, x[ 5], r[20], k[20]);
    CALC_STATE_ROUND2 (dwState3, dwState0, dwState1, dwState2, x[10], r[21], k[21]); 
    CALC_STATE_ROUND2 (dwState2, dwState3, dwState0, dwState1, x[15], r[22], k[22]); 
    CALC_STATE_ROUND2 (dwState1, dwState2, dwState3, dwState0, x[ 4], r[23], k[23]); 
    CALC_STATE_ROUND2 (dwState0, dwState1, dwState2, dwState3, x[ 9], r[24], k[24]); 
    CALC_STATE_ROUND2 (dwState3, dwState0, dwState1, dwState2, x[14], r[25], k[25]); 
    CALC_STATE_ROUND2 (dwState2, dwState3, dwState0, dwState1, x[ 3], r[26], k[26]); 
    CALC_STATE_ROUND2 (dwState1, dwState2, dwState3, dwState0, x[ 8], r[27], k[27]); 
    CALC_STATE_ROUND2 (dwState0, dwState1, dwState2, dwState3, x[13], r[28], k[28]); 
    CALC_STATE_ROUND2 (dwState3, dwState0, dwState1, dwState2, x[ 2], r[29], k[29]); 
    CALC_STATE_ROUND2 (dwState2, dwState3, dwState0, dwState1, x[ 7], r[30], k[30]); 
    CALC_STATE_ROUND2 (dwState1, dwState2, dwState3, dwState0, x[12], r[31], k[31]); 
    
    //3rd round 
    CALC_STATE_ROUND3 (dwState0, dwState1, dwState2, dwState3, x[ 5], r[32], k[32]); 
    CALC_STATE_ROUND3 (dwState3, dwState0, dwState1, dwState2, x[ 8], r[33], k[33]); 
    CALC_STATE_ROUND3 (dwState2, dwState3, dwState0, dwState1, x[11], r[34], k[34]); 
    CALC_STATE_ROUND3 (dwState1, dwState2, dwState3, dwState0, x[14], r[35], k[35]); 
    CALC_STATE_ROUND3 (dwState0, dwState1, dwState2, dwState3, x[ 1], r[36], k[36]); 
    CALC_STATE_ROUND3 (dwState3, dwState0, dwState1, dwState2, x[ 4], r[37], k[37]); 
    CALC_STATE_ROUND3 (dwState2, dwState3, dwState0, dwState1, x[ 7], r[38], k[38]); 
    CALC_STATE_ROUND3 (dwState1, dwState2, dwState3, dwState0, x[10], r[39], k[39]); 
    CALC_STATE_ROUND3 (dwState0, dwState1, dwState2, dwState3, x[13], r[40], k[40]); 
    CALC_STATE_ROUND3 (dwState3, dwState0, dwState1, dwState2, x[ 0], r[41], k[41]); 
    CALC_STATE_ROUND3 (dwState2, dwState3, dwState0, dwState1, x[ 3], r[42], k[42]); 
    CALC_STATE_ROUND3 (dwState1, dwState2, dwState3, dwState0, x[ 6], r[43], k[43]); 
    CALC_STATE_ROUND3 (dwState0, dwState1, dwState2, dwState3, x[ 9], r[44], k[44]); 
    CALC_STATE_ROUND3 (dwState3, dwState0, dwState1, dwState2, x[12], r[45], k[45]); 
    CALC_STATE_ROUND3 (dwState2, dwState3, dwState0, dwState1, x[15], r[46], k[46]); 
    CALC_STATE_ROUND3 (dwState1, dwState2, dwState3, dwState0, x[ 2], r[47], k[47]); 
    
    //4th round
    CALC_STATE_ROUND4 (dwState0, dwState1, dwState2, dwState3, x[ 0], r[48], k[48]); 
    CALC_STATE_ROUND4 (dwState3, dwState0, dwState1, dwState2, x[ 7], r[49], k[49]); 
    CALC_STATE_ROUND4 (dwState2, dwState3, dwState0, dwState1, x[14], r[50], k[50]); 
    CALC_STATE_ROUND4 (dwState1, dwState2, dwState3, dwState0, x[ 5], r[51], k[51]); 
    CALC_STATE_ROUND4 (dwState0, dwState1, dwState2, dwState3, x[12], r[52], k[52]); 
    CALC_STATE_ROUND4 (dwState3, dwState0, dwState1, dwState2, x[ 3], r[53], k[53]); 
    CALC_STATE_ROUND4 (dwState2, dwState3, dwState0, dwState1, x[10], r[54], k[54]); 
    CALC_STATE_ROUND4 (dwState1, dwState2, dwState3, dwState0, x[ 1], r[55], k[55]); 
    CALC_STATE_ROUND4 (dwState0, dwState1, dwState2, dwState3, x[ 8], r[56], k[56]); 
    CALC_STATE_ROUND4 (dwState3, dwState0, dwState1, dwState2, x[15], r[57], k[57]); 
    CALC_STATE_ROUND4 (dwState2, dwState3, dwState0, dwState1, x[ 6], r[58], k[58]);
    CALC_STATE_ROUND4 (dwState1, dwState2, dwState3, dwState0, x[13], r[59], k[59]); 
    CALC_STATE_ROUND4 (dwState0, dwState1, dwState2, dwState3, x[ 4], r[60], k[60]); 
    CALC_STATE_ROUND4 (dwState3, dwState0, dwState1, dwState2, x[11], r[61], k[61]); 
    CALC_STATE_ROUND4 (dwState2, dwState3, dwState0, dwState1, x[ 2], r[62], k[62]); 
    CALC_STATE_ROUND4 (dwState1, dwState2, dwState3, dwState0, x[ 9], r[63], k[63]); 

    //save calculated states
    pstStates->dwState0 += dwState0;
    pstStates->dwState1 += dwState1;
    pstStates->dwState2 += dwState2;
    pstStates->dwState3 += dwState3;

    }

