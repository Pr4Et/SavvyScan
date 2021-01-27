HEAD  SpcInitPCIBoards  (int16 *pnCount, int16 *pnPCIVersion);
HEAD  SpcInitBoard      (int16 nNr, int16 nTyp);
HEAD  SpcSetParam       (int16 nNr, int32 lReg, int32 lValue);
HEAD  SpcGetParam       (int16 nNr, int32 lReg, int32 *plValue);
HEAD  SpcGetData        (int16 nNr, int16 nCh, int32 lStart, int32 lLen, dataptr pvData);
HEAD  SpcSetData        (int16 nNr, int16 nCh, int32 lStart, int32 lLen, dataptr pvData);
HEAD  SpcGetVersionInfo (char *pszBuffer, int nBufferLen);

// these functions are only needed under 64 bit Linux to set FIFO adresses
HEAD  SpcSetAdr        (int16 nNr, int32 lReg, void*  pvAdr);
HEAD  SpcGetAdr        (int16 nNr, int32 lReg, void** ppvAdr);
