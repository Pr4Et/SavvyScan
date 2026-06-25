/*
**************************************************************************

spcm_lib_thread.cpp                                      (c) Spectrum GmbH

**************************************************************************

Offers data transfer threads for sample data, timestamp data and ABA date.
Each thread get's a bundle of parameters including a work function that
is used to do something with the data. These threads are used in all our
FIFO examples.

**************************************************************************
*/


// ----- include of common example librarys -----
#include "../common/spcm_lib_thread.h"

// ----- standard c include files -----
#include <stdio.h>
#include <string.h>



/*
**************************************************************************
bInitThreadStructure: initializes all mutexes and events. If the parameter
hSharedUpdate is set we do not allocate a new event for updat ebut use 
this one as a shared event
**************************************************************************
*/

bool bInitThreadStructure (ST_THREADDATA* pstThreadData, SPCM_EVENT_HANDLE* phSharedEventUpdate)
    {
    bool bSetupOk = true;

    // initialize the structre with 0
    memset (pstThreadData, 0, sizeof (ST_THREADDATA));

    // init the mutex and the events
    bSetupOk &= (bSetupOk && spcm_bCreateMutex (&pstThreadData->hMutexSharedAccess));
    bSetupOk &= (bSetupOk && spcm_bCreateEvent (&pstThreadData->hEventInitDone));
    bSetupOk &= (bSetupOk && spcm_bCreateEvent (&pstThreadData->hEventStart));
    bSetupOk &= (bSetupOk && spcm_bCreateEvent (&pstThreadData->hEventStartDone));
    bSetupOk &= (bSetupOk && spcm_bCreateEvent (&pstThreadData->hEventTermDone));

    // if a shared event is given we store this and set a matching flag
    if (phSharedEventUpdate != NULL)
        {
        pstThreadData->phSharedUpdate = phSharedEventUpdate;
        pstThreadData->bSharedUpdate =  true;
        }

    // otherwise we allocate a new event
    else
        {
        bSetupOk &= (bSetupOk && spcm_bCreateEvent (&pstThreadData->hEventUpdate));
        pstThreadData->phSharedUpdate = &pstThreadData->hEventUpdate;
        pstThreadData->bSharedUpdate = false;
        }

    return bSetupOk;
    }



/*
**************************************************************************
vCleanThreadStructure: clean up (we don't need to care for correct 
initialisation as this is checked in the close routines)
**************************************************************************
*/

void vCleanThreadStructure (ST_THREADDATA* pstThreadData)
    {
    spcm_vCloseMutex (&pstThreadData->hMutexSharedAccess);
    spcm_vCloseEvent (&pstThreadData->hEventStart);
    spcm_vCloseEvent (&pstThreadData->hEventInitDone);
    spcm_vCloseEvent (&pstThreadData->hEventStartDone);
    spcm_vCloseEvent (&pstThreadData->hEventTermDone);

    if (!pstThreadData->bSharedUpdate)
        spcm_vCloseEvent (&pstThreadData->hEventUpdate);
    }




/*
**************************************************************************
main loop. the loop cares for card, data, timestamp and aba data all 
enabled by flags in the pstBufferData structure. 

This loop can be either called directly or from within the thread function.
if it's called from the thread function it will handle the mutexes and
events also.
  
Mainly the loop does the following steps:

  Call bWorkSetup (user function to setup all details)
  Allocate memory for buffer
  Start card and transfers depending on the flags that have been set
  Loop
    Wait for interrupt that announces new data or timeout to occur
    Read Status and current positions in the buffers
    Recalculate positions and length
    Call bWorkDo (user function) with the available data
    Set the data that was processed available to card again
  Until (Abort OR Error)
  Call vWorkClose (user function to close files and clean up)
  Clean up buffers

**************************************************************************
*/


void vMainLoop (
    ST_BUFFERDATA*      pstBufferData,          // filled card info structure
    void*               pvWorkData,             // working data structure that is passed to the calback functions
    SPCM_WORK_INIT*     pfn_bWorkSetup,         // callback function for pre-run setup
    SPCM_WORK_DO*       pfn_bWorkDo,            // callback function that is called if new data is available
    SPCM_WORK_CLOSE*    pfn_vWorkClose,         // callback function at the end of run
    SPCM_ABORTCHECK*    pfn_bAbortCheck,        // callback function for abort checking on each update loop
    ST_THREADDATA*      pstThreadData)          // threaddata for synchronization. Can be NULL if non threaded version is used

// ***********************************************************************

    {
    void*   pvDataBuffer =  NULL;
    void*   pvTSBuffer =    NULL;
    void*   pvABABuffer =   NULL;
    int32   lStartCmd =     0;
    int32   lWaitCmd =      0;
    int32   lStatus =       0;
    int64   llAvailPos =     0;
    uint32  dwError =       ERR_OK;
    bool    bContMemUsed =  false;


    // secure the initialisation
    if (pstThreadData)
        spcm_vGetMutex (&pstThreadData->hMutexSharedAccess);

    // first we do the setup
    pstBufferData->qwDataTransferred =  0;
    pstBufferData->qwTSTransferred =    0;
    pstBufferData->qwABATransferred =   0;
    pstBufferData->dwError =            ERR_OK;
    if (pfn_bWorkSetup)
        if (!(*pfn_bWorkSetup) (pvWorkData, pstBufferData))
            {
            printf ("Setup failed!\n");
            dwError = ERR_INIT;
            }



    // now the setup routine has defined the buffer details and we allocate data
    if (!dwError)
        {
        if (pstBufferData->bStartData)
            {

            // try to use cont mem buffer if size matches
            uint64 qwContBufLen = 0;
            spcm_dwGetContBuf_i64 (pstBufferData->pstCard->hDrv, SPCM_BUF_DATA, &pvDataBuffer, &qwContBufLen);
            bContMemUsed = (qwContBufLen >= pstBufferData->dwDataBufLen);
            if (!bContMemUsed)
                pvDataBuffer = pvAllocMemPageAligned (pstBufferData->dwDataBufLen);
            }
        if (pstBufferData->bStartTimestamp)
            pvTSBuffer = pvAllocMemPageAligned (pstBufferData->dwTSBufLen);
        if (pstBufferData->bStartABA)
            pvABABuffer = pvAllocMemPageAligned (pstBufferData->dwABABufLen);
        if ((!pvDataBuffer && pstBufferData->bStartData) || (!pvTSBuffer && pstBufferData->bStartTimestamp) || (!pvABABuffer && pstBufferData->bStartABA))
            {
            printf ("Memory Allocation Error\n");
            dwError = ERR_MEMALLOC;
            }
        }


    // all is prepared and we define the buffers for the transfer
    if (!pstBufferData->dwError && pstBufferData->bStartData)
        pstBufferData->dwError = spcm_dwDefTransfer_i64 (pstBufferData->pstCard->hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, (uint32) pstBufferData->dwDataNotify, pvDataBuffer, 0, pstBufferData->dwDataBufLen);
    if (!pstBufferData->dwError && pstBufferData->bStartTimestamp)
        pstBufferData->dwError = spcm_dwDefTransfer_i64 (pstBufferData->pstCard->hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, (uint32) pstBufferData->dwTSNotify, pvTSBuffer, 0, pstBufferData->dwTSBufLen);
    if (!pstBufferData->dwError && pstBufferData->bStartABA)
        pstBufferData->dwError = spcm_dwDefTransfer_i64 (pstBufferData->pstCard->hDrv, SPCM_BUF_ABA, SPCM_DIR_CARDTOPC, (uint32) pstBufferData->dwABANotify, pvABABuffer, 0, pstBufferData->dwABABufLen);


    // send a message that the initialisation has been done
    if (pstThreadData)
        {
        pstThreadData->dwEventFlags |= EVENTFLAG_INITDONE;
        spcm_vSignalEvent (&pstThreadData->hEventInitDone);
        }


    // start of the card together with all transfers if enabled
    if (!dwError)
        {

        // wait for the start event from the main thread
        if (pstThreadData)
            {
            if ((pstThreadData->dwEventFlags & EVENTFLAG_START) == 0)
                spcm_bWaitEventWithMutex (&pstThreadData->hEventStart, &pstThreadData->hMutexSharedAccess);
            spcm_vReleaseMutex (&pstThreadData->hMutexSharedAccess);
            }

        if (pstBufferData->bStartCard)
            lStartCmd = M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER;

        if (pstBufferData->bStartData)
            {
            lStartCmd |= M2CMD_DATA_STARTDMA;
            lWaitCmd |=  M2CMD_DATA_WAITDMA;
            }

        if (pstBufferData->bStartExtraDMA)
            lStartCmd |= M2CMD_EXTRA_STARTDMA;

        if (pstBufferData->bStartTimestamp || pstBufferData->bStartABA)
            lWaitCmd |=  M2CMD_EXTRA_WAITDMA;

        spcm_dwSetParam_i32 (pstBufferData->pstCard->hDrv, SPC_TIMEOUT, pstBufferData->lTimeout); 
        dwError = spcm_dwSetParam_i32 (pstBufferData->pstCard->hDrv, SPC_M2CMD, lStartCmd);

        if (pstThreadData)
            {
            pstThreadData->dwEventFlags |= EVENTFLAG_STARTDONE;
            spcm_vSignalEvent (&pstThreadData->hEventStartDone);
            }
        }



    // this is the main loop that waits for data. As we defined a timeout this routine is latest called in our timeout interval
    while (!dwError)
        {
        dwError = spcm_dwSetParam_i32 (pstBufferData->pstCard->hDrv, SPC_M2CMD, lWaitCmd);
        spcm_dwGetParam_i32 (pstBufferData->pstCard->hDrv, SPC_M2STATUS, &lStatus);
        
        // Recording complete (for setups with SPC_LOOPS != 0)?
        if (lStatus & M2STAT_DATA_END)
            {
            // DATA_END can also occur if all data has been transfered from card after overrun
            if (lStatus & M2STAT_DATA_OVERRUN)
                dwError = ERR_FIFOHWOVERRUN;
            else
                if (pstBufferData->llDataAvailBytes == 0)
                    dwError = ERR_FIFOFINISHED;
            }

        // secure the following access as it's the main data loop
        if (pstThreadData)
            spcm_vGetMutex (&pstThreadData->hMutexSharedAccess);


        // read out the current position of data buffer and recalculate it to avoid rollover
        if (pstBufferData->bStartData)
            {
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_DATA_AVAIL_USER_LEN,   &pstBufferData->llDataAvailBytes);
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_DATA_AVAIL_USER_POS,   &llAvailPos);
            if ((llAvailPos + pstBufferData->llDataAvailBytes) >= pstBufferData->dwDataBufLen)
                pstBufferData->llDataAvailBytes = (uint64) (pstBufferData->dwDataBufLen - llAvailPos);
            pstBufferData->pvDataCurrentBuf = (void*) (((char*) pvDataBuffer) + llAvailPos);
            }


        // read out the timestanp positions and recalculate it to avoid rollover
        if (pstBufferData->bStartTimestamp)
            {
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_TS_AVAIL_USER_LEN,     &pstBufferData->llTSAvailBytes);
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_TS_AVAIL_USER_POS,     &llAvailPos);
            if ((llAvailPos + pstBufferData->llTSAvailBytes) >= pstBufferData->dwTSBufLen)
                pstBufferData->llTSAvailBytes = pstBufferData->dwTSBufLen - llAvailPos;
            pstBufferData->pvTSCurrentBuf = (void*) (((char*) pvTSBuffer) + llAvailPos);
            }


        // read out the ABA positions and recalculate it to avoid rollover
        if (pstBufferData->bStartABA)
            {
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_ABA_AVAIL_USER_LEN,     &pstBufferData->llABAAvailBytes);
            spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_ABA_AVAIL_USER_POS,     &llAvailPos);
            if ((llAvailPos + pstBufferData->llABAAvailBytes) >= pstBufferData->dwABABufLen)
                pstBufferData->llABAAvailBytes = pstBufferData->dwABABufLen - llAvailPos;
            pstBufferData->pvABACurrentBuf = (void*) (((char*) pvABABuffer) + llAvailPos);
            }

        // if there is some new data we call our work function
        if (pstBufferData->llDataAvailBytes || pstBufferData->llTSAvailBytes || pstBufferData->llABAAvailBytes)
            if (pfn_bWorkDo)
                if (!(*pfn_bWorkDo) (pvWorkData, pstBufferData))
                    dwError = ERR_ABORT;


        // now we have to set the data available for the card again (the number of bytes may have been modified by the work routine)
        if (pstBufferData->bStartData && pstBufferData->llDataAvailBytes)
            spcm_dwSetParam_i64 (pstBufferData->pstCard->hDrv, SPC_DATA_AVAIL_CARD_LEN, pstBufferData->llDataAvailBytes);

        if (pstBufferData->bStartTimestamp && pstBufferData->llTSAvailBytes)
            spcm_dwSetParam_i64 (pstBufferData->pstCard->hDrv, SPC_TS_AVAIL_CARD_LEN, pstBufferData->llTSAvailBytes);

        if (pstBufferData->bStartABA && pstBufferData->llABAAvailBytes)
            spcm_dwSetParam_i64 (pstBufferData->pstCard->hDrv, SPC_ABA_AVAIL_CARD_LEN, pstBufferData->llABAAvailBytes);
       

        // we continuously count the number of bytes that we have transferred so far
        pstBufferData->qwDataTransferred +=   pstBufferData->llDataAvailBytes;
        pstBufferData->qwTSTransferred +=     pstBufferData->llTSAvailBytes;
        pstBufferData->qwABATransferred +=    pstBufferData->llABAAvailBytes;


        // at last we check for abort (not in the threaded version as this is then done in the main thread)
        if (pfn_bAbortCheck && !pstThreadData)
            if ((*pfn_bAbortCheck) (pvWorkData, pstBufferData))
               dwError = ERR_ABORT;


        // signal the new status to the main thread
        if (pstThreadData)
            {
            pstBufferData->dwError = dwError;
            spcm_vSignalEvent (pstThreadData->phSharedUpdate);
            spcm_vReleaseMutex (&pstThreadData->hMutexSharedAccess);
            }

        // we don't count timeout as an error as we've forced the loop with a small timeout (non-threaded version)
        if (!pstThreadData)
            if (dwError == ERR_TIMEOUT)
                dwError = ERR_OK;

        // we have to suspend this thread for a short tiem to give other threads some working time
        if (pstThreadData)
            spcm_vSuspendThread (1);
        } // while (!dwError) (Main Loop)


    // check for the error in the loop
    printf ("\n\n");
    switch (dwError)
        {
        case ERR_FIFOHWOVERRUN: printf ("Loop: Hardware Overrun\n"); break;
        case ERR_FIFOFINISHED:  printf ("Loop: FIFO Mode finished"); break;
        case ERR_OK:            printf ("Loop: Finished normal\n"); break;
        case ERR_ABORT:         printf ("Loop: Aborted...\n"); break;
        default:                printf ("Loop: Error %4d\n", dwError); break;
        }



    // copy the error information to the buffer structure
    if (pstThreadData)
        spcm_vGetMutex (&pstThreadData->hMutexSharedAccess);

    pstBufferData->dwError = dwError;

    // at least we do one last upodate with the final error code
    if (pstThreadData)
        {
        spcm_vReleaseMutex (&pstThreadData->hMutexSharedAccess);
        spcm_vSignalEvent (pstThreadData->phSharedUpdate);
        }



    // close the work
    if (pfn_vWorkClose)
        (*pfn_vWorkClose) (pvWorkData, pstBufferData);

    // clean up
    if (pvDataBuffer && !bContMemUsed)
        vFreeMemPageAligned (pvDataBuffer, pstBufferData->dwDataBufLen);
    if (pvTSBuffer)
        vFreeMemPageAligned (pvTSBuffer, pstBufferData->dwTSBufLen);
    if (pvABABuffer)
        vFreeMemPageAligned (pvABABuffer, pstBufferData->dwABABufLen);

    // we need to leave the thread with a released mutex, otherwise the other threads may block
    if (pstThreadData)
        {
        pstThreadData->dwEventFlags |= EVENTFLAG_TERMDONE;
        spcm_vSignalEvent (&pstThreadData->hEventTermDone);
        }
    }




/*
**************************************************************************
vDoMainLoop: 

the non threaded version of the main loop, just calls the shared main loop 
and let the thread data section be NULL to disable all synchronization
**************************************************************************
*/

void vDoMainLoop (
    ST_BUFFERDATA*      pstBufferData,          // filled card info structure
    void*               pvWorkData,             // working data structure that is passed to the calback functions
    SPCM_WORK_INIT*     pfn_bWorkSetup,         // callback function for pre-run setup
    SPCM_WORK_DO*       pfn_bWorkDo,            // callback function that is called if new data is available
    SPCM_WORK_CLOSE*    pfn_vWorkClose,         // callback function at the end of run
    SPCM_ABORTCHECK*    pfn_bAbortCheck)        // callback function for abort checking on each update loop

    {
    vMainLoop (
        pstBufferData, 
        pvWorkData, 
        pfn_bWorkSetup, 
        pfn_bWorkDo, 
        pfn_vWorkClose, 
        pfn_bAbortCheck, 
        NULL);
    }



/*
**************************************************************************

DataTransfer Thread: it uses the main loop to run the specified thread
function. 

**************************************************************************
*/

SPCM_THREAD_RETURN SPCM_THREAD_CALLTYPE pvDataTransferThread (void* pvArguments)
    {
    ST_THREADDATA* pstThreadData = (ST_THREADDATA*) pvArguments;
    ST_BUFFERDATA* pstBufferData = (ST_BUFFERDATA*) pstThreadData->pstBufferData;
    void*          pvWorkData =    pstThreadData->pvWorkData;

    vMainLoop (
        pstBufferData, 
        pvWorkData, 
        pstThreadData->pfn_bWorkInit, 
        pstThreadData->pfn_bWorkDo,
        pstThreadData->pfn_vWorkClose,
        NULL,
        pstThreadData);

    return 0;
    }




/*
**************************************************************************
vDoThreadMainLoop: contains the main FIFO loop in the threaded version.
This function mainly initializes the thread structures and calls
the threads to do the work
**************************************************************************
*/

void vDoThreadMainLoop (
    ST_BUFFERDATA*      pstBufferData,          // filled card info structure
    void*               pvWorkData,             // working data structure that is passed to the calback functions

    SPCM_WORK_INIT*     pfn_bWorkInit,          // callback function for pre-run setup
    SPCM_WORK_DO*       pfn_bWorkDo,            // callback function that is called if new data is available
    SPCM_WORK_CLOSE*    pfn_vWorkClose,         // callback function at the end of run

    SPCM_ABORTCHECK*    pfn_bAbortCheck,        // callback function for abort checking on each update loop
    bool                bSeparateThread,        // use separate threads for data, timestamp and ABA

    SPCM_WORK_INIT*     pfn_bTSInit,            // callback function for timestamo pre-run setup
    SPCM_WORK_DO*       pfn_bTSDo,              // callback function that is called if new timestamps are available
    SPCM_WORK_CLOSE*    pfn_vTSClose,           // callback function at the end of timestamp run

    SPCM_WORK_INIT*     pfn_bABAInit,           // callback function for ABA pre-run setup
    SPCM_WORK_DO*       pfn_bABADo,             // callback function that is called if new ABA data is available
    SPCM_WORK_CLOSE*    pfn_vABAClose)          // callback function at the end of ABA run

// ***********************************************************************

    {
    SPCM_THREAD_HANDLE  hThread[3];             // handles for the maximum of three threads we support
    int16               nThreadCount = 0;       // number of threads we currently use
    ST_THREADDATA       stThreadData[3];        // shared data structure between main and data thread
    ST_BUFFERDATA       stBufferData[3];        // own copies of buffer data
    bool                bOk = true;
    int                 i;
    uint32              dwError;                // error code of the thread

    // initialize all to zero for later clean up
    for (i=0; i<3; i++)
        {
        hThread[i] = NULL_HANDLE;
        memset (&stThreadData[i], 0, sizeof(ST_THREADDATA));
        }


    // we set the private buffer data as a copy of the given one and disable the starting flags for separate transfers
    for (i=0; i<3; i++)
        {
        memcpy ((void*) &stBufferData[i], pstBufferData, sizeof (ST_BUFFERDATA));

        // the transfers should be started separately therefore we need to set the flags individually
        if (bSeparateThread)
            {
            stBufferData[i].bStartABA =         false;
            stBufferData[i].bStartTimestamp =   false;
            stBufferData[i].bStartData =        false;
            stBufferData[i].bStartCard =        false;
            }
        }


    // count the threads to use and set separate start flags and separate work callback functions
    if (bSeparateThread)
        {
        if (pstBufferData->bStartData)
            {
            stBufferData[nThreadCount].bStartData =         true;
            stBufferData[nThreadCount].bStartExtraDMA =     false;
            stBufferData[nThreadCount].bStartCard =         pstBufferData->bStartCard;
            nThreadCount++;
            }

        if (pstBufferData->bStartTimestamp)
            {
            stBufferData[nThreadCount].bStartTimestamp =    true;
            stBufferData[nThreadCount].bStartExtraDMA =     true;
            nThreadCount++;
            }

        // extra dma can only be started once.
        if (pstBufferData->bStartABA)
            {
            stBufferData[nThreadCount].bStartABA =          true;
            stBufferData[nThreadCount].bStartExtraDMA =     !pstBufferData->bStartTimestamp;
            nThreadCount++;
            }
        }
    else
        {
        stBufferData[nThreadCount].bStartExtraDMA = pstBufferData->bStartTimestamp || pstBufferData->bStartABA;
        nThreadCount = 1;
        }




    // setup the thread data structure we thereby use a shared update as this main thread should get all update calls
    // the work data structure is also shared between the different threads
    for (i=0; (i<nThreadCount) && bOk; i++)
        {
        bOk = bInitThreadStructure (&stThreadData[i], (i == 0 ? NULL : &stThreadData[0].hEventUpdate));
        stThreadData[i].pstBufferData =     &stBufferData[i];
        stThreadData[i].pvWorkData =        pvWorkData;

        if (stBufferData[i].bStartData)
            {
            stThreadData[i].pfn_bWorkDo =        pfn_bWorkDo;
            stThreadData[i].pfn_bWorkInit =      pfn_bWorkInit;
            stThreadData[i].pfn_vWorkClose =     pfn_vWorkClose;
            }

        if (stBufferData[i].bStartTimestamp && bSeparateThread)
            {
            stThreadData[i].pfn_bWorkDo =        pfn_bTSDo;
            stThreadData[i].pfn_bWorkInit =      pfn_bTSInit;
            stThreadData[i].pfn_vWorkClose =     pfn_vTSClose;
            }

        if (stBufferData[i].bStartABA && bSeparateThread)
            {
            stThreadData[i].pfn_bWorkDo =        pfn_bABADo;
            stThreadData[i].pfn_bWorkInit =      pfn_bABAInit;
            stThreadData[i].pfn_vWorkClose =     pfn_vABAClose;
            }

        if (!stThreadData[i].pfn_bWorkDo)
            stThreadData[i].pfn_bWorkDo =        pfn_bWorkDo;

        if (!stThreadData[i].pfn_bWorkInit)
            stThreadData[i].pfn_bWorkInit =      pfn_bWorkInit;


        if (!stThreadData[i].pfn_vWorkClose)
            stThreadData[i].pfn_vWorkClose =     pfn_vWorkClose;
        }

    if (!bOk)
        {
        printf ("An error occurred while initializing mutex and events!\n");
        return;
        }



    // now we start the threads
    for (i=0; (i<nThreadCount) && bOk; i++)
        bOk = spcm_bCreateThread (&pvDataTransferThread, &hThread[i], (void*) &stThreadData[i]);

    if (!bOk)
        printf ("An error occurred while initializing mutex and events!\n");

    // no error, we wait for initialization and check for error
    if (bOk)
        printf ("\nThread(s) are started now\n...wait for initialisation of thread\n");

    // wait for thread initialisation done
    for (i=0; (i<nThreadCount) && bOk; i++)
        {
        spcm_vGetMutex (&stThreadData[i].hMutexSharedAccess);
        if ((stThreadData[i].dwEventFlags & EVENTFLAG_INITDONE) == 0)
            spcm_bWaitEventWithMutex (&stThreadData[i].hEventInitDone, &stThreadData[i].hMutexSharedAccess);
        if (stThreadData[i].pstBufferData->dwError)
            {
            bOk = false;
            printf ("...initialisation failed with error %d\n", stThreadData[i].pstBufferData->dwError);
            dwError = stThreadData[i].pstBufferData->dwError;
            }
        }


    // initialisation ok, we now send the start event
    if (bOk)
        printf ("...initialisation done, thread is now started\n");

    // first step: all threads but not the data (and card start) thread are started
    for (i=0; (i<nThreadCount) && bOk; i++)
        if (!stThreadData[i].pstBufferData->bStartData)
            {
            stThreadData[i].dwEventFlags |= EVENTFLAG_START;
            spcm_vSignalEvent (&stThreadData[i].hEventStart);

            if ((stThreadData[i].dwEventFlags & EVENTFLAG_STARTDONE) == 0)
                spcm_bWaitEventWithMutex (&stThreadData[i].hEventStartDone, &stThreadData[i].hMutexSharedAccess);

            spcm_vReleaseMutex (&stThreadData[i].hMutexSharedAccess);
            }

    // second step: the data thread is started and with this the card is started
    for (i=0; (i<nThreadCount) && bOk; i++)
        if (stThreadData[i].pstBufferData->bStartData)
            {
            stThreadData[i].dwEventFlags |= EVENTFLAG_START;
            spcm_vSignalEvent (&stThreadData[i].hEventStart);

            if ((stThreadData[i].dwEventFlags & EVENTFLAG_STARTDONE) == 0)
                spcm_bWaitEventWithMutex (&stThreadData[i].hEventStartDone, &stThreadData[i].hMutexSharedAccess);

            spcm_vReleaseMutex (&stThreadData[i].hMutexSharedAccess);
            }


    // we now simply wait for the (shared) update command
    while (bOk)
        {
        spcm_vWaitEvent (&stThreadData[0].hEventUpdate);

        // on each update we can abort the loop
        if (pfn_bAbortCheck)
            if (pfn_bAbortCheck (pvWorkData, &stBufferData[0]))
                {
                spcm_dwSetParam_i32 (pstBufferData->pstCard->hDrv, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA | M2CMD_EXTRA_STOPDMA);
                printf ("\n...Aborted by user\n");
                break;
                }

        // we check all error codes of running threads
        for (i=0; (i<nThreadCount) && bOk; i++)
            {
            spcm_vGetMutex (&stThreadData[i].hMutexSharedAccess);
            dwError = stThreadData[i].pstBufferData->dwError;
            bOk = (dwError == ERR_OK);
            spcm_vReleaseMutex (&stThreadData[i].hMutexSharedAccess);
            }
        }

    // print the error code
    switch (dwError)
        {
        case ERR_TIMEOUT:       printf ("\nTimeout has occured"); break;
        case ERR_FIFOHWOVERRUN: printf ("\nHardware Overrun"); break;
        case ERR_FIFOFINISHED:  printf ("\nFIFO Mode finished"); break;
        case ERR_OK:            printf ("\nFinished normal"); break;
        case ERR_ABORT:         printf ("\nAborted by user"); break;
        case ERR_MEMALLOC:      printf ("\nMemory Allocation error"); break;
        default:                printf ("\nError %4d ", dwError); break;
        }
    
    pstBufferData->dwError = dwError;

    // we have to stop the transfer to force the threads to return
    spcm_dwSetParam_i32 (pstBufferData->pstCard->hDrv, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA | M2CMD_EXTRA_STOPDMA);

    // now it's time to clean up
    for (i=0; i<nThreadCount; i++)
        {
        if ((stThreadData[i].dwEventFlags & EVENTFLAG_TERMDONE) == 0)
            spcm_vWaitEvent (&stThreadData[i].hEventTermDone);
        spcm_vJoinThread (&hThread[i], 0);
        spcm_vCloseThread (&hThread[i]);
        vCleanThreadStructure (&stThreadData[i]);
        }
    }



/*
**************************************************************************
bKeyAbortCheck: simple abort function that check for escape
**************************************************************************
*/ 

bool bKeyAbortCheck (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    if (bKbhit())
        if (cGetch() == 27)
            return true;

    return false;
    }
