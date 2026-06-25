// The IMOD version of functions.cpp
#include "core_headers.h"

// Headers for time and processor count
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#include <Windows.h>
typedef BOOL (WINAPI *LPFN_GLPI)(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);
#else
#include <sys/time.h>
#endif
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#define B3DMIN(a,b) ((a) < (b) ? (a) : (b))
#define B3DMAX(a,b) ((a) > (b) ? (a) : (b))
#define B3DCLAMP(a,b,c) a = B3DMAX((b), B3DMIN((c), (a)))

// From original functions.cpp
// Return a vector with the rank of the elements of the input array
std::vector<size_t> rankSort(const std::vector<float>& v_temp) {
    std::vector<std::pair<float, size_t> > v_sort(v_temp.size());

    for (size_t i = 0U; i < v_sort.size(); ++i) {
        v_sort[i] = std::make_pair(v_temp[i], i);
    }

    sort(v_sort.begin(), v_sort.end());

    std::pair<double, size_t> rank;
    std::vector<size_t> result(v_temp.size());

    for (size_t i = 0U; i < v_sort.size(); ++i) {
        if (v_sort[i].first != rank.first) {
            rank = std::make_pair(v_sort[i].first, i);
        }
        result[v_sort[i].second] = rank.second;
    }
    return result;
}

// For setting alternate print function
static CharArgType sPrintFunc = NULL;

void internalSetPrintFunc(CharArgType func)
{
  sPrintFunc = func;
}

// Print a message 
void wxPrintf(const char *format, ...)
{
  char errorMess[512];
  va_list args;
  va_start(args, format);
  vsprintf(errorMess, format, args);
  if (sPrintFunc) {
    sPrintFunc(errorMess);
  } else {
    printf("%s", errorMess);
  }
  va_end(args);
}

/*
 * Renamed or static COPIES of functions in b3dutil.c, to avoid dependency on libcfshr
 * This whole mess is for OpenMP support
 */
#define CPUINFO_LINE 80
#define MAX_CPU_SOCKETS 64
static int fgetline(FILE *fp, char s[], int limit)
{
  int c, i, length;

  if (fp == NULL){
    return(-1);
  }

  if (limit < 3){
    return(-1);
  }
     
  for (i=0; ( ((c = getc(fp)) != EOF) && (i < (limit-1)) && (c != '\n') ); i++)
    s[i]=c;

  /* 1/25/12: Take off a return too! */
  if (i > 0 && s[i-1] == '\r')
    i--;

  /* A \n or EOF on the first character leaves i at 0, so there is nothing
     special to be handled about i being 1, 9/18/09 */
               
  s[i]='\0';
  length = i;

  if (c == EOF)
    return (-1 * (length + 2));
  else
    return (length);
}


/*!
 * Returns the number of physical processor cores in [physical] and the number of logical
 * processors in [logical].  The return value is 1 if valid information could not be 
 * obtained for both items.  It uses sysctlbyname on Mac, GetLogicalProcessorInformation
 * on Windows, and /proc/cpuinfo on Linux.
 */
static int numCoresAndLogicalProcs(int *physical, int *logical)
{
  int processorCoreCount = 0;
  int logicalProcessorCount = 0;

#ifdef __APPLE__
  int temp = 0;
  size_t lenPhys = sizeof(int);
#elif defined(_WIN32)
  LPFN_GLPI glpi;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr = NULL;
  DWORD returnLength = 0;
  DWORD byteOffset = 0;
  DWORD numBits = sizeof(ULONG_PTR) * 8;
  ULONG_PTR bitTest;
  DWORD i;
#else
  FILE *fp;
  unsigned char socketFlags[MAX_CPU_SOCKETS];
  char linebuf[CPUINFO_LINE];
  int err, len, curID, curCores;
  char *colon;
#endif

#ifdef __APPLE__
  if (!sysctlbyname("hw.physicalcpu" , &temp, &lenPhys, NULL, 0))
    processorCoreCount = temp;
  lenPhys = 4;
  if (!sysctlbyname("hw.logicalcpu" , &temp, &lenPhys, NULL, 0))
    logicalProcessorCount = temp;

#elif defined(_WIN32)
  /* This is adapted from https://msdn.microsoft.com/en-us/library/ms683194 
     On systems with more than 64 logical processors, the GetLogicalProcessorInformation
     function retrieves logical processor information about processors in the processor
     group to which the calling thread is currently assigned. Use the 
     GetLogicalProcessorInformationEx function to retrieve information about processors
     in all processor groups on the system.  (Windows 7/Server 2008 or above).
  */
  glpi = (LPFN_GLPI)GetProcAddress(GetModuleHandle(TEXT("kernel32")),
                                   "GetLogicalProcessorInformation");
  if (glpi) {
    if (!glpi(buffer, &returnLength) && GetLastError() == ERROR_INSUFFICIENT_BUFFER &&
        returnLength > 0) {
      buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)malloc(returnLength);
      if (buffer) {
        if (glpi(buffer, &returnLength)) {
          ptr = buffer;
          while (byteOffset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <=
                 returnLength) {
            if (ptr->Relationship == RelationProcessorCore) {
              processorCoreCount++;
              bitTest = (ULONG_PTR)1;
              for (i = 0; i < numBits; i++) {
                if (ptr->ProcessorMask & bitTest)
                  logicalProcessorCount++;
                bitTest *= 2;
              }
            }              
            byteOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
            ptr++;
          }
          free(buffer);
        }
      }
    }
  }
#else

  /* Linux: look at /proc/cpuinfo */
  fp = fopen("/proc/cpuinfo", "r");
  if (fp) {
    curID = -1;
    curCores = -1;
    memset(socketFlags, 0, MAX_CPU_SOCKETS);
    while (1) {
      err = 0;
      len = fgetline(fp, linebuf, CPUINFO_LINE);
      if (!len)
        continue;
      if (len == -2)
        break;
      err = 1;
      if (len == -1)
        break;

      /* Look for a "physical id :" and a "cpu cores :" in either order */
      if (strstr(linebuf, "physical id")) {

        /* Error if already got a physical id without cpu cores */
        if (curID >= 0)
          break;
        colon = strchr(linebuf, ':');
        if (colon)
          curID = atoi(colon+1);

        /* Error if no colon or ID out of range */
        if (!colon || curID < 0 || curID >= MAX_CPU_SOCKETS)
          break;
      }
      if (strstr(linebuf, "cpu cores")) {

        /* Error if already got a cpu cores without physical id  */
        if (curCores >= 0)
          break;
        colon = strchr(linebuf, ':');
        if (colon)
          curCores = atoi(colon+1);

        /* Error if no colon or core count illegal */
        if (!colon || curCores <= 0)
          break;
      }

      /* If have both ID and core count, add one logical processor and the number of
         cores the first time this ID is seen to the core count; set ID flag and reset
         the tow numbers */
      if (curID >= 0 && curCores > 0) {
        logicalProcessorCount++;
        if (!socketFlags[curID])
          processorCoreCount += curCores;
        socketFlags[curID] = 1;
        curID = -1;
        curCores = -1;
      }
      err = 0;
      if (len < 0)
        break;
    }
    if (err)
      processorCoreCount *= -1;
    fclose(fp);
  }
#endif
  *physical = processorCoreCount;
  *logical = logicalProcessorCount;
  return (processorCoreCount <= 0 || logicalProcessorCount < 0) ? 1 : 0;
}

// Get an appropriately limited number of threads
int ctfNumOMPthreads(int optimalThreads)
{
  int numThreads = optimalThreads;
  int physicalProcs = 0;
  int logicalProcessorCount = 0;
  int processorCoreCount = 0;
#ifdef _OPENMP
  static int limThreads = -1;
  static int numProcs = -1;
  static int forceThreads = -1;
  static int ompNumProcs = -1;
  char *ompNum;

  /* One-time determination of number of physical and logical cores */
  if (numProcs < 0) {
    ompNumProcs = numProcs = omp_get_num_procs();

    /* if there are legal numbers and the logical count is the OMP
     number, set the physical processor count */
    if (!numCoresAndLogicalProcs(&processorCoreCount, &logicalProcessorCount) &&
        processorCoreCount > 0 && logicalProcessorCount == numProcs)
      physicalProcs = processorCoreCount;
    if (getenv("IMOD_REPORT_CORES"))
      printf("core count = %d  logical processors = %d  OMP num = %d => physical "
             "processors = %d\n", processorCoreCount, logicalProcessorCount, numProcs,
             physicalProcs); fflush(stdout);
      
    if (physicalProcs > 0)
      numProcs = B3DMIN(numProcs, physicalProcs);
  }

  /* Limit by number of real cores */
  numThreads = B3DMAX(1, B3DMIN(numProcs, numThreads));

  /* One-time determination of the limit set by OMP_NUM_THREADS */
  if (limThreads < 0) {
    ompNum = getenv("OMP_NUM_THREADS");
    if (ompNum)
      limThreads = atoi(ompNum);
    limThreads = B3DMAX(0, limThreads);
  }

  /* Limit to number set by OMP_NUM_THREADS and to number of real cores */
  if (limThreads > 0)
    numThreads = B3DMIN(limThreads, numThreads);

  /* One-time determination of whether user wants to force a number of threads */
  if (forceThreads < 0) {
    forceThreads = 0;
    ompNum = getenv("IMOD_FORCE_OMP_THREADS");
    if (ompNum) {
      if (!strcmp(ompNum, "ALL_CORES")) {
        if (numProcs > 0)
          forceThreads = numProcs;
      } else if (!strcmp(ompNum, "ALL_HYPER")) {
        if (ompNumProcs > 0)
          forceThreads = ompNumProcs;
      } else {
        forceThreads = atoi(ompNum);
        forceThreads = B3DMAX(0, forceThreads);
      }
    }
  }

  /* Force the number if set */
  if (forceThreads > 0)
    numThreads = forceThreads;

  if (getenv("IMOD_REPORT_CORES"))
      printf("numProcs %d  limThreads %d  numThreads %d\n", numProcs,
             limThreads, numThreads); fflush(stdout);
#else
  numThreads = 1;
#endif
  return numThreads;
}

/*!
 * Returns the thread number of the current thread, numbered from 0
 */
int ctfOMPthreadNum()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/*!
 * Returns a measure of time in seconds with microsecond precision on Linux and Mac
 * and the precision of the high performance counter on Windows.
 */
double ctfWallTime(void)
{
#ifdef _WIN32
  LARGE_INTEGER freq, counts;
  QueryPerformanceFrequency(&freq);
  if (!freq.QuadPart)
    return 0.;
  QueryPerformanceCounter(&counts);
  return (((double)counts.QuadPart) / freq.QuadPart);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((double)tv.tv_sec + tv.tv_usec / 1000000.);
#endif
}
