/*  winversion.c - Simple functions to test for a windows version
 *
 * $Id$
 */
#ifdef _WIN32
#include <Windows.h>
#endif

/* All functions return 1 if true, 0 if Windows and not true, or -1 if not Windows */
int isWindows2000()
{
#ifdef _WIN32
  return isWindowsVersion(5, VER_EQUAL, 0, VER_EQUAL);
#else
  return -1;
#endif
}

int isWindowsXP()
{
#ifdef _WIN32
  return isWindowsVersion(5, VER_EQUAL, 1, VER_GREATER_EQUAL);
#else
  return -1;
#endif
}

int isWindowsVista()
{
#ifdef _WIN32
  return isWindowsVersion(6, VER_EQUAL, 0, VER_EQUAL);
#else
  return -1;
#endif
}

int isWindows7()
{
#ifdef _WIN32
  return isWindowsVersion(6, VER_EQUAL, 1, VER_EQUAL);
#else
  return -1;
#endif
}

/* This could be broken into 8.0 (6.2) and 8.1 (6.3+) */
int isWindows8()
{
#ifdef _WIN32
  return isWindowsVersion(6, VER_EQUAL, 2, VER_GREATER_EQUAL);
#else
  return -1;
#endif
}

/* Actually Windows 10 or higher */
int isWindows10()
{
#ifdef _WIN32
  return isWindowsVersion(10, VER_EQUAL, 0, VER_GREATER_EQUAL);
#else
  return -1;
#endif
}

/* The real test function, must be called with the VER_EQUAL, etc op values */
int isWindowsVersion(int major, int opMajor, int minor, int opMinor)
{
#ifdef _WIN32  
   OSVERSIONINFOEX osvi;
   DWORDLONG dwlConditionMask = 0;
   int op = VER_EQUAL;

   // Initialize the OSVERSIONINFOEX structure.
   ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
   osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
   osvi.dwMajorVersion = major;
   osvi.dwMinorVersion = minor;
   osvi.dwPlatformId = 2;

   // Initialize the condition mask.
   VER_SET_CONDITION( dwlConditionMask, VER_MAJORVERSION, opMajor);
   VER_SET_CONDITION( dwlConditionMask, VER_MINORVERSION, opMinor);
   VER_SET_CONDITION( dwlConditionMask, VER_PLATFORMID, op);

   // Perform the test.
   return VerifyVersionInfo(&osvi, VER_MAJORVERSION | VER_MINORVERSION | VER_PLATFORMID, 
                            dwlConditionMask) ? 1 : 0;
#else
  return -1;
#endif
}
