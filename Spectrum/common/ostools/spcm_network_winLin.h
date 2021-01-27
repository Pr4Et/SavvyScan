/******************************************************************************/
/* File name  : socketAPI.c                                                      */
/* Description: API for windows and linux                                      */
/*              This file contains the most needed function for using network */
/******************************************************************************/
#ifndef NETWORKWINLIN_H
#define NETWORKWINLIN_H

#if defined(WIN32) || defined (WIN64)
#   include <Winsock2.h>
#   include <ws2tcpip.h>
#   include <windows.h>
#   include <iphlpapi.h>

#   pragma comment(lib, "IPHLPAPI.lib")
#   define MALLOC(x) HeapAlloc(GetProcessHeap(), 0, (x))
#   define FREE(x) HeapFree(GetProcessHeap(), 0, (x))

typedef int ssize_t;
#else
#   include <netinet/in.h>
#   include <netinet/tcp.h>
#   include <arpa/inet.h>
#   include <sys/types.h>
#   include <sys/socket.h>
#   include <unistd.h>
#   include <netdb.h>
#   include <string.h>
#   include <sys/ioctl.h>
#   include <net/if.h>
#   include <ctype.h>

typedef int SOCKET;
#   define SOCKET_ERROR     -1
#   define INVALID_SOCKET   -1

#endif

#define MAX_CONN 16 // TODO: 16 guter Wert?

#define SIZEOF_IPADDR           16  //xxx.xxx.xxx.xxx
#define SIZEOF_IPADDR6          46  //xxx.xxx.xxx.xxx
#define SIZEOF_MACADDR          7  //xx xx xx xx

// ----- driver dll includes -----
#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"

#include <cstdio>
#include <cstdlib>

// our own namespace
namespace SPCM_NAMESPACE {

uint32 dwInitNetwork    ();
void   vShutdownNetwork ();

// ----- functions for connecting -----
SOCKET lSocket  (int32 lDomain, int32 lType, int32 lProtocol);
int32 lBind    (SOCKET lSockfd, const struct sockaddr* pstAddr, socklen_t plAddrlen);
int32 lListen  (SOCKET lSockfd, int32 lBacklog);
int32 lConnect (SOCKET lSockfd, const struct sockaddr* pstAddr, socklen_t plAddrlen);
SOCKET lAccept  (SOCKET lSockfd, struct sockaddr* pstAddr, socklen_t* plAddrlen);
int32 lClose   (SOCKET lSockfd);

//receiving data
ssize_t dwRead (SOCKET lSockfd, void* pvBuf, size_t dwLen);
ssize_t dwRecv (SOCKET lSockfd, void* pvBuf, size_t dwLen, int32 lFlags);
ssize_t dwRecvfrom (SOCKET lSockfd, void* pvBuf, size_t dwLen, int32 lFlags, struct sockaddr* pstSrcAddr, socklen_t* plAddrlen);

//sending data
ssize_t dwWrite (SOCKET lSockfd, const void* pvBuf, size_t dwLen);
ssize_t dwSend (SOCKET lSockfd, const void* pvBuf, size_t dwLen, int32 lFlags);
ssize_t dwSendto (SOCKET lSockfd, const void* pvBuf, size_t dwLen, int32 lFlags, const struct sockaddr* pstDestAddr, socklen_t lAddrlen);

int lSelect (SOCKET lSockfd, fd_set* pstReadFds, fd_set* pstWriteFds, fd_set* pstExceptFds, struct timeval* pstTimeout);

//info functions
int32 lSetsockopt (SOCKET lSockfd, int32 lLevel, int32 lOptname, const void* pvOptval, socklen_t lOptlen);
int32 lGetaddrinfo (const char* pcNode, const char* pcService, const struct addrinfo* pstHints, struct addrinfo** ppstRes);
int32 lGetpeername (SOCKET lSockfd, struct sockaddr* pstAddr, socklen_t* plAddrlen );
int32 lGetnameinfo (const struct sockaddr* pstSa, socklen_t lSalen, char* pcHost, size_t lHostlen, char* pcServ, size_t lServlen, int32 lFlags);
int32 lSetSocketBlocking (SOCKET lSockfd, bool bBlock);
const char* pcInet_ntop (int32 lAf, const void* pvSrc, char* pcDst, socklen_t lSize);
int lInet_pton (const char* szHostname, struct sockaddr_in* pstSockAddr);
int lInet_ntop (const struct sockaddr_in* pstSockAddr, char* szIP);

uint32 dwGetNumNICs ();
uint32 szFindOwnIpAdr (int32 lNICIdx, SOCKET lSockfd, char* szOwnIpAddress);
void vGetIPAndSubnetMask (int32 lNICIdx, SOCKET lSockfd, char* szIPAddress, char* szSubNetMask);

void vGetAdapterName (uint32 dwNICIdx, char* szName, uint32 dwBufLen);
int32 lGetMacAddress (char* szMacAddress);
int32 lGetMacAddress (uint64* pqwMacAddress);

int32 lGetNetworkErrorCode ();
char* szNetworkErrorMessage (int32 lErrorCode);

char* szIPFromName (const char* szName);

// used by netbox state monitor
const char* szIPFromSockAddr (const struct sockaddr_in* pstIP);
uint32 dwGetNetboxState (const char* szIP, uint32 dwTimeout_ms, char* szNetboxType, uint32* pdwState);

} // /SPCM_NAMESPACE

#endif //NETWORKWINLIN_H
