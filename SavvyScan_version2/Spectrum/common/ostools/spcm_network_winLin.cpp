/******************************************************************************/
/* File name  : spcm_network_winLin.cpp                                       */
/* Description: API for windows and linux                                     */
/*              This file contains the most needed function for using network */
/******************************************************************************/
#include "spcm_network_winLin.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#ifdef _LINUX
#   include <errno.h>
#   include <fcntl.h>
#   include <ifaddrs.h>
#endif

#ifdef _QNX
#   include <errno.h>
#   include <fcntl.h>
#endif

#ifdef WIN32
#   undef max
#   undef min
#endif

#define SPCM_MAX_BYTES_TO_COPY      (1024*1024)/*64*1024*/

// we're using this namespace in our program
namespace SPCM_NAMESPACE {

/*
    Initialize network.
*/
uint32 dwInitNetwork ()
    {
#ifdef WIN32
    WSADATA wsa;
    if (WSAStartup (MAKEWORD (2,0), &wsa) != 0)
        {
        return WSANOTINITIALISED;
        }
#else
    // nothing
#endif
    return ERR_OK;
    }

void vShutdownNetwork ()
    {
#ifdef WIN32
    WSACleanup ();
#else
    // nothing
#endif
    }


/*
    Create a socket. 
    lDomain = An address format specification. The only format currently supported is PF_INET, which is the ARPA Internet address format.
    lType = A type specification for the new socket.
    lProtocol = A particular protocol to be used with the socket, or 0 if the caller does not wish to specify a protocol.
*/
SOCKET lSocket (int32 lDomain, int32 lType, int32 lProtocol)
    {
    return socket (lDomain, lType, lProtocol);
    }


/*
    Close a socket.
    lSockfd = A (file) descriptor identifying a socket.
*/
int32 lClose (SOCKET lSockfd)
    {
#ifdef WIN32
        return closesocket (lSockfd);
#else
        return close (lSockfd);
#endif
    }

/*
    Establish a connection to a peer.
    lSockfd = A socket (file) descriptor identifying an unconnected socket.
    pstAddr = The address to assign to the socket. The address have to be casted to "struct sockaddr *".
    lAddrlen = The length of the name.
*/
int32 lConnect (SOCKET lSockfd, const struct sockaddr* pstAddr, socklen_t lAddrlen)
    {
    return connect (lSockfd, pstAddr, lAddrlen);
    }

/* 
    Read data from a socket.
    lSockfd = A socket (file) descriptor identifying a connected socket. 
    pvBuf = A buffer for the incoming data. 
    lLen = The length of buf. 
*/
ssize_t dwRead (SOCKET lSockfd, void* pvBuf, size_t dwLen)
    {
#ifdef WIN32
    return recv (lSockfd, (char*)pvBuf, (int)dwLen, 0);
#else
    return read (lSockfd, pvBuf, dwLen);
#endif
    }

/*
    Receive data from a socket.
    lSockfd = A socket (file) descriptor identifying a connected socket. 
    pvBuf = A buffer for the incoming data. 
    lLen = The length of buf. 
    lFlags = Specifies the way in which the call is made.
*/
ssize_t dwRecv (SOCKET lSockfd, void* pvBuf, size_t dwLen, int32 lFlags)
    {
#ifdef WIN32
    return recv (lSockfd, (char*)pvBuf, (int)dwLen, lFlags);
#else
    return recv (lSockfd, (char*)pvBuf, dwLen, lFlags);
#endif
    }

/*
    Receive a datagram and store the source address.
    lSockfd = A socket (file) descriptor identifying a bound socket. 
    pvBuf = A buffer for the incoming data. 
    lLen = The length of buf. 
    lFlags = Specifies the way in which the call is made. 
    pstSrc_addr = An optional pointer to a buffer which will hold the source address upon return. 
    plAddrlen = An optional pointer to the size of the from buffer.
*/
ssize_t dwRecvfrom (SOCKET lSockfd, void* pvBuf, size_t dwLen, int32 lFlags, struct sockaddr *pstSrc_addr, socklen_t* plAddrlen)
    {
#ifdef WIN32
    return recvfrom (lSockfd, (char*)pvBuf, (int)dwLen, lFlags, pstSrc_addr, plAddrlen);
#else
    return recvfrom (lSockfd, (char*)pvBuf, dwLen, lFlags, pstSrc_addr, plAddrlen);
#endif
    }

/*
    Send data on a connected socket.
    lSockfd = A socket (file) descriptor identifying a connected socket. 
    pvBuf = A buffer containing the data to be transmitted. 
    lLen = The length of the data in buf. 
*/
ssize_t dwWrite (SOCKET lSockfd, const void* pvBuf, size_t dwLen)
    { 
#ifdef WIN32
    return send (lSockfd, (char*)pvBuf, (int)dwLen, 0);
#else
        return write (lSockfd, pvBuf, dwLen);
#endif
    }

/*
    Send data on a connected socket.
    lSockfd = A socket (file) descriptor identifying a connected socket. 
    pvBuf = A buffer containing the data to be transmitted. 
    lLen = The length of the data in buf. 
    lFlags = Specifies the way in which the call is made.
*/
ssize_t dwSend (SOCKET lSockfd, const void* pvBuf, size_t dwLen, int32 lFlags)
    {
#ifdef WIN32
    return send (lSockfd, (char*)pvBuf, (int)dwLen, lFlags);
#else
    return send (lSockfd, (char*)pvBuf, dwLen, lFlags);
#endif
    }


/*
    Send data to a specific destination.
    lSockfd = A socket (file) descriptor identifying a socket. 
    pvBuf = A buffer containing the data to be transmitted. 
    lLen = The length of the data in buf. 
    lFlags = Specifies the way in which the call is made. 
    pstDest_addr = A optional pointer to the address of the target socket. 
    lAddrlen = The size of the address in to.
*/
ssize_t dwSendto (SOCKET lSockfd, const void* pvBuf, size_t dwLen, int32 lFlags, const struct sockaddr* pstDest_addr, socklen_t lAddrlen)
    {
#ifdef WIN32
    return sendto (lSockfd, (char*)pvBuf, (int)dwLen, lFlags, pstDest_addr, lAddrlen);
#else
    return sendto (lSockfd, (char*)pvBuf, dwLen, lFlags, pstDest_addr, lAddrlen);
#endif
    }

int lSelect (SOCKET lSockfd, fd_set* pstReadFds, fd_set* pstWriteFds, fd_set* pstExceptFds, struct timeval* pstTimeout)
    {
#ifdef WIN32
    return select (0/*unused*/, pstReadFds, pstWriteFds, pstExceptFds, pstTimeout);
#else
    return select (lSockfd, pstReadFds, pstWriteFds, pstExceptFds, pstTimeout);
#endif
    }

/*
    Accept a connection on a socket
    lSockfd = A socket (file) descriptor identifying a socket which is listening for connections after a listen(). 
    pstAddr = An optional pointer to a buffer which receives the address of the connecting entity, as known to the communications layer. 
        The exact format of the addr argument is determined by the address family established when the socket was created. 
    plAddrlen = A optional pointer to an integer which contains the length of the address addr.
*/

SOCKET lAccept (SOCKET lSockfd, struct sockaddr* pstAddr, socklen_t* plAddrlen)
    {
    return accept (lSockfd, pstAddr, plAddrlen);
    }


/*
    Associate a local address with a socket. 
    lSockfd = A socket (file) descriptor identifying an unbound socket.
    pstAddr = The address to assign to the socket. The address have to be casted to "struct sockaddr *".
    lAddrlen = The length of the name.
*/
int32 lBind (SOCKET lSockfd, const struct sockaddr* pstAddr, socklen_t lAddrlen)
    {
    return bind (lSockfd, pstAddr, lAddrlen);
    }


/*
    Establish a socket to listen for incoming connection. 
    lSockfd = A socket (file) descriptor identifying a bound, unconnected socket. 
    lBacklog = The maximum length to which the queue of pending connections may grow.
*/
int32 lListen (SOCKET lSockfd, int32 lBacklog)
    {
    return listen (lSockfd, lBacklog);
    }   

/*
    Set options on sockets.
    lSockfd = A socket (file) descriptor identifying a bound, unconnected socket. 
    lLevel = The level at which the option is defined (for example, SOL_SOCKET).
    lOptname = The socket option for which the value is to be set (for example, SO_BROADCAST). 
        The optname parameter must be a socket option defined within the specified level, or behavior is undefined.
    pvOptval = A pointer to the buffer in which the value for the requested option is specified.
    lOptlen = The size, in bytes, of the buffer pointed to by the optval parameter.
*/
int32 lSetsockopt (SOCKET lSockfd, int32 lLevel, int32 lOptname, const void* pvOptval, socklen_t lOptlen)
    {
#ifdef WIN32
        return setsockopt (lSockfd, lLevel, lOptname, (const char*) pvOptval, lOptlen);
#else
        return setsockopt (lSockfd, lLevel, lOptname, pvOptval, lOptlen);
#endif
    }

/*
    provides protocol-independent translation from an ANSI host name to an address.
    cpNode = A pointer to a NULL-terminated ANSI string that contains a host (node) name or a numeric host address string. 
        For the Internet protocol, the numeric host address string is a dotted-decimal IPv4 address or an IPv6 hex address.
    cpService = A pointer to a NULL-terminated ANSI string that contains either a service name or port number represented as a string.
        A service name is a string alias for a port number.
    pstHints = A pointer to an addrinfo structure that provides hints about the type of socket the caller supports. 
    ppstRes = A pointer to a linked list of one or more addrinfo structures that contains repstonse information about the host.
*/

int32 lGetaddrinfo (const char* pcNode, const char* pcService, const struct addrinfo* pstHints, struct addrinfo** ppstRes)
    {
    return getaddrinfo (pcNode, pcService, pstHints, ppstRes);
    }

/*
    retrieves the address of the peer to which a socket is connected.
    dwSockfd = A descriptor identifying a connected socket.
    sAddr = The SOCKADDR structure that receives the address of the peer.
    lAddrlen = A pointer to the size, in bytes, of the name parameter.
*/
int32 lGetpeername (SOCKET lSockfd, struct sockaddr* pstAddr, socklen_t* plAddrlen)
    {
    return getpeername (lSockfd,  pstAddr, plAddrlen);
    }

/*
    address-to-name translation in protocol-independent manner
    pstSa = A pointer to a socket address structure that contains the address and port number of the socket. 
        For IPv4, the sa parameter points to a sockaddr_in structure. For IPv6, the sa parameter points to a sockaddr_in6 structure.
    lSaLen = The length, in bytes, of the structure pointed to by the sa parameter.
    pcHost = A pointer to an ANSI string used to hold the host name. 
        On success, the host name is returned as a Fully Qualified Domain Name (FQDN) by default. 
        If the host parameter is NULL, this indicates the caller does not want to receive a host name string.
    sHostlen = The length, in bytes, of the buffer pointed to by the pcHost parameter. 
        The caller must provide a buffer large enough to hold the host name, including the terminating NULL character.
    pcServ = A pointer to an ANSI string to hold the service name. 
        On success, an ANSI string that represents the service name associated with the port number is returned. 
        If the serv parameter is NULL, this indicates the caller does not want to receive a service name string.
    lServLen = The length, in bytes, of the buffer pointed to by the pcServ parameter. 
        The caller must provide a buffer large enough to hold the service name, including the terminating NULL character.
    lFlags = A value used to customize processing of the getnameinfo function. See the Remarks section.
*/
int32 lGetnameinfo (const struct sockaddr* pstSa, socklen_t lSaLen, char* pcHost, size_t dwHostlen, char* pcServ, size_t dwServLen, int32 lFlags)
    {
#ifdef WIN32
    return getnameinfo (pstSa, lSaLen, pcHost, static_cast < DWORD > (dwHostlen), pcServ, static_cast < DWORD > (dwServLen), lFlags);
#else
    return getnameinfo (pstSa, lSaLen, pcHost, dwHostlen, pcServ, dwServLen, lFlags);
#endif
    }


// ****************************************************************************
// ***** set blocking mode of socket
// ****************************************************************************
int32 lSetSocketBlocking (SOCKET lSockfd, bool bBlock)
    {
#ifdef WIN32
    uint32 dwNonblockMode = (bBlock? 0 : 1);
    return ioctlsocket (lSockfd, FIONBIO, &dwNonblockMode);
#else
    int lArg = fcntl (lSockfd, F_GETFL, NULL); 
    if (lArg == -1)
        return SOCKET_ERROR;
    if (bBlock)
        lArg &= ~O_NONBLOCK;
    else
        lArg |= O_NONBLOCK; 
    if (fcntl (lSockfd, F_SETFL, lArg) == -1)
        return SOCKET_ERROR;
    return 0;
#endif
    }

// ----- get the number of available network interfaces from the system -----
uint32 dwGetNumNICs ()
    {
#if defined (_LINUX) || defined (_QNX)
    // ----- create temp socket for request -----
    int32 lTmpSock = socket (AF_INET, SOCK_STREAM, 0);
    if (lTmpSock == SOCKET_ERROR)
        return 0;

    // ----- point ifconf's ifc_buf to our array of interface ifreqs -----
    struct ifconf stIfconf;
    struct ifreq astIfreq[20]; // random maximum of interfaces
    stIfconf.ifc_buf = (char *) astIfreq;

    // Set ifconf's ifc_len to the length of our array of interface ifreqs.
    stIfconf.ifc_len = sizeof (astIfreq);

    //  Populate ifconf.ifc_buf (ifreq) with a list of interface names and addresses.
    if (ioctl (lTmpSock, SIOCGIFCONF, &stIfconf) == -1)
        return 0;

    close (lTmpSock);

    // Divide the length of the interface list by the size of each entry.
    // This gives us the number of interfaces on the system.
    return stIfconf.ifc_len / sizeof (struct ifreq);
#else
    uint32 dwBufLen = 0;
    GetInterfaceInfo (NULL, &dwBufLen); // get necessary buffer length
    char* acBuf = new char[dwBufLen];
    IP_INTERFACE_INFO* pstIPInterfaceInfo = reinterpret_cast < IP_INTERFACE_INFO* > (acBuf);
    GetInterfaceInfo (pstIPInterfaceInfo, &dwBufLen); // get infos

    uint32 dwNumNICs = pstIPInterfaceInfo->NumAdapters;
    delete [] acBuf;
    return dwNumNICs;
#endif
    }

uint32 szFindOwnIpAdr (int32 lNICIdx, SOCKET lSockfd, char* szOwnIpAddress)
    {
#ifdef _LINUX
        // ----- request info structure from kernel -----
        struct ifconf stIfconf;
        struct ifreq astIfreq[20]; // random maximum of interfaces
        stIfconf.ifc_buf = (char *) astIfreq;
        stIfconf.ifc_len = sizeof (astIfreq);
        if (ioctl (lSockfd, SIOCGIFCONF, &stIfconf) == -1)
            return ERR_ABORT;

        const uint32 dwNumInterfaces = stIfconf.ifc_len / sizeof (struct ifreq);
        if (static_cast < uint32 > (lNICIdx) >= dwNumInterfaces)
            return ERR_ABORT;

        struct sockaddr_in *pstAddress = (struct sockaddr_in *) &astIfreq[lNICIdx].ifr_addr;
        if (pstAddress == NULL)
            return ERR_ABORT;

        // Convert the binary IP address into a readable string.
        if (inet_ntop (AF_INET, &pstAddress->sin_addr, szOwnIpAddress, SIZEOF_IPADDR) == NULL)
            return ERR_ABORT;

#else
        char szBuf[255];
        if (gethostname (szBuf, 255) != SOCKET_ERROR)
            {
#   if defined(_MSC_VER) && (_MSC_VER < 1900) // 1900 = VS2015
            struct hostent* pstHostend = gethostbyname (szBuf);
            if (pstHostend != NULL)
                {
                for (int32 lCnt = 0; lCnt < lNICIdx + 1; ++lCnt)
                    {
                    if (pstHostend->h_addr_list[lCnt] == NULL)
                        {
                        return ERR_ABORT; // no valid address found
                        }
                    }

                if (pstHostend->h_addr_list[lNICIdx])
                    {
                    struct in_addr* pstAddr = reinterpret_cast < struct in_addr* > (pstHostend->h_addr_list[lNICIdx]);

                    char* szIP = inet_ntoa (*pstAddr);
                    if (szIP == NULL)
                        return ERR_ABORT; // TODO

                    memcpy (szOwnIpAddress, szIP, SIZEOF_IPADDR);
                    }
                }
#   else
            struct addrinfo stHints;
            memset (&stHints, 0, sizeof (stHints));
            stHints.ai_family = AF_INET; // do not show IPv6 for now

            struct addrinfo* pstAddrInfo = NULL;
            getaddrinfo (szBuf, NULL, &stHints, &pstAddrInfo);
            struct addrinfo* pstCurAddrInfo = pstAddrInfo;

            // ----- find struct for NIC -----
            for (int32 lCnt = 0; lCnt < lNICIdx && pstCurAddrInfo; ++lCnt)
                {
                if (pstCurAddrInfo->ai_next == NULL)
                    return ERR_ABORT;
                pstCurAddrInfo = pstCurAddrInfo->ai_next;
                }

            // ----- convert binary IP to string -----
            switch (pstCurAddrInfo->ai_family)
                {
                case AF_INET:
                    InetNtopA (AF_INET, &((struct sockaddr_in*)pstCurAddrInfo->ai_addr)->sin_addr, szOwnIpAddress, SIZEOF_IPADDR);
                    break;
                case AF_INET6:
                    InetNtopA (AF_INET6, &((struct sockaddr_in6*)pstCurAddrInfo->ai_addr)->sin6_addr, szOwnIpAddress, SIZEOF_IPADDR6); // should not happen with stHints set above
                    break;
                }
            freeaddrinfo (pstAddrInfo);
#   endif
            }
        
#endif

    return ERR_OK;
    }

void vGetIPAndSubnetMask (int32 lNICIdx, SOCKET lSockfd, char* szIPAddress, char* szSubNetMask)
    {
    if (szIPAddress != NULL)
        memset (szIPAddress,  0, SIZEOF_IPADDR);
    if (szSubNetMask != NULL)
        memset (szSubNetMask, 0, SIZEOF_IPADDR);

    #ifdef WIN32

    // get list of network adapters
    IP_ADAPTER_INFO* pstFirstAdapterInfo = NULL;
    bool bDelete = false;
    ULONG dwLen = 0;
    DWORD dwErr = GetAdaptersInfo (pstFirstAdapterInfo, &dwLen);
    if (dwErr == ERROR_BUFFER_OVERFLOW)
        {
        pstFirstAdapterInfo = (IP_ADAPTER_INFO*)malloc (dwLen); // use malloc because of heap corruption when using new
        dwErr = GetAdaptersInfo (pstFirstAdapterInfo, &dwLen);
        bDelete = true;
        }
    if (dwErr != ERROR_SUCCESS)
        {
        if (bDelete)
            free (pstFirstAdapterInfo);
        return;
        }

    // find the requested adapter
    IP_ADAPTER_INFO* pstCurrAdapterInfo = pstFirstAdapterInfo;
    for (uint32 dwIdx = 0; dwIdx < static_cast < uint32 > (lNICIdx) && pstCurrAdapterInfo->Next; ++dwIdx)
        {
        pstCurrAdapterInfo = pstCurrAdapterInfo->Next;
        }

    if (szIPAddress != NULL)
        memcpy (szIPAddress,  &pstCurrAdapterInfo->IpAddressList.IpAddress, SIZEOF_IPADDR);
    if (szSubNetMask != NULL)
        memcpy (szSubNetMask, &pstCurrAdapterInfo->IpAddressList.IpMask,    SIZEOF_IPADDR);

    free (pstFirstAdapterInfo);

    #else
        // using getifaddrs to get subnet mask did not work for unknown reasons, so we use ioctl calls

        // ----- request info structure from kernel -----
        struct ifconf stIfconf;
        struct ifreq astIfreq[20]; // random maximum of interfaces
        stIfconf.ifc_buf = (char *) astIfreq;
        stIfconf.ifc_len = sizeof (astIfreq);
        if (ioctl (lSockfd, SIOCGIFCONF, &stIfconf) == -1)
            return;

        const uint32 dwNumInterfaces = stIfconf.ifc_len / sizeof (struct ifreq);
        if (static_cast < uint32 > (lNICIdx) >= dwNumInterfaces)
            return;

        struct sockaddr_in* pstAddress = (struct sockaddr_in *) &astIfreq[lNICIdx].ifr_addr;
        if (pstAddress == NULL)
            return;

        // Convert the binary IP address into a readable string.
        if (inet_ntop (AF_INET, &pstAddress->sin_addr, szIPAddress, SIZEOF_IPADDR) == NULL)
            return;

        // ----- Subnetmask -----
        if (ioctl (lSockfd, SIOCGIFNETMASK, astIfreq + lNICIdx) == -1)
            return;

        pstAddress = (struct sockaddr_in *) &astIfreq[lNICIdx].ifr_addr;

        // Convert the binary subnet mask into a readable string.
        if (inet_ntop (AF_INET, &pstAddress->sin_addr, szSubNetMask, SIZEOF_IPADDR) == NULL)
            return;

    #endif
    }

void vGetAdapterName (uint32 dwNICIdx, char* szName, uint32 dwBufLen)
    {
    memset (szName, '\0', dwBufLen);
#if defined (_LINUX) || defined (_QNX)
    // ----- create temp socket for request -----
    int32 lTmpSock = socket (AF_INET, SOCK_STREAM, 0);
    if (lTmpSock == SOCKET_ERROR)
        return;

    // ----- point ifconf's ifc_buf to our array of interface ifreqs -----
    struct ifconf stIfconf;
    struct ifreq astIfreq[20]; // random maximum of interfaces
    stIfconf.ifc_buf = (char *) astIfreq;

    // Set ifconf's ifc_len to the length of our array of interface ifreqs.
    stIfconf.ifc_len = sizeof (astIfreq);

    //  Populate ifconf.ifc_buf (ifreq) with a list of interface names and addresses.
    if (ioctl (lTmpSock, SIOCGIFCONF, &stIfconf) == -1)
        return;

    close (lTmpSock);

    const uint32 dwNumAdapters = stIfconf.ifc_len / sizeof (struct ifreq);
    if (dwNICIdx >= dwNumAdapters)
        return;

    strncpy (szName, stIfconf.ifc_req[dwNICIdx].ifr_name, dwBufLen - 1);
#else
    IP_ADAPTER_INFO* pstFirstAdapterInfo = NULL;
    IP_ADAPTER_INFO* pstCurrAdapterInfo  = NULL;
    bool bDelete = false;
    ULONG dwLen = 0;
    DWORD dwErr = GetAdaptersInfo (pstFirstAdapterInfo, &dwLen);
    if (dwErr == ERROR_BUFFER_OVERFLOW)
        {
        pstFirstAdapterInfo = (IP_ADAPTER_INFO*)malloc (dwLen); // use malloc because of heap corruption when using new
        dwErr = GetAdaptersInfo (pstFirstAdapterInfo, &dwLen);
        bDelete = true;
        }
    if (dwErr != ERROR_SUCCESS)
        {
        if (bDelete)
            free (pstFirstAdapterInfo);
        return;
        }

    pstCurrAdapterInfo = pstFirstAdapterInfo;
    for (uint32 dwIdx = 0; dwIdx < dwNICIdx && pstCurrAdapterInfo->Next; ++dwIdx)
        {
        pstCurrAdapterInfo = pstCurrAdapterInfo->Next;
        }

    strncpy (szName, pstCurrAdapterInfo->Description, dwBufLen - 1);

    if (bDelete)
        free (pstFirstAdapterInfo);
#endif
    }


int32 lGetMacAddress (uint64* pqwMacAddress)
    {
    *pqwMacAddress = 0;

#ifdef WIN32
    /* Declare and initialize variables */
    uint32 dwError = 0;

    // Set the flags to pass to GetAdaptersAddresses
    uint32 dwFlags = GAA_FLAG_INCLUDE_PREFIX;

    // default to unspecified address family (both)
    uint32 dwFamily = AF_UNSPEC;
    PIP_ADAPTER_ADDRESSES pstAddresses = NULL;
    uint32 dwBufLen = 0;

    dwFamily = AF_INET;

    dwBufLen = sizeof (IP_ADAPTER_ADDRESSES);
    pstAddresses = (IP_ADAPTER_ADDRESSES *) MALLOC (dwBufLen);

    // Make an initial call to GetAdaptersAddresses to get the
    // size needed into the outBufLen variable
    if (GetAdaptersAddresses (dwFamily, dwFlags, NULL, pstAddresses, &dwBufLen) == ERROR_BUFFER_OVERFLOW) 
        {
        FREE (pstAddresses);
        pstAddresses = (IP_ADAPTER_ADDRESSES *) MALLOC (dwBufLen);
        }

    if (pstAddresses == NULL) 
        {
        printf("Memory allocation failed for IP_ADAPTER_ADDRESSES struct\n");
        exit(1);
        }

    // Make a second call to GetAdapters Addresses to get the
    // actual data we want
    dwError = GetAdaptersAddresses (dwFamily, dwFlags, NULL, pstAddresses, &dwBufLen);

    // take first adapter
    memcpy (pqwMacAddress, pstAddresses->PhysicalAddress, SIZEOF_MACADDR - 1);
    FREE(pstAddresses);

#else
    FILE* pFile = fopen ("/sys/class/net/eth0/address", "r");

    char szMAC[18]; // 00:03:2D:21:AE:AE\0
    char* szTmp __attribute__((unused)) = fgets (szMAC, sizeof (szMAC), pFile);
    fclose (pFile);

    char* pcPos = szMAC;
    uint32 dwShift = 40;
    for (int lIdx = 0; lIdx < 6; ++lIdx)
        {
        *pqwMacAddress |= static_cast < uint64 > (strtoul (pcPos, NULL, 16)) << dwShift;
        dwShift -= 8;
        pcPos += 3;
        }
#endif

    return ERR_OK;
    }

int32 lGetMacAddress (char* szMacAddress)
    {
    uint64 qwMacAddress = 0;
    lGetMacAddress (&qwMacAddress);

#ifdef WIN32
    sprintf (szMacAddress, "%2lX:%2lX:%2lX:%2lX:%2lX:%2lX",
#else
    sprintf (szMacAddress, "%2X:%2X:%2X:%2X:%2X:%2X",
#endif
        static_cast < uint32 > ((qwMacAddress >> 40) & 0xFF),
        static_cast < uint32 > ((qwMacAddress >> 32) & 0xFF),
        static_cast < uint32 > ((qwMacAddress >> 24) & 0xFF),
        static_cast < uint32 > ((qwMacAddress >> 16) & 0xFF),
        static_cast < uint32 > ((qwMacAddress >>  8) & 0xFF),
        static_cast < uint32 > ((qwMacAddress >>  0) & 0xFF));

    return ERR_OK;
    }



int32 lGetNetworkErrorCode ()
    {
#ifdef WIN32
    return WSAGetLastError ();
#else
    return errno;
#endif
    }


// ----- Returns a string with the last network error message -----
char* szNetworkErrorMessage (int32 lErrorCode)
    {
#ifdef WIN32
    static char s_szErrMsgBuffer[200];

    LPSTR szErrorText = NULL;
    FormatMessageA (FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_ALLOCATE_BUFFER, NULL, lErrorCode, MAKELANGID(LANG_NEUTRAL, SUBLANG_NEUTRAL), (LPSTR)&szErrorText, 0, NULL);

    // convert to char*
    _snprintf_s (s_szErrMsgBuffer, 100, _TRUNCATE, "%s", szErrorText);

    // release memory allocated by FormatMessage()
    LocalFree (szErrorText);

    return s_szErrMsgBuffer;
#else
    return strerror (lErrorCode);
#endif
    }


// ----- converts an IP representation to string -----
const char* szIPFromSockAddr (const struct sockaddr_in* pstIP)
    {
    static char s_szIPBuffer[16];

#ifdef WIN32
    sprintf (s_szIPBuffer, "%u.%u.%u.%u", pstIP->sin_addr.S_un.S_un_b.s_b1, pstIP->sin_addr.S_un.S_un_b.s_b2, pstIP->sin_addr.S_un.S_un_b.s_b3, pstIP->sin_addr.S_un.S_un_b.s_b4);
#else
    uint32 dwIP = ntohl (pstIP->sin_addr.s_addr);
    sprintf (s_szIPBuffer, "%u.%u.%u.%u", (dwIP & 0xFF000000)>>24, (dwIP & 0xFF0000)>>16, (dwIP & 0xFF00)>>8, (dwIP & 0xFF));
#endif
    return s_szIPBuffer;
    }

char* szIPFromName (const char* szName)
    {
    static char szIPBuffer[16]; //123.123.123.123\0

    struct addrinfo* pstResult = NULL;

    int lErr = getaddrinfo (szName, NULL, NULL, &pstResult);
    if (lErr != 0)
        {
        szIPBuffer[0] = '\0';
        return szIPBuffer;
        }

    struct sockaddr_in *ipv = (struct sockaddr_in *)pstResult->ai_addr;

    // Convert the binary IP address into a readable string.
    if (lInet_ntop (ipv, szIPBuffer) == SOCKET_ERROR)
        {
        freeaddrinfo (pstResult);
        return NULL;
        }
    freeaddrinfo (pstResult);

    return szIPBuffer;
    }

int lInet_pton (const char* szHostname, struct sockaddr_in* pstSockAddr)
    {
#if defined(_MSC_VER) && (_MSC_VER < 1900) // for compilation with older VS. 1900 = VS2015
    pstSockAddr->sin_addr.s_addr = inet_addr (szHostname);
#else
    if (inet_pton (AF_INET, szIPFromName (szHostname), &(pstSockAddr->sin_addr)) <= 0)
        {
        return SOCKET_ERROR;
        }
#endif
    return 0;
    }

int lInet_ntop (const struct sockaddr_in* pstSockAddr, char* szIP)
    {
#if defined(_MSC_VER) && (_MSC_VER < 1900) // for compilation with older VS. 1900 = VS2015
    strncpy (szIP, inet_ntoa (pstSockAddr->sin_addr), SIZEOF_IPADDR);
#else
    if (inet_ntop (AF_INET, &(const_cast < struct sockaddr_in* > (pstSockAddr)->sin_addr), szIP, SIZEOF_IPADDR) == NULL)
        {
        return SOCKET_ERROR;
        }
#endif
    return 0;
    }

} // /SPCM_NAMESPACE
