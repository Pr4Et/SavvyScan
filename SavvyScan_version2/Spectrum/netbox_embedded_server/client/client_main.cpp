#ifdef WIN32
#   include <winsock2.h> // before everything else!
#   include <Ws2tcpip.h>



#   undef min
#   undef max

typedef int ssize_t;
#else
#   include <sys/socket.h>
#   include <netinet/in.h>
#   include <arpa/inet.h>
#   include <unistd.h>
#   include <errno.h>
#   include <sys/types.h>
#   include <netdb.h>

typedef int SOCKET;
#define SOCKET_ERROR -1
#endif

#include <cstring>
#include <ctime>

#include <iostream>

#include <cstdlib>

// ----- include standard driver header from library -----
#include "../../c_header/dlltyp.h"

#include "../../common/ostools/spcm_oswrap.h"
#include "../../common/ostools/spcm_ostools.h"
#include "../../common/ostools/spcm_network_winLin.h"

#define SERVER_PORT     22927
#define START 0
#define STOP  1


int main ()
    {
    // ----- init network -----
    dwInitNetwork ();

    // ----- create network socket -----
    SOCKET sockfd = lSocket (AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        {
        std::cerr << "Could not create socket" << std::endl;
        return EXIT_FAILURE;
        }

    // ----- hostname or IP address of server -----
    const char* szHostname = "192.168.1.10";
    //const char* szHostname = "SERVER";

    // ----- convert hostname to IP address -----
    struct sockaddr_in serv_addr;
    memset (&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;

    int lErr = lInet_pton (szHostname, &serv_addr);
    if (lErr == SOCKET_ERROR)
        {
        std::cerr << "Could not convert IP address" << std::endl;
        lClose (sockfd);
        return EXIT_FAILURE;
        }
    serv_addr.sin_port = htons (SERVER_PORT);

    // ----- connect to server -----
    if (lConnect (sockfd, (struct sockaddr *)&serv_addr, sizeof (serv_addr)) == -1)
        {
        std::cerr << "Could not connect to server" << std::endl;
        lClose (sockfd);
        return EXIT_FAILURE;
        }

    // ----- send start command to server -----
    uint32 dwMagic = START;
    dwWrite (sockfd, &dwMagic, sizeof (uint32));

    char acBuffer[1024];
    while (true)
        {
        // ----- read data from sender -----
        unsigned dwReadBytes = dwRead (sockfd, acBuffer, sizeof (acBuffer) - 1);
        if (dwReadBytes == 0)
            {
            // connection closed
            std::cerr << "Connection closed" << std::endl;
            break;
            }
        else if (dwReadBytes == SOCKET_ERROR)
            {
            // error occured
            std::cerr << "NETWORK ERROR" << std::endl;
            break;
            }
        else
            {
            acBuffer[dwReadBytes] = '\0';
            std::cout << acBuffer << std::endl;
            }

        // check for escape = abort
        if (bKbhit ())
            {
            if (cGetch () == 27)
                {
                uint32 dwMagic = STOP;
                dwWrite (sockfd, &dwMagic, sizeof (uint32));
                break;
                }
            }
        }

    lClose (sockfd);

    return EXIT_SUCCESS;
    }
