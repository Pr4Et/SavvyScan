/*
**************************************************************************

spcm_netbox_sendmail.cpp

**************************************************************************

This code uses libcurl for mail transfer
  
**************************************************************************
*/

#include <iostream>
#include <curl/curl.h>
#include <cstdlib>
#include <cstring>
 
#define MAX_ATTACHMENTS 16
 
static char* g_szSMTP = NULL;
static char* g_szSMTPUser = NULL;
static char* g_szSMTPPassword = NULL;
static char* g_szSubject = NULL;
static char* g_szBody = NULL;
static char* g_aszFilename[MAX_ATTACHMENTS] = {NULL};
static unsigned g_dwNextAttachmentIdx = 0;

static const int CHARS= 76;     //Sending 54 chararcters at a time with \r , \n and \0 it becomes 57
static const int ADD_SIZE= 15 + MAX_ATTACHMENTS*4;  // 15 lines for To, From, ..., plus 4 per attachment
static const int SEND_BUF_SIZE= 54;
static char (*fileBuf)[CHARS] = NULL;
static const char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
 
using namespace std;
 
bool LARGEFILE = false; /*For Percent*/
int status = 0;   /*For Percent*/
int percent2 = 0; /*For Percent*/
int percent3 = 0; /*For Percent*/
int percent4 = 0; /*For Percent*/
int percent5 = 0; /*For Percent*/
void LargeFilePercent(int rowcount)
    {
    //This is for displaying the percentage
    //while encoding larger files
    int percent = rowcount/100;

    if(LARGEFILE == true)
        {
        status++;
        percent2++;
        if (percent2 == 18)
            {
            percent3++;
            percent2 = 0;
            }
        if (percent3 == percent)
            {
            percent4++;
            percent3 = 0;
            }
        if (percent4 == 10)
            {
            system("cls");
            cout << "Larger Files take longer to encode, Please be patient." << endl
                 << endl << "Encoding, please be patient..." << endl;
            cout << percent5 << "%";
            percent5 += 10;
            percent4 = 0;
            }
        if (status == 10000)
            {
            if(percent5 == 0)
                {
                cout << " 0%"; percent5 = 10;
                }
            cout << ".";
            status = 0;
            }
        }
    }
 
void encodeblock(unsigned char in[3], unsigned char out[4], int len)
    {
    out[0] = cb64[ in[0] >> 2 ];
    out[1] = cb64[ ((in[0] & 0x03) << 4) | ((in[1] & 0xf0) >> 4) ];
    out[2] = (unsigned char) (len > 1 ? cb64[ ((in[1] & 0x0f) << 2) | ((in[2] & 0xc0) >> 6) ] : '=');
    out[3] = (unsigned char) (len > 2 ? cb64[ in[2] & 0x3f ] : '=');
    }
 
void encode(FILE *infile, unsigned char *output_buf, int rowcount/*For Percent*/)
    {
    unsigned char in[3], out[4];
        int i, len;
        *output_buf = 0;
 
    while(!feof(infile))
        {
        len = 0;
        for(i = 0; i < 3; i++)
            {
            in[i] = (unsigned char) getc(infile);
            if(!feof(infile) )
                {
                len++;
                }
            else
                {
                in[i] = 0;
                }
            }
        if(len)
            {
            encodeblock(in, out, len);
            strncat((char*)output_buf, (char*)out, 4);
            }
            LargeFilePercent(rowcount); //Display encoded file percent /*For Percent*/
        }
    }
 
 
struct fileBuf_upload_status
    {
    int lines_read;
    };
 
FILE* g_hFileTmp = NULL;
size_t read_file (char** aszFilenames, const char* szFrom, const char* szTo)
    {
    size_t len         = 0;
    size_t buffer_size = 0;
    char key = ' ';

    unsigned dwFileIdx = 0;
    struct
        {
        std::string sEncoded;
        int lNumRows;
        } astEncodedFiles[MAX_ATTACHMENTS];
    while (aszFilenames[dwFileIdx] != NULL)
        {
        //Open the file and make sure it exsits
        FILE* hFile = fopen (aszFilenames[dwFileIdx], "rb");
        if(!hFile)
            {
            cout << "File not found!!!" << endl;
            exit (EXIT_FAILURE);
            }

        //Get filesize
        fseek (hFile, 0, SEEK_END);
        size_t fileSize = ftell (hFile);
        fseek (hFile, 0, SEEK_SET);


        //Check Filesize
        if (fileSize > 256000/*bytes*/)
            {
            cout << "Larger Files take longer to encode, Please be patient." << endl;
            LARGEFILE = true; /*For Percent*/
            }
        cout << endl << "Encoding " << aszFilenames[dwFileIdx] << " please be patient..." << endl;

        //Calculate the number of rows in Base64 encoded string
        //also calculate the size of the new char to be created
        //for the base64 encode string
        astEncodedFiles[dwFileIdx].lNumRows = fileSize / SEND_BUF_SIZE + 1;
        int charsize = (astEncodedFiles[dwFileIdx].lNumRows * 72) + (astEncodedFiles[dwFileIdx].lNumRows * 2);

        //Base64 encode image and create encoded_buf string
        unsigned char* b64encode = new unsigned char[charsize];
        *b64encode = 0;
        encode (hFile, b64encode, astEncodedFiles[dwFileIdx].lNumRows/*For Percent*/);
        fclose (hFile);
        astEncodedFiles[dwFileIdx].sEncoded = (char*)b64encode;
        delete [] b64encode;

        dwFileIdx++;
        }

    if (LARGEFILE == true)
        cout << endl << endl; /*For Percent*/

    //Create structure of email to be sent
    unsigned dwNumLines = ADD_SIZE; //ADD_SIZE for TO,FROM,SUBJECT,CONTENT-TYPE,CONTENT-TRANSFER-ENCODING,CONETNT-DISPOSITION and \r\n and content
    for (unsigned i = 0; i < dwFileIdx; ++i)
        dwNumLines += astEncodedFiles[i].lNumRows;
    fileBuf = new char[dwNumLines][CHARS];
    for (unsigned i = 0; i < dwNumLines; ++i)
        fileBuf[i][0] = '\0';

#define BOUNDARY "------------060709040305030006090208"
    buffer_size += snprintf (fileBuf[len++], CHARS, "From: %s \r\n", szFrom);
    buffer_size += snprintf (fileBuf[len++], CHARS, "To: %s \r\n", szTo);
    buffer_size += snprintf (fileBuf[len++], CHARS, "Subject: %s\r\n", g_szSubject);
    buffer_size += snprintf (fileBuf[len++], CHARS, "MIME-Version: 1.0\r\n");
    if (dwFileIdx != 0)
        {
        buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Type: multipart/mixed;\r\n");
        buffer_size += snprintf (fileBuf[len++], CHARS, " boundary=\""BOUNDARY"\"\r\n");
        buffer_size += snprintf (fileBuf[len++], CHARS, "\r\n"); // empty line to divide header and data
        buffer_size += snprintf (fileBuf[len++], CHARS, "This is a multi-part message in MIME format.\r\n");
        buffer_size += snprintf (fileBuf[len++], CHARS, "--"BOUNDARY"\r\n");
        }
    buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Type: text/plain; charset=utf-8\r\n");
    buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Transfer-Encoding: 7bit\r\n\r\n");
    buffer_size += snprintf (fileBuf[len++], CHARS, "%s\r\n", g_szBody);
    buffer_size += snprintf (fileBuf[len++], CHARS, "\r\n");
    buffer_size += snprintf (fileBuf[len++], CHARS, "\r\n");
    for (unsigned i = 0; i < dwFileIdx; ++i)
        {
        buffer_size += snprintf (fileBuf[len++], CHARS, "--"BOUNDARY"\r\n");
        buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Type: application/x-msdownload; name=\"%s\"\r\n", strrchr (aszFilenames[i], '/') + 1);
        buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Transfer-Encoding: base64\r\n");
        buffer_size += snprintf (fileBuf[len++], CHARS, "Content-Disposition: attachment; filename=\"%s\"\r\n\r\n", strrchr (aszFilenames[i], '/') + 1);

        //This part attaches the Base64 encoded string and
        //sets the Base64 linesize to 72 characters + \r\n
        size_t dwPos = 0;
        std::string sub_encoded_buf;
        for (int lRow = 0; lRow <= astEncodedFiles[i].lNumRows - 1; lRow++)
            {
            sub_encoded_buf = astEncodedFiles[i].sEncoded.substr (dwPos * 72, 72);  // Reads 72 characters at a time
            sub_encoded_buf += "\r\n";                        // and appends \r\n at the end
            strcpy (fileBuf[len++], sub_encoded_buf.c_str()); // copy the 72 characters & \r\n to email
            buffer_size += sub_encoded_buf.size();            // now increase the buffer_size  
            dwPos++;                                          // finally increase pos by 1
            }
        }

// this seems to be unncessary.
//    if (dwFileIdx != 0)
//        buffer_size += snprintf (fileBuf[len++], CHARS, BOUNDARY "--");


    return buffer_size;
    }
 
 /*
  The fileBuf_source() is a function which CURL calls when it need to obtain data that will be uploaded to the server.
  Imagine that fileBuf_source() is something similar to fread(). When called it performs any voodoo-mumbo-jumbo that is needed,
  but in the end uploadable data must be stored in *ptr buffer, which is curl's internal buffer. For your in-memory buffer
  memcpy() will do just fine as body of fileBuf_source(), so you don't need real fread() at all.
 
  size * nmemb tells you how big buffer curl has reserved for a single chunk of data. The last void* is a pointer which was
  set by CURLOPT_READDATA option - it's a do-what-ever-you-need-with-it kind of pointer, so it can point to a structure
  containing data which you're uploading and a some additional info e.g. current progress.
  */
static size_t fileBuf_source(void *ptr, size_t size, size_t nmemb, void *userp)
    {
    struct fileBuf_upload_status *upload_ctx = (struct fileBuf_upload_status *)userp;
    const char *fdata;

    if((size == 0) || (nmemb == 0) || ((size*nmemb) < 1))
        {
        return 0;
        }

    fdata = fileBuf[upload_ctx->lines_read];

    if (strcmp (fdata, ""))
        {
        //printf ("Line %d: %s\n", upload_ctx->lines_read, fdata);
        size_t len = strlen(fdata);
        memcpy(ptr, fdata, len);
        upload_ctx->lines_read++;
        return len;
        }

    return 0;
    }

// ----- defines SMTP server and login credentials -----
void vNetboxMailSetServerAndLogin (const char* szSMTPServer, const char* szUsername, const char* szPassword)
    {
    if (g_szSMTP != NULL)
        delete [] g_szSMTP;
    g_szSMTP = new char[strlen (szSMTPServer) + 1];
    strcpy (g_szSMTP, szSMTPServer);

    if (g_szSMTPUser != NULL)
        delete [] g_szSMTPUser;
    g_szSMTPUser = new char[strlen (szUsername) + 1];
    strcpy (g_szSMTPUser, szUsername);

    if (g_szSMTPPassword != NULL)
        delete [] g_szSMTPPassword;
    g_szSMTPPassword = new char[strlen (szPassword) + 1];
    strcpy (g_szSMTPPassword, szPassword);
    }

// ----- defines the subject of the mail -----
void vNetboxMailSetSubject (const char* szSubject)
    {
    if (g_szSubject != NULL)
        delete [] g_szSubject;
    g_szSubject = new char[strlen (szSubject) + 1];
    strcpy (g_szSubject, szSubject);
    }

// ----- defines the body (content) of the mail -----
void vNetboxMailSetBody (const char* szBody)
    {
    if (g_szBody != NULL)
        delete [] g_szBody;
    g_szBody = new char[strlen (szBody) + 1];
    strcpy (g_szBody, szBody);
    }

// ----- adds a file that should be attached to the mail -----
// ----- call this function multiple times for multiple attachments -----
void vNetboxMailAddAttachment (const char* szFilename)
    {
    if (g_aszFilename[g_dwNextAttachmentIdx] != NULL)
        delete [] g_aszFilename[g_dwNextAttachmentIdx];
    g_aszFilename[g_dwNextAttachmentIdx] = new char[strlen (szFilename) + 1];
    strcpy (g_aszFilename[g_dwNextAttachmentIdx], szFilename);
    g_dwNextAttachmentIdx++;
    }

// ----- defines sender and recipient of the mail and starts mail transfer -----
void vNetboxMailSendMail (const char* szFrom, const char* szTo)
    {
    CURL *curl;
    CURLcode res = CURLE_OK;
    struct curl_slist *recipients = NULL;
    struct fileBuf_upload_status file_upload_ctx;
    size_t file_size(0);

    file_upload_ctx.lines_read = 0;

    file_size = read_file (g_aszFilename, szFrom, szTo);

    curl = curl_easy_init();
    if(curl)
        {
        curl_easy_setopt(curl, CURLOPT_URL, g_szSMTP);              // SMTP server address
        curl_easy_setopt(curl, CURLOPT_USERNAME, g_szSMTPUser);     // username on the SMTP server
        curl_easy_setopt(curl, CURLOPT_PASSWORD, g_szSMTPPassword); // password for the SMTP server
/*
        curl_easy_setopt(curl, CURLOPT_USE_SSL, (long)CURLUSESSL_ALL);
*/
        //curl_easy_setopt(curl, CURLOPT_CAINFO, "google.pem");
        curl_easy_setopt(curl, CURLOPT_MAIL_FROM, szFrom);
        recipients = curl_slist_append(recipients, szTo);
        curl_easy_setopt(curl, CURLOPT_MAIL_RCPT, recipients);
        curl_easy_setopt(curl, CURLOPT_INFILESIZE, file_size);
        curl_easy_setopt(curl, CURLOPT_READFUNCTION, fileBuf_source);
        curl_easy_setopt(curl, CURLOPT_READDATA, &file_upload_ctx);
        curl_easy_setopt(curl, CURLOPT_UPLOAD, 1L);
        curl_easy_setopt(curl, CURLOPT_VERBOSE, 0); //Dont display Curl Connection data Change 1L to 0
 
        res = curl_easy_perform(curl);
 
        if(res != CURLE_OK)
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
        curl_slist_free_all(recipients);
        curl_easy_cleanup(curl);

        g_dwNextAttachmentIdx = 0;
        memset (g_aszFilename, 0, sizeof (char*) * MAX_ATTACHMENTS);
        }
    delete[] fileBuf;
    }

