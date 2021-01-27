#ifndef SPCM_NETBOX_SENDMAIL_H
#define SPCM_NETBOX_SENDMAIL_H

// These files use libcurl for mail transfer

// mail content
void vNetboxMailSetServerAndLogin (const char* szSMTPServer, const char* szUsername, const char* szPassword);
void vNetboxMailSetSubject (const char* szSubject);
void vNetboxMailSetBody   (const char* szBody);
void vNetboxMailAddAttachment (const char* szFile);

void vNetboxMailSendMail (const char* szFrom, const char* szTo);

#endif // SPCM_NETBOX_SENDMAIL_H

