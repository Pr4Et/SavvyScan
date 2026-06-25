#include <qapplication.h>
#include <qfile.h>
#include <qstringlist.h>
#include <stdlib.h>
#include "dia_qtutils.h"

static bool noInitSet = false;

/*!
 * Sets up the library path for Qt so that it loads plugins from the IMOD installation.
 * For Qt 5, this needs to be called before starting the QApplication so that platform
 * plugins can be found.
 */
void diaSetQtLibraryPath(void)
{
  QStringList strList;
  int i;
  char *plugdir;

  // Put plugin dir on the library path so image plugins can be found
  // Replace the default path so we don't get anything from our own Qt install,
  // Unless we are running an IMOD installation without an imageformats directory
  plugdir = getenv("IMOD_PLUGIN_DIR");
  if (plugdir)
    strList << plugdir;
  plugdir = getenv("IMOD_DIR");
  if (plugdir)
    strList << QString(plugdir) + QString("/lib/imodplug");

  // Add these for standalone 3dmod package
#ifdef _WIN32
  else
    strList << QString("C:/Program Files/IMOD/lib/imodplug") << 
      QString("C:/Program Files/3dmod/lib/imodplug");
#endif

  // Or this for running by clicking on it
#ifdef Q_OS_MACX
  else
    strList << QString("/Applications/IMOD/lib/imodplug");
#endif
  for (i = 0; i < strList.count(); i++)
    if (QFile::exists(strList[i] + "/imageformats"))
      break;
  if (strList.count() && i < strList.count())
    QApplication::setLibraryPaths(strList);
  diaQtNoInitOMPonFork();
}

/*!
 * Works around a hang when starting a QProcess with version 2017 or 2018 of the Intel 
 * OpenMP library for Linux by setting environment variable KMP_INIT_AT_FORK=FALSE. 
 * This is called by @diaSetQtLibraryPath but may need to be called explicitly if that 
 * function is called after forking.
 */
void diaQtNoInitOMPonFork()
{
#ifdef Q_OS_LINUX
  if (!noInitSet)
    setenv("KMP_INIT_AT_FORK", "FALSE", 1);
#endif
  noInitSet = true;
}
