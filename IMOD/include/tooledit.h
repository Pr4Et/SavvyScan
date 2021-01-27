/*   tooledit.h  -  declarations for tooledit.cpp
 *
 *   Copyright (C) 1995-2002 by Boulder Laboratory for 3-Dimensional Electron
 *   Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *   Colorado.  See implementation file for full copyright notice.
 *  $Id$
 */                                                                           

#ifndef TOOLEDIT_H
#define TOOLEDIT_H
#include <qlineedit.h>
//Added by qt3to4:
#include <QFocusEvent>
#include "dllexport.h"

class DLL_IM_EX ToolEdit : public QLineEdit
{
  Q_OBJECT

 public:
  ToolEdit( QWidget * parent, int columns = 0, const char * name = 0 );
  ~ToolEdit();
  void setColumnWidth(int columns = 0);

 signals:
  void focusLost();

 public slots:
  void doneEditing();

 protected:
  virtual void changeEvent(QEvent *e);

 private:
  int mColumns;    // Number of columns it is sized for, of 0 if not

};
#endif
