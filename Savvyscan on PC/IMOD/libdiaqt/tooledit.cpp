/*
 *  tooledit.cpp  - A subclass of QLineEdit that sends a signal when it loses 
 *                  focus.  It can also be set to a fixed width in columns.
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2006 by Boulder Laboratory for 3-Dimensional Electron
 *  Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */

#include "tooledit.h"
//Added by qt3to4:
#include <QFocusEvent>

/*!
 * A line edit widget that sends a signal, focusLost(), when it loses input 
 * focus, and that can be set to fixed column width by providing the number 
 * of characters in [columns].  [columns] defaults to 0.  [name] is ignored.
 * Use the signal editingFinished() instead of focusLost() and returnPressed().
 */
ToolEdit::ToolEdit( QWidget * parent, int columns, const char * name)
  : QLineEdit(parent)
{
  mColumns = 0;
  setColumnWidth(columns);
  connect(this, SIGNAL(editingFinished()), this, SLOT(doneEditing()));
}

ToolEdit::~ToolEdit()
{
}

void ToolEdit::changeEvent(QEvent *e)
{
  QLineEdit::changeEvent(e);
  if (e->type() == QEvent::FontChange)
    setColumnWidth();
}

// Pass on the editingFinished signal for compatibility.  It was blocked
// by catching the focusOut event
void ToolEdit::doneEditing()
{
  emit focusLost();
}

/*!
 * Set the edit box to fit the number of characters in [columns], or resize
 * it to fit a previously specified number if [columns] is zero
 */
void ToolEdit::setColumnWidth(int columns)
{
  int i, width;
  QString str;
  if (!columns)
    columns = mColumns;
  mColumns = columns;
  if (columns) {
    for (i= 0; i < columns; i++)
      str += "8";

    // Need to add 1.5 columns for right-justified edit boxes at least
    width = ((2 * columns + 3) * fontMetrics().width(str) ) / (2 * columns);
    setFixedWidth(width);
  }
}
