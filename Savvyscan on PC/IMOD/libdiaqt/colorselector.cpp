/* 
 *  colorselector.cpp       Implementation of color selector class
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2006 by Boulder Laboratory for 3-Dimensional Electron
 *  Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */
#include <qlayout.h>
#include <qlabel.h>
#include <qcolordialog.h>
//Added by qt3to4:
#include <QTimerEvent>
#include <QKeyEvent>
#include <QCloseEvent>
#include "multislider.h"
#include "colorselector.h"
#include "dia_qtutils.h"

#define HOT_SLIDER_KEYUP 0
#define HOT_SLIDER_KEYDOWN 1
#define NO_HOT_SLIDER 2

static const char *sliderLabels[] = {"Red", "Green", "Blue"};
static const char *buttonLabels[] = {"Done", "Restore", "Qt Selector"};
static const char *buttonTips[] = 
  {"Close color selector",
   "Restore to starting color when selector was opened",
   "Select color with Qt color selector box"};

/*!
 * This class provides a color selector with a sample color panel, and three
 * sliders for adjusting red, green, and blue.  [label] is used to set a
 * label at the top of the panel, and the color is initialized with [red], 
 * [green], and [blue].  It manages the color of the
 * panel continuously during changes, and emits a signal for a new color
 * if the slider is clicked.  It will also emit signals during a drag if 
 * [hotFlag] is not 2; if the key given by [hotKey] is up when [hotFlag] is 0;
 * or if that key is down when [hotFlag] is 1.  [name] defaults to NULL, and
 * the window flags default to Qt::WDestructiveClose | Qt::WType_TopLevel.
 * ^     Signals emitted are:
 * ^ void newColor(int r, int g, int b);  -  When the color changes
 * ^ void done();   -  When the Done button is pressed
 * ^ void closing();  -  When the window is closing
 * ^ void keyPress(QKeyEvent *e);   -  When a key is pressed
 * ^ void keyRelease(QKeyEvent *e);  -  When a key is released
 * ^     In addition, there is one method:
 * ^ bool hotSliding();  -  Returns true is a slider is being dragged
 */
ColorSelector::ColorSelector(QWidget *parent, QString label, int red,
                             int green, int blue, int hotFlag, int hotKey,
                             bool rounded, const char *name, Qt::WindowFlags fl)
  : DialogFrame(parent, 3, 1, buttonLabels, buttonTips, false, rounded, "test",
                "test2", name, fl)
{
  QString str;

  mOriginalRGB[0] = mCurrentRGB[0] = red;
  mOriginalRGB[1] = mCurrentRGB[1] = green;
  mOriginalRGB[2] = mCurrentRGB[2] = blue;
  mCtrlPressed = false;
  mHotKey = hotKey;
  mHotFlag = hotFlag;
  mDragging = false;

  // Get the top label
  QLabel *topLabel = new QLabel(label, this);
  mLayout->addWidget(topLabel);
  
  // 3/12/18: Switch back to frame with background coloring; GL not needed
  mColorBox = new QFrame(this);
  mColorBox->setFrameStyle(QFrame::Plain);
  mColorBox->setFixedHeight(50);
  mLayout->addWidget(mColorBox);
  mColorBox->setAutoFillBackground(true);
 
  // Get the sliders, connect them and initialize them to current color
  mSliders = new MultiSlider(this, 3, sliderLabels);
  mLayout->addLayout(mSliders->getLayout());
  connect(mSliders, SIGNAL(sliderChanged(int, int, bool)), this, 
          SLOT(sliderChanged(int, int, bool)));

  // Connect them: have to connect to release of Qt selector because the modal
  // box keeps the button from coming back up (maybe mixed X problem only)
  connect(this, SIGNAL(actionClicked(int)), this, SLOT(buttonClicked(int)));

  imposeColor(true, false);
}

ColorSelector::~ColorSelector()
{
}

void ColorSelector::buttonClicked(int which)
{
  if (which == 0)
    donePressed();
  else if (which == 1)
    restorePressed();
  else if (which == 2)
    qtSelectorPressed();
}

void ColorSelector::donePressed()
{
  emit done();
}

// Restore: restore the original color and send signal
void ColorSelector::restorePressed()
{
  for (int i = 0; i < 3; i++)
    mCurrentRGB[i] = mOriginalRGB[i];
  imposeColor(true, true);
}

// Open color selector, and take color if it is valid
void ColorSelector::qtSelectorPressed()
{
  QColor retColor = QColorDialog::getColor(QColor(mCurrentRGB[0],
                                                  mCurrentRGB[1],
                                                  mCurrentRGB[2]), this);
  if (!retColor.isValid())
    return;
  retColor.getRgb(&mCurrentRGB[0], &mCurrentRGB[1], &mCurrentRGB[2]);
  imposeColor(true, true);
}

// Slider changed: change the current value, update color box,
// and emit new color if not dragging or ctrl pressed
void ColorSelector::sliderChanged(int which, int value, bool dragging)
{
  mCurrentRGB[which] = value;
  imposeColor(false, 
	      !dragging || (mHotFlag == HOT_SLIDER_KEYDOWN && mCtrlPressed) ||
	      (mHotFlag == HOT_SLIDER_KEYUP && !mCtrlPressed));

  // Keep track of dragging state AFTER sending signal so first move of a
  // drag is still treated as non-drag and the last one resets the flag
  mDragging = dragging;
}

// Act on a new color
void ColorSelector::imposeColor(bool setSliders, bool emitSignal)
{
  if (setSliders)
    for (int i = 0; i < 3; i++)
      mSliders->setValue(i, mCurrentRGB[i]);

  diaSetWidgetColor(mColorBox, QColor(mCurrentRGB[0], mCurrentRGB[1], 
                                     mCurrentRGB[2]));
  if (emitSignal)
    emit newColor(mCurrentRGB[0], mCurrentRGB[1], mCurrentRGB[2]);
}

void ColorSelector::closeEvent ( QCloseEvent * e )
{
  emit closing();
  e->accept();
}

// watch for ctrl key; emit the key event to pass it on
void ColorSelector::keyPressEvent ( QKeyEvent * e )
{
  bool closing = e->key() == Qt::Key_Escape;
#ifdef Q_OS_MACX
  if (!closing)
    closing = e->key() == Qt::Key_W && (e->modifiers() ==Qt::ControlModifier);
#endif
  if (closing) {
    emit done();
  } else {
    
    if (mHotFlag != NO_HOT_SLIDER && e->key() == mHotKey) {
      mCtrlPressed = true;
      grabKeyboard();
    }
    emit (keyPress(e));
  }

}


void ColorSelector::keyReleaseEvent ( QKeyEvent * e )
{
  if (e->key() == mHotKey) {
    mCtrlPressed = false;
    releaseKeyboard();
  }
  emit (keyRelease(e));
}
