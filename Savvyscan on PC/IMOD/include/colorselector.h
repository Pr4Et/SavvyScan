/*   colorselector.h  -  declarations for colorselector.cpp
 *
 *   Copyright (C) 1995-2002 by  the Regents of the University of 
 *   Colorado.  See implementation file for full copyright notice.
 *
 *   $Id$
 */                                                                           

#include "dialog_frame.h"
//Added by qt3to4:
#include <QTimerEvent>
#include <QKeyEvent>
#include <QCloseEvent>
#include "dllexport.h"

class MultiSlider;
class QFrame;
class ColorSelectorGL;

class DLL_IM_EX ColorSelector : public DialogFrame
{
  Q_OBJECT

 public:
  ColorSelector(QWidget *parent, QString label, int red, int green, int blue, 
                int hotFlag, int hotKey, bool rounded, const char *name = NULL,
                Qt::WindowFlags fl = Qt::Window);
  ~ColorSelector();
  bool hotSliding() {return mDragging;};

 signals:
  void newColor(int r, int g, int b);
  void done();
  void closing();
  void keyPress( QKeyEvent * e );
  void keyRelease( QKeyEvent * e );

  public slots:
    void buttonClicked(int which);
    void sliderChanged(int which, int value, bool dragging);

 protected:
    void closeEvent ( QCloseEvent * e );
    void keyPressEvent ( QKeyEvent * e );
    void keyReleaseEvent ( QKeyEvent * e );

 private:
    void donePressed();
    void restorePressed();
    void qtSelectorPressed();
    void imposeColor(bool setSliders, bool emitSignal);
    bool mCtrlPressed;
    bool mDragging;
    int mHotKey;
    int mHotFlag;
    int mOriginalRGB[3];
    int mCurrentRGB[3];
    MultiSlider *mSliders;
    QFrame *mColorBox;
};
