#include "dialog_frame.h"
//Added by qt3to4:
#include <QKeyEvent>
#include <QCloseEvent>
class QPushButton;
class QCheckBox;
class DockingDialog;

class BeadFixer2 : public DialogFrame
{
  Q_OBJECT

 public:
  BeadFixer2(DockingDialog *parent, const char *name = NULL);
  ~BeadFixer2() {};
  void reread(int which);

  public slots:
  void buttonPressed(int which);
  void nextGap();
  void openFile();
  void rereadFile() {reread(0);};
  void nextLocal() {reread(1);};
  void nextRes();
  void backUp();
  void movePoint();
  void undoMove();
  void clearList();
  void onceToggled(bool state);
  void topCloseEvent ( QCloseEvent * e );
  void topChangeEvent(QEvent *e);

 protected:
  void keyPressEvent ( QKeyEvent * e );
  void keyReleaseEvent ( QKeyEvent * e );

 private:
  int foundgap(int obj, int cont, int ipt, int before);
  void clearExtraObj();
  DockingDialog *mTopWin;
  QPushButton *rereadBut;
  QPushButton *nextLocalBut;
  QPushButton *nextResBut;
  QPushButton *backUpBut;
  QPushButton *movePointBut;
  QPushButton *undoMoveBut;
  QPushButton *clearListBut;
};
