#ifndef diwidget_hh
#define diwidget_hh

#include <QWidget>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>

class diWidget : public QWidget
{
    Q_OBJECT

public:
    explicit diWidget(QWidget *parent = 0);
    ~diWidget(){};

private slots:
    /* void on_actionAbout_triggered(); */
    /* void on_actionDialog_with_Layout_triggered(); */
    /* void on_actionEdit_File_triggered(); */

private:
    QGridLayout * layout;
    QLabel      * label1;
    QPushButton * pushbutton1;
};

#endif // diwidget_hh
