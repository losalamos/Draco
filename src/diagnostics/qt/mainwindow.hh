#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
//#include <QDialog>
#include <QGridLayout>
#include <QLabel>
//#include <QLineEdit>
#include <QPushButton>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    /* void on_actionAbout_triggered(); */
    /* void on_actionDialog_with_Layout_triggered(); */
    /* void on_actionEdit_File_triggered(); */

private:
    QGridLayout * layout;
    QLabel      * label1;
    QPushButton * pushbutton1;

    QMenu *mainWindowMenu;
    
    void setupMenuBar();
};

#endif // MAINWINDOW_H
