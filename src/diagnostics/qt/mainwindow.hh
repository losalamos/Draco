#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
/* #include "aboutdialog.h" */
/* #include "dialog_with_layout.h" */
/* #include "editwindow.h" */
#include <QDialog>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
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
    Ui::MainWindow *ui;
    bool use_ui_file;
    /* AboutDialog *aboutdialog; */
    /* dialog_with_layout * m_dialog_with_layout; */
    /* EditWindow * editwindow; */
    QGridLayout * layout;
    QLabel      * label1;
    QLineEdit   * text1;
    QLabel      * label2;
    QLineEdit   * text2;
    QPushButton * pushbutton1;
};

#endif // MAINWINDOW_H
