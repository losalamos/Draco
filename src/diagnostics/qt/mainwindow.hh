#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "diWidget.hh"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow(){};

private slots:
    // None

private:
    // Forms
    Ui::MainWindow *ui;
    
    // Widgets
    diWidget    * diw;
};

#endif // MAINWINDOW_H
