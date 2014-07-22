#include "mainwindow.hh"

MainWindow::MainWindow(QWidget* /*parent*/)
{
    // http://qt-project.org/doc/qt-4.8/QMainWindow.html#qt-main-window-framework
    setObjectName("MainWindow");

    // Set the window title
    setWindowTitle("draco_info-gui");

    diw = new diWidget(this);
    setCentralWidget(diw);
}
