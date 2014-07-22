#include "mainwindow.hh"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("draco_info-gui");
    app.setOrganizationName("LANL CCS-2");
    MainWindow *mainWin = new MainWindow;
    mainWin->show();
    return app.exec();
}
