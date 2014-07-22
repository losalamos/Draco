#include "mainwindow.hh"
#include <QMenuBar>

#include "../draco_info.hh"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    // http://qt-project.org/doc/qt-4.8/QMainWindow.html#qt-main-window-framework
    setObjectName("MainWindow");
    
    //setupMenuBar();
    //statusBar()->showMessage(tr("Status: Ok"));

    layout = new QGridLayout(this);
    
    pushbutton1 = new QPushButton("&Ok");
    connect(pushbutton1,SIGNAL(clicked()),this,SLOT(close()));

    // Get the release information from Lib_diagnostics
    QString msg( rtt_diagnostics::DracoInfo().briefReport().c_str() );
    
    label1 = new QLabel(msg);
    label1->setWordWrap(false);
    
    layout->addWidget(label1,0,0);
    layout->addWidget(pushbutton1,1,0,Qt::AlignCenter);
    layout->columnStretch(0);
    
    QWidget *window = new QWidget();
    window->setLayout(layout);
    // window->setWindowTitle("draco_info-gui");
    setCentralWidget(window);

    // Set the window title
    this->setWindowTitle("draco_info-gui");
}

MainWindow::~MainWindow()
{
}

void MainWindow::setupMenuBar()
{
    QMenu *menu = menuBar()->addMenu(tr("&File"));
    menu->addAction(tr("&Quit"), this, SLOT(close()));
}
