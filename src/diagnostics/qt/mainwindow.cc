#include "mainwindow.hh"
//#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
    // ui(new Ui::MainWindow)
{
    // ui->setupUi(this);

    // Set the window title
    setWindowTitle("draco_info-gui");

    use_ui_file=false;
    
    if(use_ui_file)
    {
        // // Create a button called with text "Close"
        // // I used the GUI editor on mainwindow.ui to set the signal option
        // // "clicked" to close the application. 
        // ui->pushButton->setText("Close");
        
        // // Manually create a signal/slot communication between the slider and the
        // // progress bar.
        
        // connect(ui->horizontalSlider,SIGNAL(valueChanged(int)),
        //         ui->progressBar,SLOT(setValue(int)));
        
        // connect(ui->horizontalSlider,SIGNAL(valueChanged(int)),
        //         ui->progressBar_2,SLOT(setValue(int)));
    }
    else
    {
        // Play with layouts without the .ui file.

        layout = new QGridLayout;
        
        label1 = new QLabel("Name:");
        text1  = new QLineEdit;
        layout->addWidget(label1,0,0);
        layout->addWidget(text1,0,1);

        QString mStatus("YOLO");
        label2 = new QLabel("/e/qt exists:");
        text2  = new QLineEdit(mStatus);
        layout->addWidget(label2,1,0);
        layout->addWidget(text2,1,1);
        
        pushbutton1 = new QPushButton("OK");
        layout->addWidget(pushbutton1,2,0,1,2); // span 2 cols.
        
        setLayout(layout);
        
        // pressing the button will close this dialog box.
        connect(pushbutton1,SIGNAL(clicked()),
                this,SLOT(close()));
    }
    return;
}

MainWindow::~MainWindow()
{
    // if( use_ui_file )
    // delete ui;
}

// void MainWindow::on_actionAbout_triggered()
// {
//     aboutdialog = new AboutDialog(this);
//     aboutdialog->show();
// }

// void MainWindow::on_actionDialog_with_Layout_triggered()
// {
//     m_dialog_with_layout = new dialog_with_layout(this);
//     m_dialog_with_layout->show();
// }

// void MainWindow::on_actionEdit_File_triggered()
// {
//     editwindow = new EditWindow(this);
//     editwindow->show();
// }
