//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/qt/main.cc
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2014, 17:05 pm
 * \brief  Main program for Gui version of draco info.
 * \note   Copyright (C) 2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
// end of diagnostics/qt/main.cc
//---------------------------------------------------------------------------//
