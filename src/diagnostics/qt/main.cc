//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/qt/main.cc
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Main program for Gui version of draco info.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
#include "mainwindow.hh"
#include <QApplication>

int main(int argc, char *argv[]) {
  // http://qt-project.org/doc/qt-5/qtwidgets-mainwindows-mainwindow-main-cpp.html
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
