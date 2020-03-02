//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   diagnostics/qt/mainwindow.cc
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Implementation for draco info main Qt window.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
#include "mainwindow.hh"
#include <QStatusBar>
//#include "ui_mainwindow.h"

// Compiling with high warning levels will produce this warning:
// warning: base class 'class Ui_MainWindow' has a non-virtual destructor [-Weffc++]
//     class MainWindow: public Ui_MainWindow {};
// This cannot be fixed because ui_mainwindow.h is generated automatically by Qt.

//----------------------------------------------------------------------------//
//! Constructor
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
      //ui(new Ui::MainWindow),
      diw(new diWidget(this)) {
  //ui->setupUi(this);
  setCentralWidget(diw);
  statusBar()->showMessage(tr("Status Bar"));
}

//----------------------------------------------------------------------------//
// end of diagnostics/qt/mainwindow.cc
//----------------------------------------------------------------------------//
