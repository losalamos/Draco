//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/qt/mainwindow.hh
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Declarations for draco info main Qt window.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
#ifndef diagnostics_qt_mainwindow_hh
#define diagnostics_qt_mainwindow_hh

#include "diWidget.hh"
#include <QMainWindow>

//namespace Ui
//{
//   class MainWindow;
//}

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow(){};

private slots:
  // None

private:
  // disable copy construction
  MainWindow(MainWindow const &rhs);

  // disable assignment
  MainWindow &operator=(MainWindow const &rhs);

  // Forms
  //Ui::MainWindow *ui;

  // Widgets
  diWidget *diw;
};

#endif // diagnostics_qt_mainwindow_hh

//---------------------------------------------------------------------------//
// end of diagnostics/qt/mainwindow.hh
//---------------------------------------------------------------------------//
