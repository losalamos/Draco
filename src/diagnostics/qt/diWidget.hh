//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/qt/diWidget.hh
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Declaration of draco info widget.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//
#ifndef diwidget_hh
#define diwidget_hh

#include <QWidget>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>

class diWidget : public QWidget
{
    Q_OBJECT

public:
    explicit diWidget(QWidget *parent = 0);
    ~diWidget() { /* empty */ };

private slots:
    /* void on_actionAbout_triggered(); */
    /* void on_actionDialog_with_Layout_triggered(); */
    /* void on_actionEdit_File_triggered(); */

private:
    QGridLayout * layout;
    QLabel      * label1;
    QPushButton * pushbutton1;

    // disable copy construction
    diWidget(diWidget const & rhs);

    // disable assignment
    diWidget & operator=(diWidget const & rhs);
};

#endif // diwidget_hh

//---------------------------------------------------------------------------//
// end of diagnostics/qt/diWidget.hh
//---------------------------------------------------------------------------//
