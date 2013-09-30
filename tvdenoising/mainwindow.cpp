/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "nabla.h"

#include <QFileDialog>

using namespace R4R;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_f(),
    m_A(),
    m_nabla(),
    m_u()
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionExit_triggered()
{
    QApplication::exit();
}

void MainWindow::on_actionOpen_triggered()
{

    // get filename
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Open image file...",
                                                    ".",
                                                    "(*.png);;(*.bmp);;(*.jpg)");

    // open qimag
    QImage qimg(filename);

    // first convert into a grayvalue image
    CGrayValueImage myimg(qimg);

    // then into a floating point array
    m_f = matf(myimg);

    // set other member variables
    m_u = m_f.Clone();
    show_image();

    // compute derivatives for new image size
    CImageDenoising der(m_u.NCols(),m_u.NRows());
    der.ComputeGradientOperator(m_nabla);
    der.ComputeJacobian(m_A);

}

void MainWindow::show_image() {

    if(m_u.NElems()==0)
        return;

    // convert to gray value
    CGrayValueImage img = CGrayValueImage(m_u);

    // convert to Qt image
    QImage qimg = img.operator QImage();

    // set pixmap
    QPixmap pixmap = QPixmap::fromImage(qimg);

    if(3*qimg.height()>4*qimg.width() || qimg.height()==qimg.width())
        ui->imgLabel->setPixmap(pixmap.scaledToHeight(ui->imgLabel->height()));
    else
        ui->imgLabel->setPixmap(pixmap.scaledToWidth(ui->imgLabel->width()));

}

void MainWindow::on_actionSave_triggered()
{

    if(m_u.NElems()==0)
        return;

    // convert to gray value
    CGrayValueImage img = CGrayValueImage(m_u);

    // convert to Qt image
    QImage qimg = img.operator QImage();

    // get filename
    QString filename = QFileDialog::getSaveFileName(this,
                                                    "Save as...",
                                                    ".",
                                                    "(*.png);;(*.bmp);;(*.jpg)");

    qimg.save(filename);

}

void MainWindow::on_actionReload_triggered()
{

    // set other member variables
    m_u = m_f.Clone();
    show_image();

}

void MainWindow::on_iterButton_clicked()
{

    // nothing to do if there is no image
    if(m_f.NElems()==0)
        return;

    // create linear solver
    CConjugateGradientMethodLeastSquares<smatf,float> linsolver(CPreconditioner<smatf,float>(),
                                                                ui->linIterSpinBox->value(),
                                                                1e-20,
                                                                true);
    // get reshaped of f and u
    vecf u = vecf(m_u);
    vecf f = vecf(m_f);

    // create SB solver
    CSplitBregman<smatf,float> solver(m_A,
                                      m_nabla,
                                      f,
                                      u,
                                      linsolver,
                                      ui->muEdit->text().toFloat(),
                                      ui->lambdaEdit->text().toFloat(),
                                      ui->epsEdit->text().toFloat());

    // do one iteration step
    solver.Iterate(ui->niterSpinBox->value());

    // convert back to image
    m_u = matf(m_u.NRows(),m_u.NCols(),u);

    // show errors
    vector<double>& error = solver.GetTotalError();
    vector<double>& constraint = solver.GetConstraintViolation();
    ui->errorLcdNumber->display(error.back());
    ui->constraintLcdNumber->display(constraint.back());

    // show image
    show_image();

}
