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

#include <QFileDialog>

using namespace R4R;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
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


}

void MainWindow::show_image() {

    // convert to gray value
    CGrayValueImage img = CGrayValueImage(m_u);

    // convert to Qt image
    QImage qimg = img.operator QImage();

    // set pixmap
    QPixmap pixmap = QPixmap::fromImage(qimg);

    if(3*qimg.height()>4*qimg.width())
        ui->imgLabel->setPixmap(pixmap.scaledToHeight(ui->imgLabel->height()));
    else
        ui->imgLabel->setPixmap(pixmap.scaledToWidth(ui->imgLabel->width()));

}

void MainWindow::on_actionSave_triggered()
{

}
