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

#include <QFileDialog>
#include <QInputDialog>
#include <QProgressDialog>
#include <QMessageBox>

#include <sys/time.h>
#include <sys/resource.h>

#include <omp.h>

#include "mainwindow.h"
#include "ui_mainwindow.h"

using namespace cv;
using namespace std;
using namespace R4R;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_preferences(new Preferences),
    m_cap(),
    m_pyramid(),
    m_timer(this),
    m_params(),
    m_tracker()
{

    // set gui
    ui->setupUi(this);
    this->setFixedSize(this->width(),this->height());

    // connect signals and slots
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(on_stepButton_clicked()));
    connect(m_preferences,SIGNAL(params_changed(CParameters)),this,SLOT(set_params(CParameters)));
    connect(this,SIGNAL(show_memoryUsage()),this,SLOT(on_showMemoryUsage_triggered()));

    // set parameters
    m_preferences->on_applyButton_clicked();

}

MainWindow::~MainWindow()
{

    delete m_tracker;
    m_cap.release();
    delete ui;

}

void MainWindow::on_actionExit_triggered()
{
     QApplication::exit();
}

void MainWindow::show_image(const QImage& qimg) {

    QPixmap pixmap = QPixmap::fromImage(qimg.rgbSwapped());

    if(3*qimg.height()>4*qimg.width())
        ui->labelImage->setPixmap(pixmap.scaledToHeight(ui->labelImage->height()));
    else
        ui->labelImage->setPixmap(pixmap.scaledToWidth(ui->labelImage->width()));

}


void MainWindow::on_actionOpen_triggered()
{

    // check if a capture is already running
    m_timer.stop();
    if(m_cap.isOpened())
        return;

    // declare variable for first frame
    Mat img, img_gray;

    // get file name
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Open video file...",
                                                    ".",
                                                    "(*.mp4);;(*.avi);;(*.mov);;(*.png);;(*.jpg);;(*.bmp);;(*.ppm);;(*.pgm)");


    // check if it is an image or video
    if(filename.endsWith(".mov")||
       filename.endsWith(".avi")||
       filename.endsWith(".mp4")) {

        m_cap = VideoCapture(filename.toStdString().c_str());

        // check if we can properly open the video
        if(!m_cap.isOpened() && !filename.isEmpty()) {
            QMessageBox::critical(this,"Error","Could not open video stream. Check path and whether required codecs are installed.");
            return;
        }

        // grab the first frame
        if(!m_cap.grab()) {

            ui->statusBar->showMessage("Failed to get frame.");
            return;

        }

        m_cap.retrieve(img);

    } else if(filename.endsWith(".png")||
              filename.endsWith(".jpg")||
              filename.endsWith(".bmp")||
              filename.endsWith(".ppm")||
              filename.endsWith(".pgm")) {

        // open image
        img = imread(filename.toStdString().c_str());

        if(img.rows==0 || img.cols==0) {
            QMessageBox::critical(this,"Error","Could not open image.");
            return;
        }

    }
    else return;

    // convert
    cvtColor(img, img_gray, COLOR_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    // init tracker
    m_tracker = new CSimpleTracker(&m_params);

    // initial detection
    m_tracker->Init(m_pyramid);

    // initinialization of descriptors
    m_tracker->UpdateDescriptors(m_pyramid);

    // draw
    QImage qimg(img.data,img.cols,img.rows,QImage::Format_RGB888);
    m_tracker->Draw(qimg,1);

    // display
    show_image(qimg);
    emit show_memoryUsage();

}

void MainWindow::on_stepButton_clicked()
{

    if(!m_cap.isOpened()) {
        m_timer.stop();
        return;
    }

    Mat img, img_gray;
    vector<Mat> pyramid;

    for(size_t s=0; s<m_pyramid.size(); s++)
        pyramid.push_back(m_pyramid[s].clone());

    m_pyramid.clear();

    if(!m_cap.grab()) {

        m_timer.stop();
        ui->statusBar->showMessage("Failed to get frame.");
        return;

    }

    m_cap.retrieve(img);
    cvtColor(img, img_gray, COLOR_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    // start measuring time
    double t0, t1;
    t0 = omp_get_wtime();

    // update motion estimates
    m_tracker->Update(pyramid,m_pyramid);

    // add new tracks
    size_t noactive = m_tracker->GetNumberOfActiveTracks();
    cout << "Number of active tracklets: " << noactive << endl;
    if(noactive<(size_t)m_params.GetIntParameter("MIN_NO_FEATURES"))
         m_tracker->AddTracklets(m_pyramid);

    // update descriptors
    m_tracker->UpdateDescriptors(m_pyramid);

    // check validity of tracks
    m_tracker->Clean(pyramid,m_pyramid);
    //trackers[i]->DeleteInvalidTracks();

    // compute and display framerate
    t1 = omp_get_wtime();
    double fps = 1.0/(t1-t0);
    ui->speedLcdNumber->display(fps);

    // draw trails
    QImage qimg(img.data,img.cols,img.rows,QImage::Format_RGB888);
    m_tracker->Draw(qimg,5);
    show_image(qimg);

    // display frame no and mem usage
    ui->frameLcdNumber->display((int)m_tracker->GetTime());
    emit show_memoryUsage();

}

void MainWindow::on_showMemoryUsage_triggered() {

    rusage usage;
    getrusage(RUSAGE_SELF,&usage);
    double mb = (double)usage.ru_maxrss/1024;
    ui->memLcdNumber->display(mb);

}

void MainWindow::on_playButton_clicked() {

    if(m_cap.isOpened())
        m_timer.start(35);
    else
        QMessageBox::critical(this,"Error","Could not open source.");

}

void MainWindow::on_pauseButton_clicked() {

     m_timer.stop();

}

void MainWindow::on_actionSave_Tracks_triggered()
{
    m_timer.stop();

    QString dirname = QFileDialog::getExistingDirectory(this,tr("Choose folder..."), ".");

    string dir = dirname.toStdString() + string("/");
    string prefix = string("track");

    m_tracker->SaveToFile(dir.c_str(),prefix.c_str());

}

void MainWindow::on_actionPreferences_triggered() {

    m_timer.stop();
    m_preferences->show();

}

void MainWindow::on_actionSave_Descriptors_triggered()
{

    m_timer.stop();

    bool ok;
    QString comment = QInputDialog::getText(this,
                                            "Save Descriptors",
                                            "Comment:", QLineEdit::Normal,
                                            "",
                                            &ok);

    if (!ok)
        return;

    QStringList items;
    items << tr("Identity") << tr("Gradient") << tr("HoG");

    QString item = QInputDialog::getItem(this, tr("Save Descriptors"),
                                          tr("Descriptor:"), items, 0, false, &ok);


    if(!ok)
        return;

    string name;

    if(item=="Identity")
        name = string("ID");
    else if(item=="Gradient")
        name = string("GRAD");
    else if(item=="HoG")
        name = string("HOG");


    // now save file
    QString filename = QFileDialog::getSaveFileName(this,
                                                    tr("Save File"),
                                                    "",
                                                    tr("(*.txt)"));

    if(filename.isEmpty())
        return;

    // create aggregator
    CDescriptorAggregator<matf>* aggregator;

    switch(m_params.GetIntParameter("AGGREGATOR")) {

    case 1:
        cout << "Init frame." << endl;
        aggregator = new CInitFrameAggregator<matf>(m_tracker,name.c_str());
        break;
    case 2:
        aggregator = new CSubsampleAggregator<matf>(m_tracker,name.c_str(),m_params.GetIntParameter("AGG_DS"));
        break;

    case 3:
        aggregator = new CSplineInterpolationAggregator<matf>(m_tracker,name.c_str(),10,3);
        break;

    case 4:
        aggregator = new CMeanAggregator<matf>(m_tracker,name.c_str());
        break;

    default:

        aggregator = new CDescriptorAggregator<matf>(m_tracker,name.c_str());

        break;
    }

    // aggregate
    aggregator->Aggregate();
    list<imfeature> feats = aggregator->Get();

    imfeature::SaveToFile(filename.toStdString().c_str(),feats,comment.toStdString().c_str());

    delete aggregator;

}

void MainWindow::on_actionClose_triggered()
{

    m_timer.stop();
    ui->labelImage->clear();
    delete m_tracker;
    m_cap.release();

    emit show_memoryUsage();

    ui->frameLcdNumber->display(0);


}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::information(this,"About","More info: http://vision.ucla.edu");
}

