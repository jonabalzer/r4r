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

#include "viewer.h"

using namespace cv;
using namespace std;
using namespace R4R;

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_preferences(new Preferences),
    m_cap(),
    m_pyramid(),
    m_timer(this),
    m_params(),
    m_tracker(),
    m_cam(),
    m_map() {

    ui->setupUi(this);
    this->setFixedSize(this->width(),this->height());

    // connect signals and slots
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(on_stepButton_clicked()));
    connect(this,SIGNAL(show_memoryUsage()),this,SLOT(on_showMemoryUsage_triggered()));
    connect(m_preferences,SIGNAL(params_changed(CParameters)),this,SLOT(set_params(CParameters)));

    // set parameters
    m_preferences->on_applyButton_clicked();

    // delete the one we already have and replace it
    /*CView<float> view(m_cam);
    ui->openGlWidget = new CViewer(view,this);
    ui->openGlWidget->setObjectName(QString::fromUtf8("openGlWidget"));
    ui->openGlWidget->setGeometry(QRect(640, 45, 640, 480));*/

}

MainWindow::~MainWindow() {

    delete m_tracker;
    m_cap.release();
    delete ui;

}

void MainWindow::set_params(CParameters params) {

    m_params = params;
    m_cam = CPinholeCam(m_params.GetIntParameter("SU"),
                        m_params.GetIntParameter("SV"),
                        m_params.GetDoubleParameter("FU"),
                        m_params.GetDoubleParameter("FV"),
                        m_params.GetDoubleParameter("CU"),
                        m_params.GetDoubleParameter("CV"));


}

void MainWindow::on_actionQuit_triggered() {

     QApplication::exit();

}

void MainWindow::show_image(const QImage& qimg) {

    QPixmap pixmap = QPixmap::fromImage(qimg.rgbSwapped());

    if(3*qimg.height()>4*qimg.width())
        ui->labelImage->setPixmap(pixmap.scaledToHeight(ui->labelImage->height()));
    else
        ui->labelImage->setPixmap(pixmap.scaledToWidth(ui->labelImage->width()));

}

void MainWindow::on_actionPreferences_triggered() {

    m_timer.stop();
    m_preferences->show();

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

    // get map and color it
    if(ui->renderCheckBox->isChecked()) {

        CMotionTracker* tracker = dynamic_cast<CMotionTracker*>(m_tracker);
        CViewer* viewer = dynamic_cast<CViewer* >(ui->openGlWidget);
        list<pair<vec2f,vec3f> >& newpts = tracker->GetMap();
        list<pair<vec3f,rgb> >& allpts = viewer->get_map();

        // color them and send them to the viewer
        list<pair<vec2f,vec3f> >::iterator it;
        for(it=newpts.begin(); it!=newpts.end(); it++) {

            CVector<size_t,2> p = CVector<size_t,2>(it->first);

            rgb red = { img.at<Vec3b>(p.Get(1),p.Get(0))[2], img.at<Vec3b>(p.Get(1),p.Get(0))[1], img.at<Vec3b>(p.Get(1),p.Get(0))[0] };
            allpts.push_back(pair<vec3f,rgb>(it->second,red));

         }

        // update view
        CView<float> view = tracker->GetLatestView();
        viewer->update_display(view);

    }

    // update descriptors
    m_tracker->UpdateDescriptors(m_pyramid);

    // check validity of tracks
    m_tracker->Clean(pyramid,m_pyramid);
    //trackers[i]->DeleteInvalidTracks();

    // add new tracks
    size_t noactive = m_tracker->GetNumberOfActiveTracks();
    if(noactive<(size_t)m_params.GetIntParameter("MIN_NO_FEATURES"))
         m_tracker->AddTracklets(m_pyramid);

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

void MainWindow::on_pauseButton_clicked() {

     m_timer.stop();

}

void MainWindow::on_playButton_clicked() {

    if(m_cap.isOpened())
        m_timer.start(35);
    else
        QMessageBox::critical(this,"Error","Could not open source.");

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
                                                    "(*.mov);;(*.mp4);;(*.avi)");


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

    } else
        return;


    // convert
    cvtColor(img, img_gray, COLOR_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    // init tracker
    m_tracker = new CMotionTracker(&m_params,m_cam);

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

void MainWindow::on_showMemoryUsage_triggered() {

    rusage usage;
    getrusage(RUSAGE_SELF,&usage);
    double mb = (double)usage.ru_maxrss/1024;
    ui->memLcdNumber->display(mb);

}


void MainWindow::on_actionSave_Motion_triggered() {

    QString filename = QFileDialog::getSaveFileName(this, tr("Save file..."),
                               ".",
                               tr("(*.txt)"));

    CMotionTracker* tracker = dynamic_cast<CMotionTracker*>(m_tracker);

    list<vecf> motion = tracker->GetMotion();

    ofstream out(filename.toStdString().c_str());

    if(!out) {

        QMessageBox::critical(this,"Error","Could not save file.");
        return;

     }

    list<vecf>::iterator it;

    for(it=motion.begin(); it!=motion.end(); it++) {

        vecf m = (*it);
        m.Transpose();
        out << m << endl;

    }

    out.close();

}

void MainWindow::on_actionSave_Map_triggered() {

    QString filename = QFileDialog::getSaveFileName(this, tr("Save file..."),
                               ".",
                               tr("(*.ply)"));


    CViewer* viewer = dynamic_cast<CViewer* >(ui->openGlWidget);
    list<pair<vec3f,rgb> >& pts = viewer->get_map();

    list<pair<vec3f,rgb> >::iterator it;

    ofstream out(filename.toStdString().c_str());

    if(!out) {

        QMessageBox::critical(this,"Error","Could not save file.");
        return;

     }

    out << "ply" << endl;
    out <<  "format ascii 1.0" << endl;
    out <<  "comment" << endl;
    out << "element vertex " << pts.size() << endl;
    out << "property float32 x" << endl;
    out << "property float32 y" << endl;
    out << "property float32 z" << endl;
    out << "property uchar red" << endl;
    out << "property uchar green" << endl;
    out << "property uchar blue" << endl;
    out << "end_header" << endl;

    for(it=pts.begin(); it!=pts.end(); it++)
        out << it->first.Get(0) << " " << it->first.Get(1) << " " << it->first.Get(2) << " " << (unsigned int)it->second.Get(0) << " " <<  (unsigned int)it->second.Get(1) << " " <<  (unsigned int)it->second.Get(2) << endl;

    out.close();

}

void MainWindow::on_actionSave_Tracks_triggered()
{

    QString dir = QFileDialog::getExistingDirectory(this,tr("Choose directory..."),".",QFileDialog::ShowDirsOnly);

    bool error = m_tracker->SaveToFile(dir.toStdString().c_str(),"/track");

    if(error)
         QMessageBox::critical(this,"Error","Could not save tracks.");

}

void MainWindow::on_actionClose_triggered()
{

    m_timer.stop();
    ui->labelImage->clear();
    delete m_tracker;
    m_cap.release();

    CViewer* viewer = dynamic_cast<CViewer* >(ui->openGlWidget);
    viewer->get_map().clear();
    viewer->repaint();

    emit show_memoryUsage();

    ui->frameLcdNumber->display(0);

}
