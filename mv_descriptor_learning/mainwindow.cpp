#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>
#include <QProgressDialog>
#include <QMessageBox>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_preferences(new Preferences),
    m_cap(),
    m_pyramid(),
    m_timer(this),
    m_params(),
    m_trackers()
{

    // set gui
    ui->setupUi(this);

    // connect signals and slots
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(on_stepButton_clicked()));
    connect(m_preferences,SIGNAL(params_changed(CParameters)),this,SLOT(set_params(CParameters)));

    // set parameters
    m_preferences->on_applyButton_clicked();

}

MainWindow::~MainWindow()
{

    for(size_t i=0; i<m_trackers.size(); i++)
        delete m_trackers[i];

    m_cap.release();

    delete ui;

}

void MainWindow::on_actionExit_triggered()
{
     QApplication::exit();
}

void MainWindow::show_image(const Mat& img) {

    QImage qimg(img.data,img.cols,img.rows,QImage::Format_RGB888);
    ui->labelImage->setPixmap(QPixmap::fromImage(qimg.rgbSwapped()));

}


void MainWindow::on_actionOpen_triggered()
{

    // check if a capture is already running
    m_timer.stop();
    if(m_cap.isOpened())
        return;

    // open stream
    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Open video file...",
                                                    ".",
                                                    "(*.mov);;(*.avi);;(*.mp4)");

    m_cap = VideoCapture(filename.toStdString().c_str());

    if(!m_cap.isOpened() && !filename.isEmpty()) {
        QMessageBox::critical(this,"Error","Could not open file. Check path and whether required codecs are installed.");
        return;
    }

    // get initial frame
    Mat img, img_gray;

    if(!m_cap.grab()) {

        ui->statusBar->showMessage("Failed to get frame.");
        return;

    }

    m_cap.retrieve(img);
    cvtColor(img, img_gray, CV_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    // init tracker
    CSimpleTracker* tracker = new CSimpleTracker(m_params);
    m_trackers.push_back(tracker);

    for(size_t i = 0; i<m_trackers.size(); i++) {

         // initial detection
         m_trackers[i]->Init(m_pyramid);

         // initinialization of descriptors
         m_trackers[i]->UpdateDescriptors(m_pyramid);

         // draw
         m_trackers[i]->Draw(img);

    }

    show_image(img);

}

void MainWindow::on_stepButton_clicked()
{

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
    cvtColor(img, img_gray, CV_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    for(size_t i = 0; i<m_trackers.size(); i++) {

        // update motion estimates
        m_trackers[i]->Update(pyramid,m_pyramid);

        // add new tracks
        size_t noactive = m_trackers[i]->GetNumberOfActiveTracks();
        if(noactive<(size_t)m_params.GetIntParameter("MIN_NO_FEATURES"))
            m_trackers[i]->AddTracklets(m_pyramid);

        // update descriptors
        m_trackers[i]->UpdateDescriptors(m_pyramid);

        // check validity of tracks
        m_trackers[i]->Clean(pyramid,m_pyramid);
        //trackers[i]->DeleteInvalidTracks();

        // draw
        m_trackers[i]->Draw(img);
        m_trackers[i]->DrawTails(img,20);

    }

    show_image(img);

}

void MainWindow::on_playButton_clicked()
{

    if(m_cap.isOpened())
        m_timer.start(35);
    else
        QMessageBox::critical(this,"Error","Could not open source. Make sure the Kinect sensor is connected to your computer and drivers are working properly.");

}

void MainWindow::on_pauseButton_clicked()
{
     m_timer.stop();
}

void MainWindow::on_actionSave_Tracks_triggered()
{
    m_timer.stop();

    QString dirname = QFileDialog::getExistingDirectory(this,tr("Choose folder..."), ".");

    for(size_t i = 0; i<m_trackers.size(); i++) {

        stringstream ss;
        ss << i;

        string dir = dirname.toStdString() + string("/");
        string prefix = string("tracker") + ss.str() + string("_");

        //m_trackers[i]->SaveToFile(dir.c_str(),prefix.c_str());
        m_trackers[i]->SaveToFileBlockwise(dir.c_str(),61*61,1);

    }


}

void MainWindow::on_actionPreferemces_triggered()
{
    m_preferences->show();
}
