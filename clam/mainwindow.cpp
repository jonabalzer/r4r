#include <QFileDialog>
#include <QInputDialog>
#include <QProgressDialog>
#include <QMessageBox>

#include <sys/time.h>
#include <sys/resource.h>

#include <omp.h>

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_preferences(new Preferences),
    m_cap(),
    m_pyramid(),
    m_timer(this),
    m_params(),
    m_tracker(),
    m_cam() {

    ui->setupUi(this);
    this->setFixedSize(this->width(),this->height());

    // connect signals and slots
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(on_stepButton_clicked()));
    connect(this,SIGNAL(show_memoryUsage()),this,SLOT(on_showMemoryUsage_triggered()));
    connect(m_preferences,SIGNAL(params_changed(CParameters)),this,SLOT(set_params(CParameters)));

    // set parameters
    m_preferences->on_applyButton_clicked();

}

MainWindow::~MainWindow() {

    delete m_tracker;
    m_cap.release();
    delete ui;

}

void MainWindow::set_params(CParameters params) {

    m_params = params;
    m_cam = CPinholeCam(m_params.GetDoubleParameter("FU"),
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
    cvtColor(img, img_gray, CV_BGR2GRAY);
    buildPyramid(img_gray,m_pyramid,m_params.GetIntParameter("SCALE"));

    // start measuring time
    double t0, t1;
    t0 = omp_get_wtime();

    // update motion estimates
    m_tracker->Update(pyramid,m_pyramid);

    // add new tracks
    size_t noactive = m_tracker->GetNumberOfActiveTracks();
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
                                                    "(*.mp4);;(*.avi);;(*.mov)");


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
    cvtColor(img, img_gray, CV_BGR2GRAY);
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

    CMotionTracker* tracker = static_cast<CMotionTracker*>(m_tracker);

    tracker->SaveMotion(filename.toStdString().c_str());

}

void MainWindow::on_actionSave_Map_triggered() {

    QString filename = QFileDialog::getSaveFileName(this, tr("Save file..."),
                               ".",
                               tr("(*.ply)"));

    CMotionTracker* tracker = static_cast<CMotionTracker*>(m_tracker);

    tracker->SaveMap(filename.toStdString().c_str());

}
