#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>
#include <QInputDialog>
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
    m_tracker()
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

    delete m_tracker;

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

    cout << "PARAMS: " << endl;
    cout << m_params << endl;

    // init tracker
    m_tracker = new CSimpleTracker(&m_params);

    // initial detection
    m_tracker->Init(m_pyramid);

    // initinialization of descriptors
    m_tracker->UpdateDescriptors(m_pyramid);

    // draw
    m_tracker->Draw(img);


    show_image(img);

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

    // draw
    m_tracker->Draw(img);
    m_tracker->DrawTails(img,20);

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

    string dir = dirname.toStdString() + string("/");
    string prefix = string("track");

    m_tracker->SaveToFile(dir.c_str(),prefix.c_str());

}

void MainWindow::on_actionPreferences_triggered()
{
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
    items << tr("Identity") << tr("Gradient");

    QString item = QInputDialog::getItem(this, tr("Save Descriptors"),
                                          tr("Descriptor:"), items, 0, false, &ok);


    if(!ok)
        return;


    string name;

    if(item=="Identity")
        name = string("ID");
    else if(item=="Gradient")
        name = string("GRAD");

    QString filename = QFileDialog::getSaveFileName(this,
                                                    tr("Save File"),
                                                    "",
                                                    tr("(*.txt)"));

    m_tracker->SaveDescriptors(filename.toStdString().c_str(),name.c_str(),comment.toStdString().c_str());

}

void MainWindow::on_actionClose_triggered()
{
    m_timer.stop();
    ui->labelImage->clear();
    m_cap.release();

}
