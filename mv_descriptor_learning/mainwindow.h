#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QKeyEvent>
#include "preferences.h"

#include "opencv2/opencv.hpp"

#include "params.h"
#include "stracker.h"

using namespace cv;
using namespace std;
using namespace R4R;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_actionExit_triggered();

    void on_actionOpen_triggered();

    void on_stepButton_clicked();

    void on_playButton_clicked();

    void on_pauseButton_clicked();

    void on_actionSave_Tracks_triggered();

    void on_actionPreferences_triggered();

    void set_params(CParameters params) { m_params = params; };

    void on_actionSave_Descriptors_triggered();

    void on_actionClose_triggered();

private:
    Ui::MainWindow* ui;

    Preferences* m_preferences;

    VideoCapture m_cap;
    vector<cv::Mat> m_pyramid;
    QTimer m_timer;
    CParameters m_params;
    CTracker* m_tracker;

    void show_image(const Mat& img);

};

#endif // MAINWINDOW_H
