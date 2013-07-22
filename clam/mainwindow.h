#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QKeyEvent>

#include "opencv2/opencv.hpp"

#include "params.h"
#include "mtracker.h"
#include "dagg.h"

#include "preferences.h"


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

    void on_actionQuit_triggered();

    void on_actionPreferences_triggered();

    void on_stepButton_clicked();

    void on_pauseButton_clicked();

    void on_playButton_clicked();

    void on_actionOpen_triggered();

    void on_showMemoryUsage_triggered();

    void set_params(CParameters params);

    void on_actionSave_Motion_triggered();

    void on_actionSave_Map_triggered();

    void on_actionSave_Tracks_triggered();

signals:

    void show_memoryUsage();

private:

    Ui::MainWindow *ui;
    Preferences* m_preferences;

    VideoCapture m_cap;
    vector<cv::Mat> m_pyramid;
    QTimer m_timer;
    CParameters m_params;
    CTracker* m_tracker;
    CPinholeCam m_cam;
    std::map<R4R::vec3f,R4R::rgb> m_map;

    void show_image(const QImage& qimg);

};

#endif // MAINWINDOW_H
