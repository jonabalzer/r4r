#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QKeyEvent>

#include "opencv2/opencv.hpp"

#include "params.h"
#include "stracker.h"

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

private:
    Ui::MainWindow *ui;
    cv::VideoCapture m_cap;
    std::vector<cv::Mat> m_pyramid;
    QTimer m_timer;
    R4R::CParameters m_params;
    std::vector<R4R::CTracker*> m_trackers;



    void show_image(const cv::Mat& img);

};

#endif // MAINWINDOW_H
