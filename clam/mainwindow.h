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

    void on_actionClose_triggered();

signals:

    void show_memoryUsage();

private:

    Ui::MainWindow *ui;
    Preferences* m_preferences;

    cv::VideoCapture m_cap;
    vector<cv::Mat> m_pyramid;
    QTimer m_timer;
    CParameters m_params;
    CTracker* m_tracker;
    CPinholeCam m_cam;
    std::map<R4R::vec3f,R4R::rgb> m_map;

    void show_image(const QImage& qimg);

};

#endif // MAINWINDOW_H
