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

#ifndef VIEWER_H
#define VIEWER_H

#ifdef QT_GUI_LIB

#include <QGLWidget>
#include <QWheelEvent>

#include "cam.h"

class CViewer:public QGLWidget {

    Q_OBJECT

public:

    //! Constructor.
    explicit CViewer(const R4R::CView<double>& view, QWidget* parent = 0);

signals:

public slots:

    //! A slot that processes signals that indicate a change in view point.
    void updateView(const R4R::CView<double>& view) { m_view = view; loadView(view); updateGL(); }

protected:

    //! Initializes OpenGL context.
    void initializeGL();

    //! Triggers drawing of the scene.
    void paintGL();

    //! Translating wheel event into zooming.
    virtual void wheelEvent(QWheelEvent* event);

    //! Handling of mouse clicks.
    virtual void mousePressEvent(QMouseEvent* event);

    //! Handling of mouse drags.
    virtual void mouseMoveEvent(QMouseEvent* event);

    //! Handling of keyboard inputs.
    virtual void keyPressEvent(QKeyEvent* event);

private:

    R4R::CView<double> m_view;                  //!< camera view
    double m_znear;                             //!< near clipping plane
    double m_zfar;                              //!< far clipping plane
    QPoint m_last_point;                        //!< auxiliary variable to store mouse pointer locations
    R4R::vec3 m_center;                         //!< center of camera rotations

protected:

    //! Sends a view to OpenGL.
    void loadView(const R4R::CView<double>& view);

};

#endif

#endif // VIEWER_H
