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

#ifndef CVIEWER_H
#define CVIEWER_H

#include <QWheelEvent>
#include <QGLWidget>

#include "cam.h"
#include <list>

class CViewer : public QGLWidget
{
    Q_OBJECT
public:

    explicit CViewer(QWidget *parent = 0);

    explicit CViewer(const R4R::CView<float>& view, QWidget *parent = 0);

    std::list<std::pair<R4R::vec3f,R4R::rgb> >& get_map() { return m_map; }

signals:
    
public slots:
    
    void update_display(const R4R::CView<float>& view);

    void load_camera(const R4R::CPinholeCam& cam);

protected:

    void initializeGL();

    void resizeGL(int w, int h);

    void paintGL();

    void load_view(const R4R::CView<float>& view);

    virtual void wheelEvent(QWheelEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);

private:

    void translate(float dx, float dy, float dz);

    R4R::CPinholeCam m_cam;
    R4R::CView<float> m_view;                               //!< current view
    std::list<std::pair<R4R::vec3f,R4R::rgb> > m_map;       //!< colored map
    QPoint m_last_point;                                    //!< last clicked point

};

#endif // CVIEWER_H
