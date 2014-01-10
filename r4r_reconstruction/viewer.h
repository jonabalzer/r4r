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
#include "trimesh.h"

namespace R4R {

class CViewer:public QGLWidget {

    Q_OBJECT

public:

    //! Constructor.
    explicit CViewer(const CView<double>& view, QWidget* parent = 0);

    //! Computes depth map for a view.
    matf getDepthMap(const CView<double>& view);

signals:

public slots:

    //! A slot that processes signals that indicate a change in view point.
    void updateView(const R4R::CView<double>& view) { m_view = view; loadView(view); updateGL(); }

    //! Use to trigger changes of the color flag from outside.
    void onShowErrorChanged(int newstate) { m_show_color = bool(newstate);  updateGL(); }

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

protected:

    CView<double> m_view;                  //!< camera view
    double m_znear;                        //!< near clipping plane
    double m_zfar;                         //!< far clipping plane
    QPoint m_last_point;                   //!< auxiliary variable to store mouse pointer locations
    vec3 m_center;                         //!< center of camera rotations
    bool m_show_color;                     //!< color flag

    //! Sends a view to OpenGL.
    virtual void loadView(const CView<double>& view);

    //! Updates projection matrix (e.e. if clip depths changed).
    void loadProjectionMatrix();

};


class CTriMeshViewer : public CViewer {

    Q_OBJECT

public:

    //! Constructor.
    explicit CTriMeshViewer(const R4R::CView<double>& view);

    //! Updates mesh and recompute its bounding box.
    void setMesh(CTriangleMesh* mesh);

    /*! \brief Triggers update of bounding box without changing the pointer to the mesh.
     *
     * This comes in handy e.g. when new polygons are added to the mesh or the mesh is deformed.
     */
    void updateBoundingBox();

public slots:

protected:

    //! \copydoc CViewer::paintGL()
    void paintGL();

    /*! \brief Sends a view to OpenGL.
     *
     * This is overloaded here because if there is contents, we also need to update the clip depths
     * everytime the vantage point changes, and hence also reload the intrinsics.
     *
     */
    void loadView(const CView<double>& view);

    //! Handling of mouse release event.
    virtual void mouseReleaseEvent(QMouseEvent* event);

private:

    //! Update clip depth based on the bounding box approximation of the current mesh.
    void updateClipDepth(const CView<double>& view);

    CTriangleMesh* m_mesh;                      //!< pointer to a triangle mesh
    CBoundingBox<float> m_bbox;                 //!< bounding box of the mesh

};

}

#endif

#endif // VIEWER_H
