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

#ifndef R4RVIEWER_H
#define R4RVIEWER_H

#ifdef QT_GUI_LIB

#include <QGLWidget>
#include <QWheelEvent>

#include "cam.h"
#include "trimesh.h"
#include "pcl.h"

namespace R4R {

class CViewer:public QGLWidget {

    Q_OBJECT

public:

    //! Constructor.
    explicit CViewer(const CView<float>& view, QWidget* parent = 0);

    //! Computes depth map for a view.
    CDenseArray<float> getDepthMap(const CView<float>& view);

    //! Method stump.
    virtual void updateBoundingBox() {}

signals:

public slots:

    //! A slot that processes signals that indicate a change in view point.
    virtual void updateView(const R4R::CView<float>& view) { m_view = view; loadView(view); updateGL(); }

    //! Use to trigger changes of the color flag from outside.
    void onShowErrorChanged(int newstate) { m_show_color = bool(newstate);  updateGL(); }

protected:

    //! Initializes OpenGL context.
    void initializeGL();

    //! Triggers drawing of the scene.
    void paintGL();

    //! Translating wheel event into zooming.
    void wheelEvent(QWheelEvent* event);

    //! Handling of mouse clicks.
    void mousePressEvent(QMouseEvent* event);

    //! Handling of mouse drags.
    void mouseMoveEvent(QMouseEvent* event);

    //! Handling of keyboard inputs.
    void keyPressEvent(QKeyEvent* event);

    //! Method stump.
    virtual void updateClipDepth(const CView<float>& view, double tolerance = 1.0) {}

protected:

    CView<float> m_view;                   //!< camera view
    float m_znear;                         //!< near clipping plane
    float m_zfar;                          //!< far clipping plane
    QPoint m_last_point;                   //!< auxiliary variable to store mouse pointer locations
    CVector<float,3> m_center;             //!< center of camera rotations
    CBoundingBox<float> m_bbox;            //!< bounding box of the scene
    bool m_show_color;                     //!< color flag

    //! Sends a view to OpenGL.
    void loadView(const CView<float>& view);

    //! Updates projection matrix (e.e. if clip depths changed).
    void loadProjectionMatrix();


};

#ifdef HAVE_OM

class CTriMeshViewer:public CViewer {

    Q_OBJECT

public:

    //! Constructor.
    explicit CTriMeshViewer(const R4R::CView<float>& view, const CTriangleMesh* mesh = nullptr, QWidget* parent = nullptr);

    //! Updates mesh and recompute its bounding box.
    void setMesh(const CTriangleMesh* mesh);

    /*! \brief Triggers update of bounding box without changing the pointer to the mesh.
     *
     * This comes in handy e.g. when new polygons are added to the mesh or the mesh is deformed.
     */
    void updateBoundingBox();

    //! Assigns to each pixel in the result the handle of the face it sees.
    CDenseArray<int> getFaceMap(const CView<float>& view);

public slots:

    /*! \copybrief CViewer::updateView(const R4R::CView<double>&)
     *
     * This needs to be overridden because we have to adjust the clip depths if the
     * vantage point changes.
     *
     */
    void updateView(const R4R::CView<float>& view);

protected:

    //! \copydoc CViewer::paintGL()
    void paintGL();

    //! Handling of mouse release event.
    void mouseReleaseEvent(QMouseEvent* event);

    //! Update clip depth based on the bounding box approximation of the current mesh.
    void updateClipDepth(const CView<float>& view, float tolerance = 1.0);

private:

    const CTriangleMesh* m_mesh;                //!< pointer to the point cloud

};

#endif

class CPointCloudViewer:public CViewer {

    typedef CPointCloud<std::list,float,3> C3dPointCloud;

    Q_OBJECT

public:

    //! Constructor.
    explicit CPointCloudViewer(const R4R::CView<float>& view, const C3dPointCloud* pcl=nullptr, QWidget* parent = nullptr);

    //! Updates point cloud and recompute its bounding box.
    void setPointCloud(const C3dPointCloud* pcl);

    /*! \brief Triggers update of the point cloud without changing the pointer to it.
     *
     * This comes in handy e.g. when new points are added to the cloud.
     */
    void updateBoundingBox();

public slots:


protected:

    //! \copydoc CViewer::paintGL()
    void paintGL();

    //! Update clip depth based on the bounding box approximation of the current point cloud.
    void updateClipDepth(const CView<float>& view, double tolerance = 1.0);

private:

    const C3dPointCloud* m_pcl;     //!< pointer to a point cloud


};

}

#endif

#endif // VIEWER_H
