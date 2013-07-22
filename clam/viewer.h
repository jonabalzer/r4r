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
