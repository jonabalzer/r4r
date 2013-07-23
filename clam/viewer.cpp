#include "viewer.h"

using namespace R4R;
using namespace std;

CViewer::CViewer(QWidget *parent):
    QGLWidget(parent),
    m_cam(),
    m_view(m_cam),
    m_last_point() {}


CViewer::CViewer(const CView<float>& view, QWidget *parent):
    QGLWidget(parent),
    m_cam(),
    m_view(view),
    m_last_point() {}

void CViewer::initializeGL(){

    glClearColor(0, 0, 0, 1.0);
    glClearDepth(1);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glDepthFunc(GL_LEQUAL);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    GLfloat pos1[] = { 0.1,  0.1, -0.02, 0.0};
    GLfloat pos2[] = {-0.1,  0.1, -0.02, 0.0};
    GLfloat pos3[] = { 0.0,  0.0,  0.1,  0.0};
    GLfloat col1[] = { 0.7,  0.7,  0.8,  1.0};
    GLfloat col2[] = { 0.8,  0.7,  0.7,  1.0};
    GLfloat col3[] = { 0.5,  0.5,  0.5,  1.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_POSITION, pos1);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
    glLightfv(GL_LIGHT0,GL_SPECULAR, col1);

    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT1,GL_POSITION, pos2);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
    glLightfv(GL_LIGHT1,GL_SPECULAR, col2);

    glEnable(GL_LIGHT2);
    glLightfv(GL_LIGHT2,GL_POSITION, pos3);
    glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
    glLightfv(GL_LIGHT2,GL_SPECULAR, col3);

    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);

    load_view(m_view);
    load_camera(m_cam);

}

void CViewer::resizeGL(int w, int h){

    glViewport(0,0,w,h);
    updateGL();

}

void CViewer::load_camera(const R4R::CPinholeCam& cam) {

    m_cam = cam;

    CVector<size_t,2> sizes = cam.GetSize();

    cout << "Sizes: " << sizes.Get(0) << " " << sizes.Get(1) << endl;
    double proj[16];

    mat K = cam.GetProjectionMatrix();

    proj[0] = 1.0; //2*K.Data().get()[0] / sizes.Get(0);
    proj[5] = 1.3; //2*K.Data().get()[4] / sizes.Get(1);

    double znear = 0.001;
    double zfar = 1000;

    proj[10] = (-zfar - znear)/(zfar - znear);
    proj[14] = -2*zfar*znear/(zfar - znear);

    proj[11] = -1.0;
    proj[15] = 0.0;

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(&proj[0]);

}


void CViewer::paintGL(){

    // set projection and modelview matrices from m_view
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    list<pair<vec3f,rgb> >::iterator it;
    glPointSize(2);

    glBegin(GL_POINTS);

    for(it=m_map.begin(); it!=m_map.end(); it++) {

        vec3f x = it->first;
        rgb color = it->second;

        glVertex3f(x.Get(0),x.Get(1),x.Get(2));
        glColor3ub(color.Get(0),color.Get(1),color.Get(2));

    }

    glEnd();

}

void CViewer::load_view(const R4R::CView<float>& view) {

    const CRigidMotion<float,3> F = view.GetTransformation();
    matf mF = matf(F);

    float mv[16];

    mv[0] = mF.Data().get()[0];
    mv[1] = -mF.Data().get()[1];
    mv[2] = -mF.Data().get()[2];
    mv[3] = mF.Data().get()[3];
    mv[4] = mF.Data().get()[4];
    mv[5] = -mF.Data().get()[5];
    mv[6] = -mF.Data().get()[6];
    mv[7] = mF.Data().get()[7];
    mv[8] = mF.Data().get()[8];
    mv[9] = -mF.Data().get()[9];
    mv[10] = -mF.Data().get()[10];
    mv[11] = mF.Data().get()[11];
    mv[12] = mF.Data().get()[12];
    mv[13] = -mF.Data().get()[13];
    mv[14] = -mF.Data().get()[14];
    mv[15] = mF.Data().get()[15];

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(&mv[0]);

}

void CViewer::update_display(const R4R::CView<float>& view) {

    m_view = view;

    load_view(m_view);

    updateGL();

}

void CViewer::translate(float dx, float dy, float dz) {

    // new cam coordinates in old
    CRigidMotion<float,3> F(dx,dy,dz,0,0,0);
    F.Invert();

    const CRigidMotion<float,3>& Fa = m_view.GetTransformation();
    CTransformation<float,3> result = F*Fa;

    m_view.SetTransformation(reinterpret_cast<CRigidMotion<float,3>& >(result));

    update_display(m_view);


}

void CViewer::wheelEvent(QWheelEvent *event) {

    float d = -(float)event->delta() / 120.0 * 0.2 * 10;
    translate(0,0,d);
    event->accept();

}


void CViewer::mousePressEvent(QMouseEvent* event) {

    if(event->button()==Qt::RightButton)
        m_last_point = event->pos();

    event->accept();

}

void CViewer::mouseMoveEvent(QMouseEvent* event) {

    QPoint pt = event->pos();

    float dx = pt.x() - m_last_point.x();
    float dy = pt.y() - m_last_point.y();

    translate(-dx*0.05,-dy*0.05,0);

    m_last_point = pt;

    event->accept();

}


void CViewer::mouseReleaseEvent(QMouseEvent *event) {

    event->accept();

}

