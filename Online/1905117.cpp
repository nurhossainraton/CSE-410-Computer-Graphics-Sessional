#include <windows.h>
#include <GL/glut.h>
#include <cmath>
#define PI (2*acos(0.0))
#include<stdio.h>
using namespace std;
struct point
{
    double x, y, z;
};
double radius=1.0;
struct point pos,l,r,u;

float triangleVertex;
float triangleMax;
float triangleMin;
float sphereMax;
float sphereMin;
float triangleCur;
float sphereCur;
float angle_rotation = 0.0;


int drawaxes=1;

void drawAxes()
{
    glLineWidth(1);
    if(drawaxes)
    {
        glBegin(GL_LINES);

        glVertex3f(-1,0,0);
        glVertex3f(1,0,0);

        glVertex3f(0,-1,0);
        glVertex3f(0,1,0);

        glVertex3f(0,0,-1);
        glVertex3f(0,0,1);
        glEnd();
    }
}

void drawCircle()
{

   glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 360; i++) {
        float theta = i * 3.14159265358979323846 / 180.0;
        float x = radius * cos(theta);
        float y = radius * sin(theta);
        glVertex3f(x, y, 0.0f);
    }
    glEnd();

}

void drawLine()
{
    glBegin(GL_LINES);

        glVertex3f(1,-1,0);
        glVertex3f(0,0,1);

        glVertex3f(1.5,-1,0);
        glVertex3f(0,0,1);

        glVertex3f(2,-1,0);
        glVertex3f(0,0,1);

        glVertex3f(-1.5,-1,0);
        glVertex3f(0,0,1);

        glVertex3f(-2,-1,0);
        glVertex3f(0,0,1);

        glVertex3f(-1,-1,0);
        glVertex3f(0,0,1);
        glEnd();
}

void drawrect()
{
    glBegin(GL_QUADS);
    glVertex3f(-2,-0.5,0);
    glVertex3f(2.5,-1,0);
    glVertex3f(-2,0.5,0);
    glVertex3f(2.5,-1,0);
     glEnd();


}
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix


    gluLookAt(pos.x,pos.y,pos.z,
              0,0,0,
              u.x,u.y,u.z);
    // draw
    glRotatef(angle_rotation, 0, 1, 0);
    drawAxes();
    drawCircle();
    drawLine();
    drawrect();


    glutSwapBuffers();  // Render now
}
void lookleftandright(double rate,int c1)
{
    r.x = r.x*cos(rate)+l.x*sin(c1*rate);
    r.y = r.y*cos(rate)+l.y*sin(c1*rate);
    r.z = r.z*cos(rate)+l.z*sin(c1*rate);

    l.x = l.x*cos(rate)-r.x*sin(c1*rate);
    l.y = l.y*cos(rate)-r.y*sin(c1*rate);
    l.z = l.z*cos(rate)-r.z*sin(c1*rate);
}

void lookupanddown(double rate,int c)
{
    u.x = u.x*cos(rate)-l.x*sin(c*rate);
    u.y = u.y*cos(rate)-l.y*sin(c*rate);
    u.z = u.z*cos(rate)-l.z*sin(c*rate);

    l.x = l.x*cos(rate)+u.x*sin(c*rate);
    l.y = l.y*cos(rate)+u.y*sin(c*rate);
    l.z = l.z*cos(rate)+u.z*sin(c*rate);
}

void tiltclock(double rate,int c)
{
    r.x = r.x*cos(rate)-u.x*sin(c*rate);
    r.y = r.y*cos(rate)-u.y*sin(c*rate);
    r.z = r.z*cos(rate)-u.z*sin(c*rate);

    u.x = u.x*cos(rate)+r.x*sin(c*rate);
    u.y = u.y*cos(rate)+r.y*sin(c*rate);
    u.z = u.z*cos(rate)+r.z*sin(c*rate);
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y)
{
    double rate = 0.01,v = 0.25,m = 0.1,s;

    // Control eye (location of the eye)
    // control eyex
    if(key=='1')
    {
        lookleftandright(rate,1);
    }


    else if(key=='2')
    {
        lookleftandright(rate,-1);

    }


    else if(key=='3')
    {
        lookupanddown(rate,1);
    }
    else if(key=='4')
    {
        lookupanddown(rate,-1);

    }

    else if(key=='5')
    {
        tiltclock(rate,1);

    }
    else if(key=='6')
    {
        tiltclock(rate,-1);

    }
    else if(key=='a')
    {
        angle_rotation =angle_rotation+2.0f;

    }
    else if(key=='d')
        angle_rotation =angle_rotation-2.0f;

    else if(key=='w')
        pos.z =pos.z+m;

    else if(key=='s')
        pos.z =pos.z-m;

    else if(key==',')
    {
        if(triangleCur >= 0.0625) triangleCur=triangleCur-triangleMax/16.0;
        if(sphereCur <= sphereMax-0.036084)sphereCur =sphereCur+sphereMax/(16.0);
        drawaxes = 0;
    }
    else if(key=='.')
    {
        if(triangleCur <= triangleMax-0.0625) triangleCur=triangleCur+triangleMax/16.0;
        if(sphereCur >= 0.036084)sphereCur =sphereCur-sphereMax/16.0;
    }

    else if(key==27)
    {
        exit(0);    // Exit window
    }

    glutPostRedisplay();
}

void shiftRightLeft(int c)
{
    pos.x =pos.x + c* r.x;
    pos.y =pos.y + c* r.y;
    pos.z =pos.z + c* r.z;
}
void shiftpageUpDown(int c)
{
    pos.x =pos.x + c*u.x;
    pos.y =pos.y + c*u.y;
    pos.z =pos.z + c*u.z;
}
void shiftUpDown(int c)
{
    pos.x =pos.x + c*l.x;
    pos.y =pos.y + c*l.y;
    pos.z =pos.z + c*l.z;
}
void specialKeyListener(int key, int x, int y)
{
    if (key == GLUT_KEY_RIGHT)
    {
        shiftRightLeft(1);
    }
    else if (key == GLUT_KEY_LEFT)
    {
        shiftRightLeft(-1);
    }
    else if (key == GLUT_KEY_PAGE_UP)
    {
        shiftpageUpDown(1);

    }
    else if (key == GLUT_KEY_PAGE_DOWN)
    {
        shiftpageUpDown(-1);
    }
    else if (key == GLUT_KEY_UP)
    {
        shiftUpDown(1);
    }
    else if (key == GLUT_KEY_DOWN)
    {
        shiftUpDown(-1);
    }

    glutPostRedisplay();
}


void init()
{
    sphereCur = 0.0;
    sphereMax = 1.0/sqrt(3);
    sphereMin = 0.0;

    pos.x=0;
    pos.y=0;
    pos.z=5;

    l.x=0;
    l.y=0;
    l.z=-1;
    u.x=0;
    u.y=1;
    u.z=0;
    r.x=1;
    r.y=0;
    r.z=0;

    triangleMin = 0.0;
    triangleMax = 1.0;
    triangleCur = 1.0;

    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix

    gluPerspective(45.0f, 1, 0.1f, 100.0f);

}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(50, 50);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("OpenGL 3D Drawing");
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();
    return 0;
}

