#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include <cmath>
#define PI (2*acos(0.0))
#include<stdio.h>
using namespace std;

#include <GL/freeglut.h>


struct point
{
    double x, y, z;
};
GLfloat angleX = 0.0;
GLfloat angleY = 0.0;

const int numDivisions = 30;
double directionx = 5.0,directiony = 7.0,positionx=5.0,positiony=7.0;
struct point pos,l,r,u,direction;
float radius = 10.0;
float angle_rotation = 0.0;
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
    double rate = 0.05,v = 0.25,m = 0.2;

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
     else if(key=='i')
     {
        positionx+=2;
        if(positionx >= 50 && positionx<= -50)
        {
            positionx = 0;
        }
     }

     else if(key=='k')
     {
         positionx-=2;
         if(positionx == 50 )
            positionx =0;
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

void drawSphere() {
    glPushMatrix();
    glTranslatef(positionx,positiony,radius);
    for (int i = 0; i < numDivisions; ++i) {
        double theta1 = (i * 2 * PI) / numDivisions;
        double theta2 = ((i + 1) * 2 * PI) / numDivisions;

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= numDivisions; ++j) {
            if ((i + j) % 2 == 0)
            {
                glColor3f(1.0, 0.0, 0.0);  // Black color
            }
            else
            {
                glColor3f(0.0, 1.0, 0.0);  // White color
            }
            double phi = (j * PI) / numDivisions;

            double x1 = radius * sin(theta1) * cos(phi);
            double y1 = radius * sin(theta1) * sin(phi);
            double z1 = radius * cos(theta1);

            double x2 = radius * sin(theta2) * cos(phi);
            double y2 = radius * sin(theta2) * sin(phi);
            double z2 = radius * cos(theta2);

            glVertex3d(x1, y1, z1);
            glVertex3d(x2, y2, z2);
        }
        glEnd();
    }
    glPopMatrix();
}
void drawrectangle()
{
    for(int i=0; i<4; i++)
    {
        glRotatef(i*90,0,0,1);
        glColor3f(1,0,0);
        glBegin(GL_QUADS);
        {

            glVertex3f(-50,-50,0);
            glVertex3f(50,-50,0);
            glVertex3f(50,-50,10);
            glVertex3f(-50,-50,10);

        }
        glEnd();

    }


}
void drawCheckerboard()
{
    int numRows = 20;
    int numCols = 20;
    float squareSize = 20.0;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            // Alternate between black and white squares
            if ((i + j) % 2 == 0)
            {
                glColor3f(0.0, 0.0, 0.0);  // Black color
            }
            else
            {
                glColor3f(1.0, 1.0, 1.0);  // White color
            }

            // Draw a square
            glBegin(GL_QUADS);
            glVertex3f(j * squareSize - 200.0, i * squareSize - 200.0, 0);
            glVertex3f((j + 1) * squareSize - 200.0, i * squareSize - 200.0, 0);
            glVertex3f((j + 1) * squareSize - 200.0, (i + 1) * squareSize - 200.0, 0);
            glVertex3f(j * squareSize - 200.0, (i + 1) * squareSize - 200.0, 0);
            glEnd();
        }
    }

    drawrectangle();
    drawSphere();

    glutSwapBuffers();
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (GLfloat)width / (GLfloat)height, 1.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}
void updateCamera()
{
    glLoadIdentity();
    gluLookAt(pos.x,pos.y,pos.z,
              l.x,l.y,l.z,
              u.x,u.y,u.z);
}
void init()
{

    pos.x=-200;
    pos.y=15;
    pos.z=50;

    l.x=1;
    l.y=1;
    l.z=1;

    u.x=0;
    u.y=0;
    u.z=1;

    r.x=1;
    r.y=0;
    r.z=0;


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, 1.1, 0.1f, 100.0f);

}
int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutCreateWindow("3D Checkerboard");
    glutFullScreen();
    glEnable(GL_DEPTH_TEST);

    init();
    glutDisplayFunc(drawCheckerboard);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutIdleFunc(updateCamera);

    glutMainLoop();
    return 0;
}


