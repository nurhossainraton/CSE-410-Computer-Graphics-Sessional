#include <windows.h>
#include <GL/glut.h>
#include <cmath>
#define PI (2*acos(0.0))
#include<stdio.h>
struct point
{
    double x, y, z;
};

struct point pos,l,r,u;

float triangleVertex;
float triangleMax;
float triangleMin;
float sphereMax;
float sphereMin;
float triangleCur;
float sphereCur;
float angle_rotation = 0.0;
float radius = 0.0;

int drawaxes=1;

void drawAxes()
{
    glLineWidth(3);
    if(drawaxes)
    {
        glBegin(GL_LINES);
        glColor3f(1,0,0);
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);

        glColor3f(0,0,1);

        glVertex3f(0,0,0);
        glVertex3f(0,1,0);

        glColor3f(0,1,0);
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
        glEnd();
    }
}
void cubemaking(struct point points[][101])
{
    glBegin(GL_QUADS);
    for (int j = 0; j < 100; j++)
    {
        for (int i = 0; i < 100; i++)
        {
            glVertex3f(points[j][i].x, points[j][i].y, points[j][i].z);
            glVertex3f(points[j][i+1].x, points[j][i+1].y, points[j][i+1].z);
            glVertex3f(points[j+1][i+1].x, points[j+1][i+1].y, points[j+1][i+1].z);
            glVertex3f(points[j+1][i].x, points[j+1][i].y, points[j+1][i].z);
        }
    }
    glEnd();
}
void drawCube()
{
    struct point points[101][101];
    for (int j = 0; j <= 100; j++)
    {
        for (int i = 0; i <= 100; i++)
        {
            points[j][i].x = -1.0+i*0.02;
            points[j][i].y = -1.0+j*0.02;
            points[j][i].z = 1.0;

            float len = points[j][i].x*points[j][i].x+points[j][i].y*points[j][i].y+points[j][i].z*points[j][i].z;
            len = sqrt(len);

            points[j][i].x = points[j][i].x/len;
            points[j][i].y = points[j][i].y/len;
            points[j][i].z = points[j][i].z/len;

            points[j][i].x =points[j][i].x * sphereCur;
            points[j][i].y =points[j][i].y * sphereCur;
            points[j][i].z =points[j][i].z * sphereCur;

        }
    }

    cubemaking(points);
}

void drawCubeSphere()
{

    for(int i=0; i<4; i++)
    {

        glPushMatrix();
        {
            glColor3f(0, i%2, (i+1)%2);
            glRotatef(90*i,0,1,0);
            glTranslatef(0,0,triangleCur);
            drawCube();

        }
        glPopMatrix();

    }

    for(int i=0; i<2; i++)
    {

        glPushMatrix();
        {
            glColor3f(1.0f, 0.0f, 0.0f);
            glRotatef(90+180*i,1,0,0);
            glTranslatef(0,0,triangleCur);
            drawCube();

        }
        glPopMatrix();

    }

}
void trianglehelper()
{
    glTranslatef((1.0-triangleCur)/3.0, (1.0-triangleCur)/3.0, (1.0-triangleCur)/3.0);
    glScaled(triangleCur, triangleCur, triangleCur);
    glBegin(GL_TRIANGLES);
    glVertex3f( 0.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
}
/* Draw an octahedron centered at the origin */
void drawOctahedron()
{
    // Front
    glPushMatrix();
    glColor3f(0.0f, 1.0f, 1.0f);
    trianglehelper();
    glPopMatrix();

    // Right
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 1.0f);
    glRotatef(90, 0, 0, 1);
    trianglehelper();
    glPopMatrix();

    // Back
    glPushMatrix();
    glColor3f(0.0f, 1.0f, 1.0f);
    glRotatef(180, 0, 0, 1);
    trianglehelper();
    glPopMatrix();

    // Left
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 1.0f);
    glRotatef(-90, 0, 0, 1);
    trianglehelper();
    glPopMatrix();

    // lower front
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 1.0f);
    glRotatef(180, 1, 1, 0);
    trianglehelper();
    glPopMatrix();

    // lower Right
    glPushMatrix();
    glColor3f(0.0f, 1.0f, 1.0f);
    glRotatef(90, 0, 0, 1);
    glRotatef(180, 1, 1, 0);
    trianglehelper();
    glPopMatrix();

    //lower back
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 1.0f);
    glRotatef(180, 0, 0, 1);
    glRotatef(180, 1, 1, 0);
    trianglehelper();
    glPopMatrix();

    //lower left
    glPushMatrix();
    glColor3f(0.0f, 1.0f, 1.0f);
    glRotatef(-90, 0, 0, 1);
    glRotatef(180, 1, 1, 0);
    trianglehelper();
    glPopMatrix();
}

void drawCylinder(double height, double radius, int segments)
{
    struct point points[segments+1];
    double phi = 70.5287794*PI /180; //70.5287794 degrees to radians
    for (int i = 0; i < segments+1; i++)
    {
        double theta = (i * phi / segments) - phi/2;
        points[i].x = radius * cos(theta);
        points[i].y = radius * sin(theta);
    }

    glBegin(GL_QUADS);
    for (int i = 0; i < segments; i++)
    {
        glVertex3f(points[i].x, points[i].y, height/2);
        glVertex3f(points[i].x, points[i].y, -height/2);
        glVertex3f(points[i+1].x, points[i+1].y, -height/2);
        glVertex3f(points[i+1].x, points[i+1].y, height/2);
    }
    glEnd();
}
void cylindercaller()
{
    glTranslatef(triangleCur/sqrt(2.0), 0.0, 0.0);
    drawCylinder(triangleCur*sqrt(2), sphereCur, 100);
}
void drawInternalCyl()
{
    //upper front
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //upper left
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(90, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //upper back
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //upper right
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(-90, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //middle front
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(90, 1, 0, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //middle left
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(-90, 1, 0, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //middle back
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(-90, 1, 0, 0);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //middle right
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(90, 1, 0, 0);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //lower front
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(180, 1, 0, 1);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //lower left
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(180, 1, 0, 1);
    glRotatef(90, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //lower back
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(180, 1, 0, 1);
    glRotatef(180, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

    //lower right
    glPushMatrix();
    glColor3f(1.0, 1.0, 1.0);
    glRotatef(180, 1, 0, 1);
    glRotatef(-90, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    glRotatef(45, 0, 1, 0);
    cylindercaller();
    glPopMatrix();

}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix


    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);
    // draw
    glRotatef(angle_rotation, 0, 1, 0);
    drawAxes();
    drawOctahedron();
    drawCubeSphere();

    drawInternalCyl();

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

    gluPerspective(45.0f, 1.1, 0.1f, 100.0f);

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

