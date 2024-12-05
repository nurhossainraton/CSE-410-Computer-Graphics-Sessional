#include<windows.h>
#include<GL/glut.h>
#include "bitmap_image.hpp"
#include "1905117_classes.h"



int drawaxes=1;
int pixelsAlongBothDimensions = 0;
int noOfObjects = 0;
int noOfPointLights = 0;
int noOfSpotLights = 0;
int imageCount=0;

Vector3D normalize(const Vector3D& v)
{
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

    if (len == 0) return Vector3D{0, 0, 0};
    return Vector3D{v.x / len, v.y / len, v.z / len};
}

Vector3D crossProduct(const Vector3D& a, const Vector3D& b)
{
    return Vector3D
    {
        a.y * b.z - b.y * a.z,  // x component
        a.z * b.x - b.z * a.x,  // y component
        a.x * b.y - b.x * a.y   // z component
    };
}




Vector3D rotateVector(const Vector3D& v, const Vector3D& reference, double rotationAngle)
{
    Vector3D vPerp = crossProduct(reference, v);
    double rotationAngleRadians = rotationAngle * PI / 180;
    double cosAngle = cos(rotationAngleRadians);
    double sinAngle = sin(rotationAngleRadians);


    Vector3D rotated
    {
        v.x * cosAngle + vPerp.x * sinAngle,
        v.y * cosAngle + vPerp.y * sinAngle,
        v.z * cosAngle + vPerp.z * sinAngle
    };

    return normalize(rotated);
}

void drawAxes()
{
    if(drawaxes == 1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);

            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);

            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }
        glEnd();
    }
}

void capture() {
    cout << "Bitmap Image Capturing" << endl;
    int imageWidth = pixelsAlongBothDimensions,imageHeight = pixelsAlongBothDimensions;

    //initialize bitmap image and set background color
    bitmap_image image(imageWidth, imageHeight);

    for(int column=0; column<imageWidth; column++) {
        for(int row=0; row<imageHeight; row++) {
            image.set_pixel(column, row, 0, 0, 0);  // color = black
        }
    }

    double planeDistance = 500/(2.0*tan(80 / 2.0 * PI / 180.0));
    Vector3D topLeft = pos + l * planeDistance - r * (500 / 2.0) + u * (500 / 2.0);

    double du = ((double) 500/imageWidth);
    double dv = ((double) 500/imageHeight);

    // Choose middle of the grid cell
    topLeft = topLeft + r * (du / 2.0) - u * (dv / 2.0);

    for(int column=0; column<imageWidth; column++) {
        for(int row=0; row<imageHeight; row++) {
            //calculate curPixel using topleft,r,u,i,j,du,dv
            Vector3D curPixel = topLeft + r * (column * du) - u * (row * dv);

            //cast ray from eye to (curPixel-eye) direction
            Ray* ray = new Ray(pos, curPixel-pos);

            int nearest = INF;
            double t, tMin=INF;

            for(int i=0; i<objects.size(); i++) {
                Color color;  // color = black
                t = objects[i]->intersect(ray, color, 0);

                if(t>0.0 && t<tMin) {
                    tMin = t;
                    nearest = i;
                }
            }

            /* finding color for current pixel */
            if(nearest != INF) {
                Color color;  // color = black
                tMin = objects[nearest]->intersect(ray, color, 1);
                int red = round(color.red * 255.0);
                int green = round(color.green * 255.0);
                int blue = round(color.blue * 255.0);
                image.set_pixel(column, row, red, green, blue);
            }
        }
    }
    imageCount++;
    string imageName = "Output_1" + to_string(imageCount) + ".bmp";
    image.save_image("../" + imageName);
    cout << "Bitmap Image Captured | Image Name - " + imageName << endl;
}

void moveForwardBackward(int c)
{
    pos.x = pos.x + (c*l.x);
    pos.y = pos.y + (c*l.y);
    pos.z = pos.z + (c*l.z);
}

void moveLeftRight(int c)
{
    pos.x += c * r.x;
    pos.y += c * r.y;
    pos.z += c * r.z;
}

void moveUpDown(int c)
{
    pos.x += c * u.x;
    pos.y += c * u.y;
    pos.z += c * u.z;
}

void lookLeftRight(int angle)
{
    l = rotateVector(l, u, angle);
    r = rotateVector(r, u, angle);
}

void lookUpDown(int angle)
{
    u = rotateVector(u, r, angle);
    l = rotateVector(l, r, angle);
}

void tilt(int angle)
{
    u = rotateVector(u, l, angle);
    r = rotateVector(r, l, angle);
}


void keyboardListener(unsigned char key, int x, int y)
{
    if (key == '0')
    {
        capture();
    }
    else if (key == '1')
    {
        lookLeftRight(-1);
    }
    else if (key == '2')
    {
        lookLeftRight(1);
    }
    else if (key == '3')
    {
        lookUpDown(1);
    }
    else if (key == '4')
    {
        lookUpDown(-1);
    }
    else if (key == '5')
    {
        tilt(-2);
    }
    else if (key == '6')
    {
        tilt(2);
    }
}

void specialKeyListener(int key, int x, int y)
{
    if (key == GLUT_KEY_DOWN)            // Down arrow key
    {
        moveForwardBackward(-1);
    }
    else if (key == GLUT_KEY_UP)         // Up arrow key
    {
        moveForwardBackward(1);
    }
    else if (key == GLUT_KEY_RIGHT)      // Right arrow key
    {
        moveLeftRight(1);
    }
    else if (key == GLUT_KEY_LEFT)       // Left arrow key
    {
        moveLeftRight(-1);  // Corrected from '=1' to '-1'
    }
    else if (key == GLUT_KEY_PAGE_UP)     // Page Up key
    {
        moveUpDown(1);
    }
    else if (key == GLUT_KEY_PAGE_DOWN)     // Page Down key
    {
        moveUpDown(-1);
    }
    // No default action needed as all cases are covered
}



void drawObjects()
{
    for (auto& object : objects)
    {
        object->draw();
    }
}


void drawLights()
{
    for (auto& light : lights)
    {
        light.draw();
    }

    for (auto& spotlight : spotlights)
    {
        spotlight.draw();
    }

}

void display()
{
    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    gluLookAt(pos.x,pos.y,pos.z,	pos.x + l.x,pos.y + l.y,pos.z + l.z,	u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    /* adding axes */
    drawAxes();

    /* adding objects */
    drawObjects();

    /* adding lights */
    drawLights();

    /* ADD this line in the end: if you use double buffer (i.e. GL_DOUBLE) */
    glutSwapBuffers();
}

void animate()
{
    /* codes for any changes in Models, Camera */
    glutPostRedisplay();
}

void init()
{

    pos.x = 80,
    pos.y = 80,
    pos.z = 80;

    u.x = 0,
    u.y = 0,
    u.z = 1;

    r.x = -1/ sqrt(2),
    r.y = 1/ sqrt(2),
    r.z = 0;

    l.x = -1/ sqrt(2),
    l.y = -1/ sqrt(2),
    l.z = 0;

    glClearColor(0, 0, 0, 0);  // color = black

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    gluPerspective(80.0, 1.0, 1.0, 1000.0);
}

void loadData()
{
    ifstream fin;
    string type;
    Object* temp;
    fin.open("scene.txt");

    if(!fin.is_open())
    {
        cout << "Error Opening File" << endl;
    }
    fin >> levelOfRecursion >> pixelsAlongBothDimensions;
    fin >> noOfObjects;
    for(int i=0; i < noOfObjects; i++)
    {
        fin >> type;
        if(type == "sphere")
        {
            Vector3D center;
            double radius;

            fin >> center.x >> center.y >> center.z;
            fin >> radius;

            temp = new Sphere(center, radius);
        }
        else if(type == "triangle")
        {
            Vector3D a, b, c;

            fin >> a.x >> a.y >> a.z;
            fin >> b.x >> b.y >> b.z;
            fin >> c.x >> c.y >> c.z;

            temp = new Triangle(a, b, c);
        }
        else if(type == "general")
        {
            GeneralQuadricSurfaceCoefficient coefficient;
            Vector3D cubeReferencePoint;
            double length, width, height;

            fin >> coefficient.a >> coefficient.b >> coefficient.c >> coefficient.d >> coefficient.e >> coefficient.f >> coefficient.g >> coefficient.h >> coefficient.i >> coefficient.j;
            fin >> cubeReferencePoint.x >> cubeReferencePoint.y >> cubeReferencePoint.z;
            fin >> length >> width >> height;

            temp = new GeneralQuadricSurface(coefficient, cubeReferencePoint, length, width, height);
        }
        else
        {
            cout << type << " Shape Type error" << endl;
        }

        Color color;
        ReflectionCoefficient reflectionCoefficient;
        int shine;

        fin >> color.red >> color.green >> color.blue;
        fin >> reflectionCoefficient.ambient >> reflectionCoefficient.diffuse >> reflectionCoefficient.specular >> reflectionCoefficient.recursive;
        fin >> shine;

        temp->setColor(color);
        temp->setReflectionCoefficient(reflectionCoefficient);
        temp->setShine(shine);

        objects.push_back(temp);
    }
    temp = NULL;

    fin >> noOfPointLights;
    for(int i=0; i < noOfPointLights; i++)
    {
        Vector3D position;
        Color color;

        fin >> position.x >> position.y >> position.z;
        fin >> color.red >> color.green >> color.blue;

        lights.push_back(PointLight(position, color, 1.0));
    }
    fin >> noOfSpotLights;
    for(int i=0; i < noOfSpotLights; i++)
    {
        Vector3D position,direction;
        Color color;
        double cutoffAngle;

        fin >> position.x >> position.y >> position.z;
        fin >> color.red >> color.green >> color.blue;
        fin >> direction.x >> direction.y >> direction.z;
        fin >> cutoffAngle;

        spotlights.push_back(SpotLight(cutoffAngle, 2.0,position, color, direction));
    }
    fin.close();

    temp = new Floor(1000.0, 20.0);

    temp->setColor(Color(1.0, 1.0, 1.0));
    temp->setReflectionCoefficient(ReflectionCoefficient(0.25, 0.25, 0.25, 0.25));
    temp->setShine(15);

    objects.push_back(temp);
    temp = NULL;
}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);


    loadData();

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
