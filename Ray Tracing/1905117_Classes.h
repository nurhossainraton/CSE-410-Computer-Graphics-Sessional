#include<bits/stdc++.h>

using namespace std;

#define PI 2*acos(0.0)
#define INF 1e9


class Vector3D
{
public:
    double x, y, z;


    Vector3D();


    Vector3D(double x, double y, double z);


    void normalize();


    Vector3D operator+(const Vector3D& v);
    Vector3D operator-(const Vector3D& v);
    Vector3D operator*(const double scalar);
};
Vector3D::Vector3D() : x(0.0), y(0.0), z(0.0) {}
Vector3D::Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

void Vector3D::normalize()
{
    double len = sqrt(x*x + y*y + z*z);
    if (len > 0)
    {
        x /= len;
        y /= len;
        z /= len;
    }
}


Vector3D Vector3D::operator+(const Vector3D& v)
{
    return Vector3D(x + v.x, y + v.y, z + v.z);
}
Vector3D Vector3D::operator-(const Vector3D& v)
{
    return Vector3D(x - v.x, y - v.y, z - v.z);
}

Vector3D Vector3D::operator*(const double scalar)
{
    return Vector3D(x * scalar, y * scalar, z * scalar);
}

double DOT(Vector3D v1, Vector3D v2)
{
    double x=v1.x*v2.x;
    double y=v1.y*v2.y;
    double z=v1.z*v2.z;
    return  x+y+z ;
}
Vector3D CROSS(Vector3D v1, Vector3D v2)
{
    return Vector3D(
               v1.y * v2.z - v1.z * v2.y,
               v1.z * v2.x - v1.x * v2.z,
               v1.x * v2.y - v1.y * v2.x
           );
}

double ValueOfVector(Vector3D v)
{
    double squaredMagnitude = pow(v.x, 2.0) + pow(v.y, 2.0) + pow(v.z, 2.0);
    double magnitude = sqrt(squaredMagnitude);
    return magnitude;
}

double distanceBetweenPoints(Vector3D p1, Vector3D p2)
{

    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;

    double squaredDistance = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);

    double distance = sqrt(squaredDistance);

    return distance;
}
class Ray
{
public:
    Vector3D R0;
    Vector3D Rd;

    Ray(Vector3D R0, Vector3D Rd); // Constructor declaration
};

Ray::Ray(Vector3D R0, Vector3D Rd) : R0(R0), Rd(Rd)
{
    this->Rd.normalize();
}

class Color
{
public:
    double red, green, blue;

    Color();

    Color(double red, double green, double blue);

    void clipColor();
};

Color::Color() : red(0), green(0), blue(0) {}

Color::Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}

void Color::clipColor()
{
    if (red > 1.0) red = 1.0;
    if (green > 1.0) green = 1.0;
    if (blue > 1.0) blue = 1.0;
}

class PointLight
{
public:
    Vector3D position;
    Color color;
    double radius;

    PointLight();
    PointLight(Vector3D position, Color color, double radius);
    Vector3D getPosition();
    Color getColor();
    void draw();
};

PointLight::PointLight()
{
    radius = 0.0;
}

PointLight::PointLight(Vector3D position, Color color, double radius) : position(position), color(color), radius(radius) {}

Vector3D PointLight::getPosition()
{
    return position;
}


Color PointLight::getColor()
{
    return color;
}


void PointLight::draw()
{
    glColor3f(color.red,color.green,color.blue);
    glPointSize(10);
    glBegin(GL_POINTS);
    {
        glVertex3f(position.x,position.y,position.z);
    }
    glEnd();

}

class SpotLight
{
public:
    Vector3D position;
    Color color;
    Vector3D direction;
    double cutoffAngle;
    double radius;

    SpotLight();
    SpotLight(double cutoffAngle, double radius, Vector3D position, Color color, Vector3D direction);
    Vector3D getPosition();
    Color getColor();
    void draw();
};

// Constructor definition for SpotLight
SpotLight::SpotLight()
{
    radius = 0.0;
}

SpotLight::SpotLight(double cutoffAngle, double radius, Vector3D position, Color color, Vector3D direction)
    : position(position), color(color), direction(direction), cutoffAngle(cutoffAngle), radius(radius) {}

// Function to get the position of the spot light
Vector3D SpotLight::getPosition()
{
    return position;
}

// Function to get the color of the spot light
Color SpotLight::getColor()
{
    return color;
}

// Function to draw the spot light
void SpotLight::draw()
{
    glColor3f(color.red, color.green, color.blue);
    glPointSize(10);
    glBegin(GL_POINTS);
    {
        glVertex3f(position.x,position.y,position.z);
    }
    glEnd();
}


class ReflectionCoefficient
{
public:
    double ambient, diffuse, specular, recursive;

    ReflectionCoefficient();
    ReflectionCoefficient(double ambient, double diffuse, double specular, double recursive);
};
ReflectionCoefficient::ReflectionCoefficient() : ambient(0.0), diffuse(0.0), specular(0.0), recursive(0.0) {}

ReflectionCoefficient::ReflectionCoefficient(double ambient, double diffuse, double specular, double recursive)
    : ambient(ambient), diffuse(diffuse), specular(specular), recursive(recursive) {}


class Object
{
public:
    Vector3D reference_point;
    Color color;
    double length;
    ReflectionCoefficient reflectionCoefficient;
    int shine;

    Object();

    void setColor(Color color);
    void setReflectionCoefficient(ReflectionCoefficient reflectionCoefficient);
    void setShine(int shine);
    Color getColor();
    virtual double intersect(Ray* ray, Color& color, int depth) = 0;
    virtual void draw() = 0;
};

Object::Object() : shine(0) {}

void Object::setColor(Color color)
{
    this->color = color;
}

void Object::setReflectionCoefficient(ReflectionCoefficient reflectionCoefficient)
{
    this->reflectionCoefficient = reflectionCoefficient;
}

void Object::setShine(int shine)
{
    this->shine = shine;
}
Color Object::getColor()
{
    return color;
}

int levelOfRecursion = 0;
vector<Object*> objects;
vector<PointLight> lights;
vector<SpotLight> spotlights;
Vector3D pos;
Vector3D u;
Vector3D r;
Vector3D l;


void pointLightDiffuseAndSpecular(int shine,Ray* ray, Vector3D& intersectionPoint, Vector3D& normal, Color& intersectionPointColor, Color& color, ReflectionCoefficient reflectionCoefficient)
{
    for(int i=0; i<lights.size(); i++)
    {

        Ray* incidentRay = new Ray(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        double t, tMinimum = INF;
        for(int j=0; j<objects.size(); j++)
        {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum)
            {
                tMinimum = t;
            }
        }
        double epsilon = 0.0000001;  // for tuning light effect
        Vector3D shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd * tMinimum;


        // if intersectionPoint is in shadow, the diffuse
        // and specular components need not be calculated
        if( distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)<distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7))
        {
            continue;
        }

        //calculate lambertValue using normal, rayl
        //find reflected ray, rayr for rayl
        //calculate phongValue using r, rayr
        double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
        lambertValue = max(lambertValue, 0.0);


        Vector3D reflectedDirection = incidentRay->Rd - normal * (DOT(incidentRay->Rd, normal) * 2.0);
        Ray* reflectedRay = new Ray(intersectionPoint, reflectedDirection);

        double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);
        phongValue = max(phongValue, 0.0);

        double diffuseContribution = reflectionCoefficient.diffuse * lambertValue;

        color.red = color.red + lights[i].getColor().red * intersectionPointColor.red * diffuseContribution;
        color.green = color.green + lights[i].getColor().green * intersectionPointColor.green * diffuseContribution;
        color.blue = color.blue + lights[i].getColor().blue * intersectionPointColor.blue * diffuseContribution;

        double specularContribution = reflectionCoefficient.specular * pow(phongValue, shine);
        color.red = color.red + lights[i].getColor().red * intersectionPointColor.red * specularContribution;
        color.green =color.green + lights[i].getColor().green * intersectionPointColor.green * specularContribution;
        color.blue =color.blue + lights[i].getColor().blue * intersectionPointColor.blue * specularContribution;

        color.clipColor();
    }
}

void spotLightDiffuseAndSpecular(int shine,Ray* ray, Vector3D& intersectionPoint, Vector3D& normal, Color& intersectionPointColor, Color& color, ReflectionCoefficient reflectionCoefficient)
{
    for (int i = 0; i < spotlights.size(); i++)
    {
        Vector3D direction = intersectionPoint - spotlights[i].getPosition();
        Ray* incidentRay = new Ray(spotlights[i].getPosition(), direction);
        Vector3D A,B;
        A = intersectionPoint - spotlights[i].getPosition();
        B = spotlights[i].direction;
        double angle = acos(DOT(A, B)/(ValueOfVector(A)*ValueOfVector(B))) * 180 / PI;

        if (spotlights[i].cutoffAngle < angle) continue;
        /* checking if intersection point is in shadow */
        double t, tMinimum = INF;
        for (int j = 0; j < objects.size(); j++)
        {
            Color dummyColor;  // light_color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);
            if (t > 0 && t < tMinimum)
            {
                tMinimum = t;
            }
        }
        // Scale the direction vector by tMinimum
        Vector3D scaledDirection = incidentRay->Rd * tMinimum;

        Vector3D shadowIntersectionPoint = incidentRay->R0 + scaledDirection;

        if( distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0) < distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7))
        {
            continue;
        }
        //calculate lambertValue using normal, rayl
        //find reflected ray, rayr for rayl
        //calculate phongValue using r, rayr
        // Calculate the reflection direction vector
        Vector3D incidentDirection = incidentRay->Rd;
        Vector3D scaledNormal = normal * (DOT(incidentDirection, normal) * 2.0);
        Vector3D reflectionDirection = incidentDirection - scaledNormal;


        Ray* reflectedRay = new Ray(intersectionPoint, reflectionDirection);
        Vector3D viewDirection = ray->Rd * (-1.0); // Reverse the view direction
        double phongValue = DOT(viewDirection, reflectedRay->Rd);

// Ensure phongValue is non-negative
        phongValue = max(phongValue, 0.0);

// Calculate the negative incident direction
        Vector3D negativeIncidentDirection = incidentRay->Rd * (-1.0);

// Compute the dot product for the Lambertian reflection (cosine of the angle)
        double lambertValue = DOT(negativeIncidentDirection, normal);

// Ensure the Lambertian value is non-negative
        lambertValue = max(lambertValue, 0.0);


        float colorContributionFactorRed = spotlights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * lambertValue;
        float colorContributionFactorGreen = spotlights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.diffuse * lambertValue;
        float colorContributionFactorBlue = spotlights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.diffuse * lambertValue;

// Explicitly updating each color component with the calculated color contribution factor
        color.red = color.red + colorContributionFactorRed;
        color.green = color.green + colorContributionFactorGreen;
        color.blue = color.blue + colorContributionFactorBlue;


        // Calculate the common multiplier for the specular reflection
        double specularMultiplier = reflectionCoefficient.specular * pow(phongValue, shine);

// Explicitly update each color component without using shorthand operators
        color.red = color.red + (spotlights[i].getColor().red * intersectionPointColor.red * specularMultiplier);
        color.green = color.green + (spotlights[i].getColor().green * intersectionPointColor.green * specularMultiplier);
        color.blue = color.blue + (spotlights[i].getColor().blue * intersectionPointColor.blue * specularMultiplier);

        color.clipColor();
    }
}

int findNearestObject(Ray* reflectedRay, double& tMinimum)
{
    int nearest = INF;
    double t;

    for(int i=0; i<objects.size(); i++)
    {
        Color dummyColor;
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 and t<tMinimum)
        {
            tMinimum = t;
            nearest = i;
        }
    }
    return nearest;
}


class Sphere: public Object
{
public:
    Vector3D center;
    double radius;

    Sphere(Vector3D center, double radius); // Constructor declaration
    void draw(); // draw method declaration
    double intersect(Ray* ray, Color& color, int level); // intersect method declaration
};

Sphere::Sphere(Vector3D center, double radius)
    : center(center), radius(radius)
{
    // Constructor body is empty since member initialization is done above
}



void Sphere::draw()
{
    glColor3f(getColor().red, getColor().green, getColor().blue);
    glPushMatrix();
    glTranslatef(center.x, center.y, center.z);
    glutSolidSphere(radius, 100, 100);
    glPopMatrix();
}

double Sphere::intersect(Ray* ray, Color& color, int level)
{
    //find intersect value t minimum of at^2+bt+c=0
    double a, b, c, tMin;

    a = DOT(ray->Rd, ray->Rd);

    double rayOriginDotRayDirection = DOT(ray->R0, ray->Rd);
    double rayDirectionDotCenter = DOT(ray->Rd, center);
    double subtractionResult = rayOriginDotRayDirection - rayDirectionDotCenter;

    b = subtractionResult * 2;

    double rayOriginDot = DOT(ray->R0, ray->R0);
    double centerDot = DOT(center, center);
    double rayOriginCenterDot = DOT(ray->R0, center);
    double radiusSquared = radius * radius;

    c = rayOriginDot + centerDot - 2.0 * rayOriginCenterDot - radiusSquared;


    double discriminant = b * b - 4.0 * a * c;

    if(discriminant < 0.0)
    {
        tMin = INF;
    }
    else
    {
        double sqrtDiscriminant = sqrt(discriminant);
        double denominator = 2.0 * a;

        double t1 = (-b - sqrtDiscriminant) / denominator;
        double t2 = (-b + sqrtDiscriminant) / denominator;

//            tMin = min(t1, t2);
        if(t1 < 0)
            tMin = t2;
        else
            tMin = t1;
    }

    if(level == 0)
    {
        return tMin;
    }

    //Illumination with the Phong Lighting Model
    Color intersectionPointColor = getColor();
    Vector3D intersectionPoint =  (ray->Rd * tMin) + ray->R0 ;


    //ambient component of reflected ray
    // Store the common ambient reflection coefficient in a variable
    double ambientReflection = reflectionCoefficient.ambient;

// Apply the ambient reflection to each color component
    color.red = intersectionPointColor.red * ambientReflection;
    color.green = intersectionPointColor.green * ambientReflection;
    color.blue = intersectionPointColor.blue * ambientReflection;


    //calculate normal at intersectionPoint
    Vector3D normal = intersectionPoint - center;
    normal.normalize();
    if(radius > distanceBetweenPoints(pos, center) )
    {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
    }

    pointLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);
    spotLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);

    //if level ≥ recursion_level, return tmin
    if( levelOfRecursion <= level )
    {
        return tMin;
    }

    //construct reflected ray from intersection point
    double dotProduct = DOT(ray->Rd, normal);
double scaledDotProduct = dotProduct * 2.0;
Vector3D scaledNormal = normal * scaledDotProduct;
Vector3D reflectionDirection = ray->Rd - scaledNormal;

    reflectionDirection.normalize();
    Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

    //find tmin from the nearest intersecting object, using
    //intersect() method, as done in the capture() method
    //if found, call intersect(rreflected, colorreflected, level+1)
    double tMinimum = INF;
    int nearest = findNearestObject(reflectedRay, tMinimum);

    // colorreflected will be updated while in the subsequent call
    // update color using the impact of reflection
    Color reflectedColor;
    if(nearest != INF)
    {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

// Calculate the contribution for each color channel
    double redContribution = reflectedColor.red * reflectionCoefficient.recursive;
    double greenContribution = reflectedColor.green * reflectionCoefficient.recursive;
    double blueContribution = reflectedColor.blue * reflectionCoefficient.recursive;

// Update each color component explicitly without using shorthand operators
    color.red = color.red + redContribution;
    color.green = color.green + greenContribution;
    color.blue = color.blue + blueContribution;


    color.clipColor();

    return tMin;
}


class Triangle: public Object
{
public:
    Vector3D a, b, c;

    Triangle(Vector3D a, Vector3D b, Vector3D c);
    void draw();
    double intersect(Ray* ray, Color& color, int level);
};

// Constructor definition
Triangle::Triangle(Vector3D a, Vector3D b, Vector3D c) : a(a), b(b), c(c) {}

// Draw method definition
void Triangle::draw()
{
    glColor3f(getColor().red, getColor().green, getColor().blue);
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
    }
    glEnd();
}
double Triangle::intersect(Ray* ray, Color& color, int level)
{
    double detA, detBeta, detGhama, detT, tMin,R0x, R0y, R0z,Rdx, Rdy, Rdz;
    R0x = ray->R0.x;
    R0y = ray->R0.y;
    R0z = ray->R0.z;
    Rdx = ray->Rd.x;
    Rdy = ray->Rd.y;
    Rdz = ray->Rd.z;

    detA = (a.x - b.x) * ((a.y - c.y) * Rdz - (a.z - c.z) * Rdy)
           +(a.x - c.x) * ((a.z - b.z) * Rdy - (a.y - b.y) * Rdz)
           +Rdx * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));

    detBeta = (a.x - R0x) * ((a.y - c.y) * Rdz - (a.z - c.z) * Rdy)
              +(a.x - c.x) * ((a.z - R0z) * Rdy - (a.y - R0y) * Rdz)
              +Rdx * ((a.y - R0y) * (a.z - c.z) - (a.z - R0z) * (a.y - c.y));

    detGhama = (a.x - b.x) * ((a.y - R0y) * Rdz - (a.z - R0z) * Rdy)
               +(a.x - R0x) * ((a.z - b.z) * Rdy - (a.y - b.y) * Rdz)
               +Rdx * ((a.y - b.y) * (a.z - R0z) - (a.z - b.z) * (a.y - R0y));

    detT = (a.x - b.x) * ((a.y - c.y) * (a.z - R0z) - (a.z - c.z) * (a.y - R0y))
           +(a.x - c.x) * ((a.z - b.z) * (a.y - R0y) - (a.y - b.y) * (a.z - R0z))
           +(a.x - R0x) * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));

    tMin = INF;
    if(detA != 0.0)
    {
        double beta = detBeta / detA;
        double ghama = detGhama / detA;
        double t = detT / detA;
        if(beta + ghama < 1.0 and beta > 0.0 && ghama > 0.0)
        {
            tMin = t;
        }
    }

    if(level == 0)
    {
        return tMin;
    }

    //Illumination with the Phong Lighting Model
    Vector3D intersectionPoint = ray->R0 + ray->Rd * tMin;
    Color intersectionPointColor = getColor();

    //ambient component of reflected ray
    // Calculate the common ambient reflection coefficient multiplier
    double ambientMultiplier = reflectionCoefficient.ambient;

// Update each color component explicitly without using shorthand operators
    color.red = intersectionPointColor.red * ambientMultiplier;
    color.green = intersectionPointColor.green * ambientMultiplier;
    color.blue = intersectionPointColor.blue * ambientMultiplier;


    //calculate normal
    Vector3D normal = CROSS((b - a), (c - a));
    normal.normalize();
    if(DOT(ray->Rd, normal) > 0.0)
    {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
    }

    pointLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);
    spotLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);

    //if level ≥ recursion_level, return tmin
    if(levelOfRecursion <= level)
    {
        return tMin;
    }

    //construct reflected ray from intersection point
    Vector3D reflectionDirection = ray->Rd - normal * (DOT(ray->Rd, normal) * 2.0);
    reflectionDirection.normalize();
    Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

    //find tmin from the nearest intersecting object, using
    //intersect() method, as done in the capture() method
    //if found, call intersect(rreflected, colorreflected, level+1)
    double tMinimum = INF;
    int nearest = findNearestObject(reflectedRay, tMinimum);

    // colorreflected will be updated while in the subsequent call
    // update color using the impact of reflection
    Color reflectedColor;

    if(nearest != INF)
    {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    // Calculate the reflection contribution for each color component
    double reflectedRedContribution = reflectedColor.red * reflectionCoefficient.recursive;
    double reflectedGreenContribution = reflectedColor.green * reflectionCoefficient.recursive;
    double reflectedBlueContribution = reflectedColor.blue * reflectionCoefficient.recursive;

// Explicitly update each color component without using shorthand operators
    color.red = color.red + reflectedRedContribution;
    color.green = color.green + reflectedGreenContribution;
    color.blue = color.blue + reflectedBlueContribution;


    color.clipColor();

    return tMin;
}


class GeneralQuadricSurfaceCoefficient
{
public:
    double a, b, c, d, e, f, g, h, i, j;
};

class GeneralQuadricSurface: public Object
{
public:
    GeneralQuadricSurfaceCoefficient coefficient;
    Vector3D cubeReferencePoint;
    double length,width,height;
    GeneralQuadricSurface()
    {
        coefficient.a =0.0;
        coefficient.b =0.0;
        coefficient.c =0.0;
        coefficient.d =0.0;
        coefficient.e = 0.0;
        coefficient.f = 0.0;
        coefficient.g = 0.0;
        coefficient.h = 0.0;
        coefficient.i = 0.0;
        coefficient.j = 0.0;
        length =0.0;
        width = 0.0;
        height = 0.0;
    }

    GeneralQuadricSurface(GeneralQuadricSurfaceCoefficient coefficient, Vector3D cubeReferencePoint, double length, double width, double height)
        : coefficient(coefficient), cubeReferencePoint(cubeReferencePoint), length(length), width(width), height(height)
    {

    }


    void draw()
    {
        ;
    }

    double clipping(double tMin, double tMax, Ray* ray)
    {
        if(tMin < INF)
        {
            if(tMax < INF)
            {
                if(tMin > 0.0)
                {
                    Vector3D intersectionPoint = ray->R0 + ray->Rd * tMin;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x or intersectionPoint.x>cubeReferencePoint.x+length)))
                    {
                        tMin = INF;
                    }
                    else if(width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y or intersectionPoint.y>cubeReferencePoint.y+width))
                    {
                        tMin = INF;
                    }
                    else if(height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z or intersectionPoint.z>cubeReferencePoint.z+height))
                    {
                        tMin = INF;
                    }
                }
                if(tMax > 0.0)
                {
                    Vector3D intersectionPoint;
                    intersectionPoint = ray->R0 + ray->Rd * tMax;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x or intersectionPoint.x>cubeReferencePoint.x+length)))
                    {
                        tMax = INF;
                    }
                    else if(width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y or intersectionPoint.y>cubeReferencePoint.y+width))
                    {
                        tMax = INF;
                    }
                    else if(height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z or intersectionPoint.z>cubeReferencePoint.z+height))
                    {
                        tMax = INF;
                    }
                }
                if(tMin > 0.0 and tMin < tMax)
                {
                    ;
                }
                else
                {
                    tMin = tMax;
                }
            }
            else
            {
                if(tMin > 0.0)
                {
                    Vector3D intersectionPoint = ray->R0 + ray->Rd * tMin;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x or intersectionPoint.x>cubeReferencePoint.x+length)))
                    {
                        tMin = INF;
                    }
                    else if((width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y or intersectionPoint.y>cubeReferencePoint.y+width)))
                    {
                        tMin = INF;
                    }
                    else if((height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z or intersectionPoint.z>cubeReferencePoint.z+height)))
                    {
                        tMin = INF;
                    }
                }
            }
        }
        return tMin;
    }

    double intersect(Ray* ray, Color& color, int level)
    {
        double a, b, c;

        double xSquaredTerm = coefficient.a * ray->Rd.x * ray->Rd.x;
        double ySquaredTerm = coefficient.b * ray->Rd.y * ray->Rd.y;
        double zSquaredTerm = coefficient.c * ray->Rd.z * ray->Rd.z;

        double xyTerm = coefficient.d * ray->Rd.x * ray->Rd.y;
        double xzTerm = coefficient.e * ray->Rd.x * ray->Rd.z;
        double yzTerm = coefficient.f * ray->Rd.y * ray->Rd.z;

        a = xSquaredTerm + ySquaredTerm + zSquaredTerm + xyTerm + xzTerm + yzTerm;

        double axTerm1 = 2.0 * coefficient.a * ray->R0.x * ray->Rd.x;
        double byTerm1 = 2.0 * coefficient.b * ray->R0.y * ray->Rd.y;
        double czTerm1= 2.0 * coefficient.c * ray->R0.z * ray->Rd.z;

        double dxyTerm1 = coefficient.d * (ray->R0.x * ray->Rd.y + ray->Rd.x * ray->R0.y);
        double exzTerm1 = coefficient.e * (ray->R0.x * ray->Rd.z + ray->Rd.x * ray->R0.z);
        double fyzTerm1 = coefficient.f * (ray->R0.y * ray->Rd.z + ray->Rd.y * ray->R0.z);

        double gxTerm1 = coefficient.g * ray->Rd.x;
        double hyTerm1 = coefficient.h * ray->Rd.y;
        double izTerm1 = coefficient.i * ray->Rd.z;

        b = axTerm1 + byTerm1 + czTerm1 + dxyTerm1 + exzTerm1 + fyzTerm1 + gxTerm1 + hyTerm1 + izTerm1;

        double axSquaredTerm = coefficient.a * ray->R0.x * ray->R0.x;
        double bySquaredTerm = coefficient.b * ray->R0.y * ray->R0.y;
        double czSquaredTerm = coefficient.c * ray->R0.z * ray->R0.z;

        double dxyTerm = coefficient.d * ray->R0.x * ray->R0.y;
        double exzTerm = coefficient.e * ray->R0.x * ray->R0.z;
        double fyzTerm = coefficient.f * ray->R0.y * ray->R0.z;

        double gxTerm = coefficient.g * ray->R0.x;
        double hyTerm = coefficient.h * ray->R0.y;
        double izTerm = coefficient.i * ray->R0.z;

        double constantTerm = coefficient.j;

        c = axSquaredTerm + bySquaredTerm + czSquaredTerm + dxyTerm + exzTerm + fyzTerm + gxTerm + hyTerm + izTerm + constantTerm;


        double tMin, tMax;
        if (a == 0.0)
        {
            tMin = (b == 0.0) ? INF : -c / b;
            tMax = INF;
        }
        else
        {
            double discriminant = b * b - 4.0 * a * c;

            if (discriminant < 0.0)
            {
                tMin = tMax = INF;
            }
            else
            {
                double sqrtDiscriminant = sqrt(discriminant);
                double divisor = 2.0 * a;
                tMax = (-b + sqrtDiscriminant) / divisor;
                tMin = (-b - sqrtDiscriminant) / divisor;
            }
        }


        tMin = clipping(tMin, tMax, ray);

        if(level == 0)
        {
            return tMin;
        }
        Color intersectionPointColor = getColor();
        Vector3D intersectionPoint = ray->R0 + ray->Rd * tMin;


        //ambient component
        color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
        color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
        color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;

        double xn, yn, zn;
        // Contributions to xn
        double axXContribution = 2.0 * coefficient.a * intersectionPoint.x;
        double dYContributionForX = coefficient.d * intersectionPoint.y;
        double eZContributionForX = coefficient.e * intersectionPoint.z;
        double gConstantForX = coefficient.g;
        xn = axXContribution + dYContributionForX + eZContributionForX + gConstantForX;

// Contributions to yn
        double byYContribution = 2.0 * coefficient.b * intersectionPoint.y;
        double dXContributionForY = coefficient.d * intersectionPoint.x;
        double fZContributionForY = coefficient.f * intersectionPoint.z;
        double hConstantForY = coefficient.h;
        yn = byYContribution + dXContributionForY + fZContributionForY + hConstantForY;

// Contributions to zn
        double czZContribution = 2.0 * coefficient.c * intersectionPoint.z;
        double eXContributionForZ = coefficient.e * intersectionPoint.x;
        double fYContributionForZ = coefficient.f * intersectionPoint.y;
        double iConstantForZ = coefficient.i;
        zn = czZContribution + eXContributionForZ + fYContributionForZ + iConstantForZ;


        Vector3D normal(xn, yn, zn);
        normal.normalize();

        pointLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);
        spotLightDiffuseAndSpecular(shine,ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);

        if(level >= levelOfRecursion)
        {
            return tMin;
        }

        double dotProduct = DOT(ray->Rd, normal);
        double scale = 2.0 * dotProduct;
        Vector3D scaledNormal = normal * scale;
        Vector3D reflectionDirection = ray->Rd - scaledNormal;

        reflectionDirection.normalize();
        Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

        double tMinimum = INF;
        int nearest = findNearestObject(reflectedRay, tMinimum);

        Color reflectedColor;

        if(nearest != INF)
        {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
        }

        color.red = color.red + reflectedColor.red*reflectionCoefficient.recursive;
        color.green = color.green + reflectedColor.green*reflectionCoefficient.recursive;
        color.blue = color.blue + reflectedColor.blue*reflectionCoefficient.recursive;

        color.clipColor();

        return tMin;
    }
};


class Floor : public Object
{
public:
    double floorWidth;
    double tileWidth;
    Color deepColor;

    // Constructors
    Floor();
    Floor(double floorWidth, double tileWidth);

    // Member functions
    void draw();
    double intersect(Ray* ray, Color& color, int level);
};

Floor::Floor() : floorWidth(0.0), tileWidth(0.0) {}

// Parameterized constructor
Floor::Floor(double floorWidth, double tileWidth)
    : floorWidth(floorWidth), tileWidth(tileWidth) {}



void Floor::draw()
{
    int tiles = floorWidth / tileWidth;
    double length = tileWidth;
    Vector3D reference_point = Vector3D(-(floorWidth / 2), -(floorWidth / 2), 0);

    for (int i = 0; i < tiles; i++)
    {
        for (int j = 0; j < tiles; j++)
        {
            if (((i + j) % 2) == 0) glColor3f(1, 1, 1);
            else glColor3f(0, 0, 0);

            glBegin(GL_QUADS);
            glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0);
            glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + j * length, 0);
            glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + (j + 1) * length, 0);
            glVertex3f(reference_point.x + i * length, reference_point.y + (j + 1) * length, 0);
            glEnd();
        }
    }
}

// intersect function
double Floor::intersect(Ray* ray, Color& color, int level)
{
    double xn = 0.0;
    double yn = 0.0;
    double zn = 1.0;
    if (pos.z < 0) zn = -1;
    Vector3D normal(xn, yn, zn);

    double tMin = INF,D = 0,nR0,nRd;
    nR0 = DOT(normal, ray->R0);
    nRd = DOT(normal, ray->Rd);
    if (nRd != 0.0)
    {
        tMin = -(D + nR0) / nRd;
    }

    if (level == 0)
    {
        return tMin;
    }

    Color intersectionPointColor;
    Vector3D intersectionPoint = ray->R0 + ray->Rd * tMin;
    Vector3D referencePosition = intersectionPoint - Vector3D(-floorWidth / 2.0, -floorWidth / 2.0, 0.0);
    // Calculate the condition outside the if statement
    int tileSum = (int)(floor(referencePosition.x / tileWidth) + floor(referencePosition.y / tileWidth));
    bool isEvenTile = (tileSum % 2) == 0;

// Use the pre-calculated condition in the if statement
    if (isEvenTile)
    {
        intersectionPointColor = getColor();
    }
    else
    {
        intersectionPointColor = deepColor;
    }


    // Ambient component
    color.red = intersectionPointColor.red * reflectionCoefficient.ambient;
    color.green = intersectionPointColor.green * reflectionCoefficient.ambient;
    color.blue = intersectionPointColor.blue * reflectionCoefficient.ambient;

    pointLightDiffuseAndSpecular(shine, ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);
    spotLightDiffuseAndSpecular(shine, ray, intersectionPoint, normal, intersectionPointColor, color, reflectionCoefficient);

    if (level >= levelOfRecursion)
    {
        return tMin;
    }

    Vector3D reflectionDirection = ray->Rd - normal * (DOT(ray->Rd, normal) * 2.0);
    reflectionDirection.normalize();
    Ray* reflectedRay = new Ray(intersectionPoint + reflectionDirection, reflectionDirection);

    double tMinimum = INF;
    int nearest = findNearestObject(reflectedRay, tMinimum);


    Color reflectedColor;
    if (nearest != INF)
    {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level + 1);
    }
// Calculate the reflection contribution for each color component
    double reflectedRedContribution = reflectedColor.red * reflectionCoefficient.recursive;
    double reflectedGreenContribution = reflectedColor.green * reflectionCoefficient.recursive;
    double reflectedBlueContribution = reflectedColor.blue * reflectionCoefficient.recursive;

// Explicitly update each color component without using shorthand operators
    color.red = color.red + reflectedRedContribution;
    color.green = color.green + reflectedGreenContribution;
    color.blue = color.blue + reflectedBlueContribution;


    color.clipColor();

    return tMin;
}
