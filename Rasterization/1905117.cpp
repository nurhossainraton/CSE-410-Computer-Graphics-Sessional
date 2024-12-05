#include <bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;

#define pb push_back
#define PI 2.0*acos(0.0)

class Matrix
{

public:
    double m[4][4];
    double triangle[4];
    Matrix();
    void printMatrix();
    void multiply(Matrix a);
    void multiply(Matrix a,Matrix b);
    vector<double>multiply(vector<double> v,Matrix a);
    Matrix multiply2(Matrix a,Matrix b);
    vector<double> RodriguesValue(vector<double>x,vector<double>a,double angle);
    vector<double>normalize(vector<double> v);
};
Matrix::Matrix()
{
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            if(i!=j)
                m[i][j]=0;
            else
                m[i][j]=1;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        triangle[i] = 0;
    }
}
vector<double> Matrix:: normalize(vector<double> v)
{
    double val =0;
    vector<double>res(3,0);
    val=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    res[0]=v[0]/val;
    res[1]=v[1]/val;
    res[2]=v[2]/val;

    return res;

}

vector<double> Matrix:: RodriguesValue(vector<double>x,vector<double>a,double angle)
{
    vector<double>result(3,0);
    result[0]=x[0]*cos(angle)+(1-cos(angle))*(a[0]*a[0]*x[0]+a[0]*a[1]*x[1]+a[0]*a[2]*x[2])+sin(angle)*(a[1]*x[2]-a[2]*x[1]);
    result[1]=x[1]*cos(angle)+(1-cos(angle))*(a[0]*a[1]*x[0]+a[1]*a[1]*x[1]+a[1]*a[2]*x[2])+sin(angle)*(a[2]*x[0]-a[0]*x[2]);
    result[2]=x[2]*cos(angle)+(1-cos(angle))*(a[0]*a[2]*x[0]+a[1]*a[2]*x[1]+a[2]*a[2]*x[2])+sin(angle)*(a[0]*x[1]-a[1]*x[0]);
    return result;
}

void Matrix::printMatrix()
{
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            cout<<m[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void Matrix::multiply(Matrix a)
{
    double temp[4],sum=0;
    for (int i = 0; i < 4; i++)
    {
        sum = 0;
        for (int j = 0; j < 4; j++)
        {
            sum = sum+(triangle[j]*a.m[i][j]);
        }
        temp[i] = sum;
    }

    for (int i = 0; i < 4; i++)
    {
        triangle[i] = temp[i];
    }
}
Matrix identity,translation,scaling,rotation;

void Matrix::multiply(Matrix a,Matrix b)
{
    double sum=0;
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            sum = 0;
            for(int k=0; k<4; k++)
            {
                sum = sum+(a.m[i][k]*b.m[k][j]);
            }
            identity.m[i][j] = sum;
        }
    }
}
vector<double> Matrix::multiply(vector<double> v,Matrix a)
{
    vector<double> d(4,0);
    double sum = 0;
    for(int i=0; i<4; i++)
    {
        sum=0;
        for(int j=0; j<4; j++)
        {
            sum = sum+(v[j]* a.m[i][j]);
        }
        d[i]=sum;
    }
    return d;

}

Matrix Matrix::multiply2(Matrix a,Matrix b)
{
    double sum=0;
    Matrix temp;
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            sum = 0;
            for(int k=0; k<4; k++)
            {
                sum = sum+(a.m[i][k]*b.m[k][j]);
            }
            temp.m[i][j] = sum;
        }
    }
    return temp;
}
std::string readAndTrimTrailingSpaces(std::string const & filePath)
{
    std::ifstream file(filePath);
    std::string   buffer(std::istreambuf_iterator<char> {file}, {});

    while (!buffer.empty() && std::isspace(buffer.back()))
        buffer.pop_back();

    return buffer;
}

class Triangle
{
public:
    double colors[3];
    double points[3][3];
};

static unsigned long int g_seed = 1;
inline int random()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

int main()
{
    Matrix mat;
    stack<Matrix>stk;
    cout<<fixed<<setprecision(7);
    double fovY,aspectRatio,near,far;
    vector<double>look(3,0),up(3,0),eye(3,0);
    string command;
    double x1, y1, z1,tx, ty, tz,sx, sy, sz,rx, ry, rz, angle;


    ifstream fp("scene.txt");
    fp>>eye[0]>>eye[1]>>eye[2];
    fp>>look[0]>>look[1]>>look[2];
    fp>>up[0]>>up[1]>>up[2];
    fp>>fovY>>aspectRatio>>near>>far;
    ofstream fp1("temp.txt");
    fp1<<fixed<<setprecision(7);

    while(!fp.eof())
    {
        fp>>command;
        if(command == "triangle")
        {
            for(int i=0; i<3; i++)
            {
                fp>>x1>>y1>>z1;

                mat.triangle[0]=x1;
                mat.triangle[1]=y1;
                mat.triangle[2]=z1;
                mat.triangle[3]=1;

                mat.multiply(identity);

                if(mat.triangle[3]!=1)
                {
                    mat.triangle[0]=mat.triangle[0]/mat.triangle[3];
                    mat.triangle[1]=mat.triangle[1]/mat.triangle[3];
                    mat.triangle[2]=mat.triangle[2]/mat.triangle[3];

                }
                for(int j=0; j<3; j++)
                {
                    fp1<<mat.triangle[j]<<" ";
                }
                fp1<<endl;

            }
            fp1<<endl;

        }
        else if (command == "translate")
        {
            fp>>tx>>ty>>tz;
            translation.m[0][3]=tx;
            translation.m[1][3]=ty;
            translation.m[2][3]=tz;

            mat.multiply(identity,translation);
            //translation.printMatrix();
        }

        else if (command == "scale")
        {
            fp>>sx>>sy>>sz;
            scaling.m[0][0]=sx;
            scaling.m[1][1]=sy;
            scaling.m[2][2]=sz;
            mat.multiply(identity,scaling);
            //scaling.printMatrix();
        }
        else if(command == "rotate")
        {
            fp>>angle>>rx>>ry>>rz;
            angle = PI/180*angle;
            double nx,ny,nz;
            vector<double>a,i,j,k;

            nx=rx/sqrt(rx*rx+ry*ry+rz*rz);
            ny=ry/sqrt(rx*rx+ry*ry+rz*rz);
            nz=rz/sqrt(rx*rx+ry*ry+rz*rz);

            a.pb(nx);
            a.pb(ny);
            a.pb(nz);

            i.pb(1);
            i.pb(0);
            i.pb(0);

            vector<double>c1= mat.RodriguesValue(i,a,angle);

            j.pb(0);
            j.pb(1);
            j.pb(0);

            vector<double>c2= mat.RodriguesValue(j,a,angle);

            k.pb(0);
            k.pb(0);
            k.pb(1);
            vector<double>c3= mat.RodriguesValue(k,a,angle);

            rotation.m[0][0] = c1[0];
            rotation.m[1][0] = c1[1];
            rotation.m[2][0] = c1[2];

            rotation.m[0][1] = c2[0];
            rotation.m[1][1] = c2[1];
            rotation.m[2][1] = c2[2];

            rotation.m[0][2] = c3[0];
            rotation.m[1][2] = c3[1];
            rotation.m[2][2] = c3[2];

            mat.multiply(identity,rotation);

            //rotation.printMatrix();

        }
        else if(command == "push")
        {
            stk.push(identity);
        }
        else if(command == "pop")
        {
            identity = stk.top();
            stk.pop();
        }
        else if(command == "end")
        {
            break;
        }


    }
    fp.close();
    fp1.close();

    string temporary = readAndTrimTrailingSpaces("temp.txt");
    ofstream fpt;
    fpt.open("stage1.txt");
    fpt << fixed << setprecision(7);
    fpt << temporary;
    fpt.close();
    remove("temp.txt");

    vector<double>L(3,0);
    L[0]=look[0]-eye[0];
    L[1]=look[1]-eye[1];
    L[2]=look[2]-eye[2];
    L=mat.normalize(L);

    vector<double>R(3,0);
    R[0]=L[1]*up[2]-L[2]*up[1];
    R[1]=L[2]*up[0]-L[0]*up[2];
    R[2]=L[0]*up[1]-L[1]*up[0];
    R =mat.normalize(R);

    vector<double>U(3,0);
    U[0]=R[1]*L[2]-R[2]*L[1];
    U[1]=R[2]*L[0]-R[0]*L[2];
    U[2]=R[0]*L[1]-R[1]*L[0];
    U = mat.normalize(U);

    Matrix T,RotationMat;

    RotationMat.m[0][0] = R[0];
    RotationMat.m[0][1] = R[1];
    RotationMat.m[0][2] = R[2];
    RotationMat.m[1][0] = U[0];
    RotationMat.m[1][1] = U[1];
    RotationMat.m[1][2] = U[2];
    RotationMat.m[2][0] = -L[0];
    RotationMat.m[2][1] = -L[1];
    RotationMat.m[2][2] = -L[2];

    T.m[0][3] = -eye[0];
    T.m[1][3] = -eye[1];
    T.m[2][3] = -eye[2];

    Matrix v = mat.multiply2(RotationMat,T);

    ifstream fp2("stage1.txt");
    ofstream fp3("temp.txt");

    fp3<<fixed << setprecision(7);

    int linecount=0;
    vector<double> arr(4,0);

    while(!fp2.eof())
    {
        vector<double>newV;
        double x1,y1,z1;
        fp2>>x1>>y1>>z1;

        arr[0]=x1;
        arr[1]=y1;
        arr[2]=z1;
        arr[3]=1;

        linecount++;

        newV= mat.multiply(arr,v);
        double first = newV[0]/newV[3];
        double s = newV[1]/newV[3];
        double th = newV[2]/newV[3];

        fp3<<first<<" "<<s<<" "<<th<<endl;

        if(linecount==3)
        {
            linecount=0;
            fp3<<endl;
        }


    }
    fp2.close();
    fp3.close();

    temporary = readAndTrimTrailingSpaces("temp.txt");
    fpt.open("stage2.txt");
    fpt << fixed << setprecision(7);
    fpt << temporary;
    fpt.close();
    remove("temp.txt");

    ifstream fp4("stage2.txt");
    ofstream fp5("stage3.txt");
    fp5 << fixed << setprecision(7);

    double fovX,tt,rr;

    fovX=fovY*aspectRatio;
    tt= near*tan(fovY*PI/360.0);
    rr= near*tan(fovX*PI/360.0);

    Matrix P;

    P.m[0][0] = near / rr;
    P.m[1][1] = near / tt;
    P.m[2][2] = -(far + near) / (far - near);
    P.m[2][3] = -(2 * far * near) / (far - near);
    P.m[3][2] = -1;
    P.m[3][3] = 0;

    int linecount1=0;
    vector<double> arr1(4,0);

    while(!fp4.eof())
    {
        vector<double>newV;
        double x1,y1,z1;
        fp4>>x1>>y1>>z1;

        arr1[0]=x1;
        arr1[1]=y1;
        arr1[2]=z1;
        arr1[3]=1;

        linecount1++;

        newV= mat.multiply(arr1,P);
        double first = newV[0]/newV[3];
        double s = newV[1]/newV[3];
        double th = newV[2]/newV[3];

        fp5<<first<<" "<<s<<" "<<th<<endl;

        if(linecount1==3)
        {
            linecount1=0;
            fp5<<endl;
        }


    }
    fp4.close();
    fp5.close();

    int width,height;
    ifstream fp6("config.txt");
    ofstream fp8("z_buffer.txt");
    fp8 << fixed << setprecision(7);
    fp6>>width>>height;
    fp6.close();
    fp6.open("stage3.txt");

    double boxR=1,boxL=-1,boxT=1,boxB=-1,zmax=1.0;

    double dx=(boxR-boxL)/width;
    double dy= (boxT-boxB)/height;
    double topY = boxT - (dy / 2);
    double bottomY = boxB + (dy / 2);
    double leftX = boxL + (dx / 2);
    double rightX = boxR - (dx / 2);

    bitmap_image *image;
    double **z_buffer = new double *[height];
    image = new bitmap_image(width,height);

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            image->set_pixel(j,i,0,0,0);

        }
    }

    for(int i=0; i<height; i++)
    {
        z_buffer[i]= new double [width];
    }

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            z_buffer[i][j] = zmax;

        }
    }
    while(!fp6.eof())
    {
        Triangle t;
        for(int i=0; i<3; i++)
        {
            t.colors[i]= random()%256;
        }
        for(int i=0; i<3; i++)
        {
            fp6>>t.points[i][0]>>t.points[i][1]>>t.points[i][2];

        }

        double xmax,xmin,ymax,ymin;
        xmax = max(t.points[0][0],max(t.points[1][0],t.points[2][0]));
        xmin = min(t.points[0][0],min(t.points[1][0],t.points[2][0]));
        ymax = max(t.points[0][1],max(t.points[1][1],t.points[2][1]));
        ymin = min(t.points[0][1],min(t.points[1][1],t.points[2][1]));

        xmin = max(xmin,leftX);
        ymin = max(ymin,bottomY);
        ymax = min(ymax,topY);
        xmax = min(xmax,rightX);

        int startY = round((topY-ymin)/dy);
        int endY = round((topY-ymax)/dy);

        double a,b,c,d,a1,b1,c1,a2,b2,c2;
        a1 = t.points[1][0] - t.points[0][0];
        b1 = t.points[1][1] - t.points[0][1];
        c1 = t.points[1][2] - t.points[0][2];
        a2 = t.points[2][0] - t.points[0][0];
        b2 = t.points[2][1] - t.points[0][1];
        c2 = t.points[2][2] - t.points[0][2];
        a = b1 * c2 - b2 * c1;
        b = a2 * c1 - a1 * c2;
        c = a1 * b2 - b1 * a2;
        d = (- a * t.points[0][0] - b * t.points[0][1] - c * t.points[0][2]);

        for(int scanY = 1+endY; scanY<=startY; scanY++)
        {
            double y,ara[2];
            int count = 0;
            y=topY-(scanY*dy);

            int i,j;
            for(i=0; i<3; i++)
            {
                j=(i+1)%3;
                if(t.points[i][1] == t.points[j][1])
                {
                    if(y==t.points[i][1])
                    {
                        double minval= min(t.points[i][0], t.points[j][0]);
                        double maxval = min(t.points[i][0], t.points[j][0]);
                        ara[count] = minval;
                        ara[count+1] = maxval;
                    }
                    continue;
                }

                double minim = min(t.points[i][1], t.points[j][1]);
                double maxm = max(t.points[i][1], t.points[j][1]);
                if(y >=minim  && y <= maxm)
                {
                    ara[count] = t.points[i][0] + (y - t.points[i][1]) * (t.points[j][0] - t.points[i][0]) / (t.points[j][1] - t.points[i][1]);
                    count++;
                }
            }
            for(int i=0; i<2; i++)
            {
                if(ara[i]<xmin) ara[i] = xmin;
                if(ara[i]>xmax) ara[i] = xmax;
            }
            if(ara[0]>=ara[1])
            {
                swap(ara[0],ara[1]);
            }
            int startX,endX;
            startX = round((ara[0] - leftX) / dx);
            endX = round((ara[1] - leftX) / dx);
            int scanX;
            for(scanX = startX; scanX<endX; scanX++)
            {
                double x,z;
                z=(-d-(a*x)-(b*y))/c;
                x=(scanX*dx)+leftX;

                if(z<-1) continue;
                if(z<z_buffer[scanY][scanX])
                {
                    z_buffer[scanY][scanX] = z;
                    image->set_pixel(scanX,scanY,t.colors[0],t.colors[1],t.colors[2]);


                }


            }
        }


    }
    fp6.close();
    image->save_image("out.bmp");
    int j;
    for(int i=0; i<height; i++)
    {
        for(j=0; j<width; j++)
        {
            if(z_buffer[i][j]<zmax)
            {
                fp8<<z_buffer[i][j]<<"\t";
            }

        }
        fp8<<endl;
    }
    fp8.close();

    for(int i=0; i<height; i++)
    {
        delete[] z_buffer[i];
    }
    delete []z_buffer;
    delete image;
    return 0;



}
