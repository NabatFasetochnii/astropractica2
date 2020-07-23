#include "MyForm.h"
#include <Windows.h>
#include <Math.h>
#include <iostream>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <Math.h>


using namespace astropractica;

using namespace System::Drawing;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Collections::Generic;
using namespace System;
using namespace Eigen;
using namespace System::Drawing::Drawing2D;


int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
    Application::EnableVisualStyles();
    Application::SetCompatibleTextRenderingDefault(false);
    Application::Run(gcnew MyForm);
    return 0;
}


//ref class BigCircle //класс большого круга
//{
//public:
//    int getCountOfPoint(); //получить количество точек отрисовки круга 
//    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis); //конструкторы
//    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots);
//    BigCircle(Graphics^ graphics, Pen^ pen);
//    BigCircle(Graphics^ graphics);
//    ~BigCircle(); //дискриптор 
//    void rotationX(double u); //ворочаем по осям на угол u (в радианах)
//    void rotationY(double u);
//    void rotationZ(double u);
//    void onDrowAll(); //нарисовать всё
//    array<Vector3d*, 1>^ getCircleCoords(); //получить массив векторов, соответствующих точкам круга
//    array<Vector3d*, 1>^ getAxis(); //получить массив векторов, соответствующих точкам оси круга
//    Vector3d getNormal(); //получить вектор нормали 
//    Vector3d getZenit();  //получить точку пересечения оси с небесной сферой 1
//    Vector3d getNadir();  //получить точку пересечения оси с небесной сферой 2
//    bool addDot(Vector3d v); //добавить на круг точку 
//    bool isDotOnCircle(Vector3d v); //проверить находится ли точка на круге
//    int getCountOfDots(); //получить счётчик пользовательских точек 
//    array<Vector3d*, 1>^ getDots(); //получить массив пользовательских точек 
//    void drowCircle(); // нарисовать по отдельности круг
//    void drowAxis(); // ось
//    void drowDots();// точки
//    void drowZenit(); //точку пересечения 1
//    void drowNadir(); //точку пересечения 2
//    double dihedralAngle_rad(BigCircle b); //двугранный угол между двумя плоскостями больших кругов, в радианах
//    double dihedralAngle_deg(BigCircle b);//двугранный угол между двумя плоскостями больших кругов, в градусах
//    Vector3d intersectionDots(BigCircle b); //получить точку пересечения двух больших кругов, вторая точка является противоположной этой, то есть просто можно домножить на -1
//
//private:
//
//    Vector3d MX(Vector3d v, double u); //оперрации поворота вокруг соотвествующей оси
//    Vector3d MY(Vector3d v, double u); //"здесь поворот совершается против часовой стрелки
//    Vector3d MZ(Vector3d v, double u); //а объект по сути поворачивается в другую сторону, по часовой стрелке" (с) Аналитик
//    //на вход идёт вектор и угол в радианах, на который нужно повернуть
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Graphics^ graphics; // объект для отрисовки
//    Vector3d* normal; //нормаль к поверхности
//    array<Vector3d*, 1>^ coords = gcnew array<Vector3d*>(629); // массив координат точек круга //исправил
//    array<Vector3d*, 1>^ axis = gcnew array<Vector3d*>(629); //массив координат точек оси //исправил
//    array<Vector3d*, 1>^ dots = gcnew array<Vector3d*>(100); //массив координат пользовательских точек // исправил 
//    List<PointF>^ points = gcnew List<PointF>(); //коллекция точек большого круга, для отрисовки
//    List<PointF>^ pointsAxis = gcnew List<PointF>(); //коллекция точек оси, для отрисовки
//    double s = 0.01; //шаг по углу построения точек круга 
//    double alpha = 0; // угол поворота, нужен для построения круга
//    int tochki = 0; // колличество точек круга
//    int unit = 200; //массштаб
//    Pen^ pen; //ручка для круга
//    Pen^ penForDots; //ручка для точек
//    Pen^ penForAxis; //ручка для оси
//    Vector3d* zenit; //точка пересечения оси с небсферой 1
//    Vector3d* nadir;//точка пересечения оси с небсферой 2
//    int countOfDots = 0; //счётчик пользовательских точек
//    void init(); // построение круга
//    const double EPS = 0.1; //бесконечно малое эпсилон
//
//};
array<Vector3d*, 1>^ BigCircle::getCircleCoords() {
    return coords;
}
array<Vector3d*, 1>^ BigCircle::getAxis()
{
    return axis;
}
Vector3d BigCircle::getNormal()
{
    return *normal;
}
Vector3d BigCircle::getZenit()
{
    return *zenit;
}
Vector3d BigCircle::getNadir()
{
    return *nadir;
}
bool BigCircle::addDot(int x, int y)
{

    if (isDotOnCircle(x - pW / 2, -y + pH / 2)) {
        dots[countOfDots] = buf;
        countOfDots++;
        return true;
    }
    return false;
}

bool BigCircle::isDotOnCircle(int x, int y) {

    for each (Vector3d * c in coords) {

        if ((abs(x - c->data()[0])) < EPS && 
            (abs(y - c->data()[1])) < EPS) {
            buf = c;
            return true;
        }
    }

    return false;
}
bool BigCircle::isDotOnCircle(Vector3d v)
{
    for each (Vector3d *c in coords) {

        if ((*c - v).norm() < EPS) {
            return true;
        }
    }

    return false;
}
int BigCircle::getCountOfDots()
{
    return countOfDots;
}
array<Vector3d*, 1>^ BigCircle::getDots()
{
    return dots;
}
void BigCircle::drowCircle()
{
    for (int i = 0; i < COUNT_OF_ARR_POINTS-1; i++) {
    
    
        if (coords[i]->data()[2] > 0) {
            graphics->DrawLine(penForBack,
                PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0),
                PointF(coords[i + 1]->data()[0] + pW / 2.0, -coords[i + 1]->data()[1] + pH / 2.0));
        }
        else
        {
            graphics->DrawLine(pen,
                PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0),
                PointF(coords[i + 1]->data()[0] + pW / 2.0, -coords[i + 1]->data()[1] + pH / 2.0));
        }
    }
    //graphics->DrawLines(pen, points->ToArray());

}
void BigCircle::drowAxis()
{
    graphics->DrawLines(penForAxis, pointsAxis->ToArray());
}
void BigCircle::drowDots()
{
    if (countOfDots > 0) {
        for (int i = 0; i < countOfDots; i++) {
                graphics->DrawEllipse(penForDots, dots[i]->data()[0] + pW / 2.0, -dots[i]->data()[1] + pH / 2.0, 5, 5);
        }
    }
    
}
void BigCircle::drowZenit()
{
    graphics->DrawEllipse(penForDots, zenit->data()[0] + pW/2.0, -zenit->data()[1] + pH/2.0, 5, 5);
}
void BigCircle::drowNadir()
{
    graphics->DrawEllipse(penForDots, nadir->data()[0] + pW / 2.0, -nadir->data()[1] + pH / 2.0, 5, 5);
}

double BigCircle::dihedralAngle_deg(BigCircle b)
{
    return (acos(normal->dot(b.getNormal())) * 57.2958);
}
Vector3d BigCircle::intersectionDots(BigCircle b)
{
    for each (Vector3d *v in coords)
    {
        for each (Vector3d *v2 in b.getCircleCoords())
        {
            if ((*v - *v2).norm() < EPS) {
                
                return *v;
            }
        }
    }
    Vector3d a;
    return a;
}

double BigCircle::dihedralAngle_rad(BigCircle b)
{
    return acos(normal->dot(b.getNormal()));
}
void BigCircle::onDrowAll() {
    
    drowCircle();
    drowAxis();
    drowDots();
    drowZenit();
    drowNadir();

}
int BigCircle::getCountOfPoint() {
    return tochki;
}

void astropractica::BigCircle::drowArc()
{
    if (isArcSet)
    {
        for (int i = arcStatPoint; i < arcEndPoint; i++) {

            graphics->DrawLine(penForDots, 
                PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0), 
                PointF(coords[i+1]->data()[0] + pW / 2.0, -coords[i+1]->data()[1] + pH / 2.0));
        }

        isArcSet = false;
        /*graphics->DrawEllipse(penForDots, 
            coords[arcStatPoint]->data()[0] + pW / 2.0, -coords[arcStatPoint]->data()[1] + pH / 2.0,
            5, 5);
        graphics->DrawEllipse(penForDots, 
            coords[arcEndPoint]->data()[0] + pW / 2.0, -coords[arcEndPoint]->data()[1] + pH / 2.0,
            5, 5);*/
    }
}


void BigCircle::rotationX(double u) {

    //points->Clear();
    pointsAxis->Clear();
    
    for (int i = 0; i < COUNT_OF_ARR_POINTS; i++) {

        *coords[i] = MX(*coords[i], u);
        *axis[i] = MX(*axis[i], u);
        if (i < 100) {
            *dots[i] = MX(*dots[i], u);
        }
        
        //points->Add(PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0));
        pointsAxis->Add(PointF(axis[i]->data()[0] + pW / 2.0, -axis[i]->data()[1] + pH / 2.0));
        
    }

    *normal = MX(*normal, u);
    *zenit = MX(*zenit, u);
    *nadir = MX(*nadir, u);
}
void BigCircle::rotationY(double u) {

   // points->Clear();
    pointsAxis->Clear();

    for (int i = 0; i < COUNT_OF_ARR_POINTS; i++) {

        *coords[i] = MY(*coords[i], u);
        *axis[i] = MY(*axis[i], u);
        if (i < 100) {
            *dots[i] = MY(*dots[i], u);
        }
       // points->Add(PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0));
        pointsAxis->Add(PointF(axis[i]->data()[0] + pW / 2.0, -axis[i]->data()[1] + pH / 2.0));

    }

    *normal = MY(*normal, u);
    *zenit = MY(*zenit, u);
    *nadir = MY(*nadir, u);
}
void BigCircle::rotationZ(double u) {

    //points->Clear();
    pointsAxis->Clear();

    for (int i = 0; i < COUNT_OF_ARR_POINTS; i++) {

        *coords[i] = MZ(*coords[i], u);
        *axis[i] = MZ(*axis[i], u);
        if (i < 100) {
            *dots[i] = MZ(*dots[i], u);
        }
        
        //points->Add(PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0));
        pointsAxis->Add(PointF(axis[i]->data()[0] + pW / 2.0, -axis[i]->data()[1] + pH / 2.0));
    }

    *normal = MZ(*normal, u);
    *zenit = MZ(*zenit, u);
    *nadir = MZ(*nadir, u);
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis, int pW, int pH)
{
    this->graphics = graphics;
    this->pen = pen;
    this->penForDots = penForDots;
    this->penForAxis = penForAxis;
    this->pW = pW;
    this->pH = pH;

    init();
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, int pW, int pH)
{
    this->graphics = graphics;
    this->pen = pen;
    this->penForDots = penForDots;
    this->pW = pW;
    this->pH = pH;

    penForAxis = gcnew Pen(Color::Black);
    penForAxis->Width = 2;
    penForAxis->LineJoin = LineJoin::Round;
    
    init();
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, int pW, int pH)
{
    this->graphics = graphics;
    this->pen = pen;
    this->graphics = graphics;
    this->pW = pW;
    this->pH = pH;

    penForDots = gcnew Pen(Color::Red);
    penForDots->Width = 4;
    penForDots->LineJoin = LineJoin::Round;

    penForAxis = gcnew Pen(Color::Black);
    penForAxis->Width = 2;
    penForAxis->LineJoin = LineJoin::Round;
    init();
}
BigCircle::BigCircle(Graphics^ graphics, int pW, int pH)
{
    this->graphics = graphics;
    this->pW = pW;
    this->pH = pH;


    pen = gcnew Pen(Color::Black);
    pen->Width = 3;
    pen->LineJoin = LineJoin::Round;

    penForDots = gcnew Pen(Color::Red);
    penForDots->Width = 4;
    penForDots->LineJoin = LineJoin::Round;

    penForAxis = gcnew Pen(Color::Black);
    penForAxis->Width = 2;
    penForAxis->LineJoin = LineJoin::Round;

    init();
}
astropractica::BigCircle::BigCircle(Graphics^ graphics, int pW, int pH, int x1, int y1, int x2, int y2)
{
    int X1 = x1 - pW / 2;
    int X2 = x2 - pW / 2;
    int Y1 = (pH / 2) - y1;
    int Y2 = (pH / 2) - y2;

    this->graphics = graphics;
        this->pW = pW;
        this->pH = pH;

        pen = gcnew Pen(Color::Black);
        pen->Width = 3;
        pen->LineJoin = LineJoin::Round;

        penForDots = gcnew Pen(Color::Red);
        penForDots->Width = 4;
        penForDots->LineJoin = LineJoin::Round;

        penForAxis = gcnew Pen(Color::Black);
        penForAxis->Width = 2;
        penForAxis->LineJoin = LineJoin::Round;

        init();

        //metod1(X1, Y1, X2, Y2);

        if (sqrt(pow(X1, 2) + pow(Y1, 2)) < unit && sqrt(pow(X2, 2) + pow(Y2, 2)) < unit) {

            Vector3d* v1 = new Vector3d();
            Vector3d* v2 = new Vector3d();
            Vector3d* v3 = new Vector3d();

            *v1 << X1, Y1, -sqrt(pow(unit, 2.0) - pow(X1, 2.0) - pow(Y1, 2.0));
            *v2 << X2, Y2, -sqrt(pow(unit, 2.0) - pow(X2, 2.0) - pow(Y2, 2.0));

            *v3 = (v1->cross(*v2)).normalized();//нормаль

          /*  Vector3d* s1 = new Vector3d();

            *s1 = plane(*normal, *v3);*/
            double a = cordTETA(*v3);
            rotationY(-a);
            a = cordFI(*v3);
            rotationZ(-a);
            
        }
}

void astropractica::BigCircle::metod1(int X1, int Y1, int X2, int Y2)
{
if (sqrt(pow(X1, 2) + pow(Y1, 2)) < unit && sqrt(pow(X2, 2) + pow(Y2, 2)) < unit) {

        Vector3d * v1 = new Vector3d();
        Vector3d * v2 = new Vector3d();
        Vector3d * v3 = new Vector3d();

        *v1 << X1, Y1, -sqrt(pow(unit, 2.0) - pow(X1, 2.0) - pow(Y1, 2.0));
        *v2 << X2, Y2, -sqrt(pow(unit, 2.0) - pow(X2, 2.0) - pow(Y2, 2.0));

        *v3 = (v1->cross(*v2)).normalized();//нормаль

        Vector3d* s1 = new Vector3d();
        Vector3d* s2 = new Vector3d();

        if (v3->data()[0] == 0 && v3->data()[1] == 0) {
            
        }
        else if(v3->data()[1] == 0 && v3->data()[2] == 0)
        {
            rotationY(M_PI_2);
        }
        else if (v3->data()[0] == 0 && v3->data()[2] == 0) {
            rotationX(M_PI_2);
        }
        else
        {
        *s1 << v3->data()[0], 0, v3->data()[2];
        *s2 << normal->data()[0], 0, normal->data()[2];
        double y = (s1->dot(*s2));
        double t = acos(y);
        //double t = s1->dot(*s2);
        rotationY(t);

        *s1 << v3->data()[0], v3->data()[1], 0;
        *s2 << normal->data()[0], normal->data()[1], 0;
        t = acos(s1->dot(*s2));
        rotationZ(t);

        *s1 << 0, v3->data()[1], v3->data()[2];
        *s2 << 0, normal->data()[1], normal->data()[2];
        t = acos(s1->dot(*s2));
        rotationX(t);
        }

    } 
}
void astropractica::BigCircle::setArc(int X1, int Y1, int X2, int Y2)
{
    //Vector3d *v1 = new Vector3d();
    //Vector3d *v2 = new Vector3d();
    int x1 = X1 - pW / 2;
    int x2 = X2 - pW / 2;
    int y1 = (pH / 2) - Y1;
    int y2 = (pH / 2) - Y2;

    //*v1 << x1, y1, -sqrt(pow(unit, 2.0) - pow(x1, 2.0) - pow(y1, 2.0));//надо повернуть 
    //*v2 << x2, y2, -sqrt(pow(unit, 2.0) - pow(x2, 2.0) - pow(y2, 2.0));
    int a, b;
    bool q = false;
    bool w = false;

    for (int i = 0; i < COUNT_OF_ARR_POINTS; i++)
    {

        

        if (abs((coords[i]->data()[0] - x1)) < EPS && 
            abs(coords[i]->data()[1] - y1) < EPS)
        {
            a = i;
            q = true;
            continue;
        }

        if (abs((coords[i]->data()[0] - x2)) < EPS &&
            abs((coords[i]->data()[1] - y2)) < EPS)
        {
            b = i;
            w = true;
            continue;
        }
    }
    
    if (q && w)
            {
            isArcSet = true;
            if (a<b)
                {
                    arcStatPoint = a;
                    arcEndPoint = b;
                }
            else
                {
                    arcStatPoint = b;
                    arcEndPoint = a;
            }
           
    }
}
Vector3d BigCircle::MX(Vector3d v, double u) {

    Matrix3d* a = new Matrix3d;
    *a << 1, 0, 0,
        0, cos(u), sin(u),
        0, -sin(u), cos(u);
    return *a * v;
}
Vector3d BigCircle::MY(Vector3d v, double u) {

    Matrix3d *a = new Matrix3d;
    *a << cos(u), 0, -sin(u),
        0, 1, 0,
        sin(u), 0, cos(u);
    return *a * v;
}
Vector3d BigCircle::MZ(Vector3d v, double u) {

    Matrix3d *a = new Matrix3d;
    *a << cos(u), sin(u), 0,
        -sin(u), cos(u), 0,
        0, 0, 1;
    return *a * v;
}
void BigCircle::init()
{
    penForBack = gcnew Pen(Color::Black);
    
    penForBack->Width = 1;
    penForBack->DashStyle = DashStyle::Dash;

    //penForBack->LineJoin = LineJoin::Round;

    normal = new Vector3d();
    zenit = new Vector3d();
    nadir = new Vector3d();

    for (int i = 0; i< COUNT_OF_ARR_POINTS; i++, alpha += s, tochki++) {

        if (i < 100) {
            dots[i] = new Vector3d();
        }

        coords[i] = new Vector3d();
        coords[i]->data()[0] = cos(alpha) * unit;
        coords[i]->data()[1] = -sin(alpha) * unit;
        coords[i]->data()[2] = 0;
        //points->Add(PointF(coords[i]->data()[0] + pW / 2.0, -coords[i]->data()[1] + pH / 2.0));
        //добавляем точку в коллекцию. полученные координаты сразу переводим в экранные единицы}


        axis[i] = new Vector3d();
        axis[i]->data()[0] = 0;
        axis[i]->data()[1] = 0;
        
        if (i == 0) {
            axis[i]->data()[2] = 0;
        }
        else if (i <= 500) {
            
        axis[i]->data()[2] = -unit + unit/(double)i;
        }
        else {
            axis[i]->data()[2] = unit - unit / (double)(i - 500);
        }

        pointsAxis->Add(PointF(axis[i]->data()[0] + pW / 2.0, -axis[i]->data()[1] + pH / 2.0));
    }
    
    (*normal) << 0, 0, 1;
    (*zenit) <<0, 0, (double)unit;
    *nadir << 0, 0, -(double)unit;
}
BigCircle::~BigCircle()
{
}

double BigCircle::cordFI(Vector3d v) {
    
    if (v.data()[1] > 0) {
        return acos(v.data()[0] / (sqrt(pow(v.data()[0], 2.0) + pow(v.data()[1], 2.0))));
    }
    else
    {
        return -acos(v.data()[0] / (sqrt(pow(v.data()[0], 2.0) + pow(v.data()[1], 2.0))));
    }

}

double BigCircle::cordTETA(Vector3d v) {
    return acos(v.data()[2]);
}
//Vector3d astropractica::BigCircle::plane(Vector3d A, Vector3d B)
//{
//    Vector3d *a = new Vector3d();
//
//    *a << A.data()[1] * B.data()[2] - A.data()[2] * B.data()[1], //коэф при х
//          A.data()[2] * B.data()[0] - A.data()[0] * B.data()[2], //коэф при y
//          A.data()[0] * B.data()[1] - A.data()[1] * B.data()[0];  //коэф при z
//    return *a;
//}


/////////////////////////////////////old//////////////////////////////////////////////////////////////
//Vector3d toScreen(Vector3d v) {
//
//    Matrix3d C;
//    C << sqrt(3), 0, -sqrt(3),
//        1, 2, 1,
//        sqrt(2), -sqrt(2), sqrt(2);
//
//    return (1 / sqrt(6)) * C * v;
//}
//Vector3d MX(Vector3d v, double u) {
//        
//    Matrix3d a;
//    a << 1,   0,      0,
//         0,  cos(u), sin(u),
//         0, -sin(u), cos(u);
//    return a * v;
//}
//Vector3d MY(Vector3d v, double u) {
//
//    Matrix3d a;
//    a << cos(u), 0, -sin(u),
//           0,    1,    0,
//         sin(u), 0, cos(u);
//    return a * v;
//}
//Vector3d MZ(Vector3d v, double u) {
//
//    Matrix3d a;
//    a <<  cos(u), sin(u), 0,
//         -sin(u), cos(u), 0,
//            0,      0,    1;
//    return a * v;
//}
//array<double>^ rotation(array<double>^ A) {
//
//    double b = acos(A[0] / sqrt(pow(A[0], 2) + pow(A[1], 2)));
//    double c = acos(A[2]/sqrt(pow(A[0], 2)+ pow(A[1], 2)+ pow(A[2], 2)));
//
//    array<double>^ a = { b, c };
//    array<double>^ d = { -b, c };
//    return (A[1] > 0) ? a : d;
//}
//array<double>^ plane(array<double>^ A, array<double>^ B) {
//
//    array<double>^ a = {
//        A[1]*B[2]-A[2]*B[1], //коэф при х
//        A[2]*B[0]-A[0]*B[2], //коэф при y
//        A[0]*B[1]-A[1]*B[0]  //коэф при z
//    };
//    return a;
//}
//double reneY(array<double>^ Array) {
//    Vector3d a = reneXYZ(Array[0], Array[1]);
//    return a[1];
//}
//double reneX(array<double>^ Array) {
//    Vector3d a = reneXYZ(Array[0], Array[1]);
//    return a[0];
//}
//double reneY(double x, double y) {
//    Vector3d a = reneXYZ(x, y);
//    return a[1];
//}
//double reneX(double x, double y) {
//    Vector3d a = reneXYZ(x, y);
//    return a[0];
//}
//double reneZ(double x, double y) {
//    Vector3d a = reneXYZ(x, y);
//    return a[2];
//}
//double reneZ(array<double>^ Array) {
//    Vector3d a = reneXYZ(Array[0], Array[1]);
//    return a[2];
//}
//Vector3d reneXYZ(double x, double y) {
//
//    Matrix3d C;
//    Vector3d A;
//    C << sqrt(3),   0,    -sqrt(3),
//            1,      2,        1,
//        sqrt(2), -sqrt(2), sqrt(2);
//
//    C.inverse();
//    A << x, y, 0; //TODO доделать, траблы с зет штрих
//
//    A = sqrt(6) * C * A;
//
//    return A;
//}
//double cordX(double fi, double teta) {
//
//    double a = sqrt(0.5) * (cos(fi) * sin(teta) - cos(teta));
//    return a;
//}
//double cordX(array<double>^ Array) {
//    return cordX(Array[0], Array[1]);
//}
//double cordY(double fi, double teta) {
//    return sqrt(1.0 / 6.0) * (cos(fi) * sin(teta) + 2 * sin(fi) * sin(teta) + cos(teta));
//}
//double cordY(array<double>^ Array) {
//    return cordY(Array[0], Array[1]);
//}
//array<double>^ cordXY(double fi, double teta) {
//    array<double>^ a = { cordX(fi, teta), cordY(fi, teta) };
//
//    return a;
//}
//array<double>^ cordXY(array<double>^ Array) {
//    array<double>^ a = { cordX(Array[0], Array[1]), cordY(Array[0], Array[1]) };
//
//    return a;
//}
//double cordFI(double x, double y) {
//    return y/x;
//}
//double cordTETA(double x, double y) {
//    return atan(sqrt((pow(x, 2.0) + pow(y, 2)) / (sqrt(1 - pow(x, 2.0) - pow(y, 2)))));
//}
//double cordFI(array<double>^ Array) {
//    return cordFI(Array[0], Array[1]);
//}
//double cordTETA(array<double>^ Array) {
//    return cordTETA(Array[0], Array[1]);
//}
//array<double>^ cordFI_TETA(double x, double y) {
//    array<double>^ a = { cordFI(x, y), cordTETA(x, y) };
//    return a;
//}
//array<double>^ cordFI_TETA(array<double>^ Array) {
//    array<double>^ a = { cordFI(Array[0], Array[1]), cordTETA(Array[0], Array[1]) };
//
//    return a;
//}