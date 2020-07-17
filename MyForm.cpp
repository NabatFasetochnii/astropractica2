#include "MyForm.h"
#include <Math.h>
#include <iostream>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <Math.h>


using namespace System::Drawing;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Collections::Generic;
using namespace System;
using namespace Eigen;
using namespace System::Drawing::Drawing2D;

//формулы перехода от точкин на поверхности сферы к точке на экране
double cordX(double fi, double teta);
double cordX(array<double>^ Array); // в массиве на вход первый элемент фи, второй - тета
double cordY(double fi, double teta);
double cordY(array<double>^ Array);
array<double>^ cordXY(double fi, double teta); // в массиве на выход первый элемент х, второй - у
array<double>^ cordXY(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//координаты сферические по кооридинатам декартовым (решения у этой штуки должно быть два) - слова аналитика
double cordFI(double x, double y); //на вход подаются декартовы координаты точки    
double cordFI(array<double>^ Array); 
double cordTETA(double x, double y);
double cordTETA(array<double>^ Array);
array<double>^ cordFI_TETA(double x, double y);
array<double>^ cordFI_TETA(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//координаты точки в декартовых через координаты на экране
double reneX(double x, double y); // на вход подаются координаты на экране, на выход декартовы координаты точки на сфере
double reneX(array<double>^ Array);
double reneY(double x, double y);
double reneY(array<double>^ Array);
Vector3d reneXYZ(double x, double y);
double reneZ(double x, double y);
double reneZ(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Коэффициенты для уравнения плоскости по двум точкам (+ начало координат) в декартовой
array<double>^ plane(array<double>^ A, array<double>^ B); //На вход подаются два массива,
//которые характеризуют точки большого круга, то есть содеоржат координаты x, y,z
//на выход подаётся массив с коэффициентами
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//углы поворота большого круга (фи и тета)
array<double>^ rotation(array<double>^ A); //
//координаты нормали N - на входе
//углы поворота(фи - вокруг оси зед, тета - вокруг оси игрек) для того чтобы повернуть нормаль(0.0.1) 
//так чтобы она стала параллельно нормали N - на выходе
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Vector3d MX(Vector3d v, double u); //оперрации поворота вокруг соотвествующей оси
Vector3d MY(Vector3d v, double u); //"здесь поворот совершается против часовой стрелки
Vector3d MZ(Vector3d v, double u); //а объект по сути поворачивается в другую сторону, по часовой стрелке"
//на вход идёт вектор и угол, на который нужно повернуть
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Vector3d toScreen(Vector3d v);

int main() {

    array<Vector2cd*>^ a = gcnew array<Vector2cd*>(20);

    Vector3d v, r; //объявляем трёх мерный вектор, который имеет значения координат типа дабл
    v << 1, 1, 1; //заполняем его единицами

    //надо повернуть на -0.785398 радиан по Y и на 0.61540309 радиан по X

    r = MY(v, -0.785398); //поворачиваем по Y точку с координатами (1,1,1)
    r = MX(r, 0.61540309); //поворачиваем по X

    Vector3d f = toScreen(r); //используя функцию получения координат экрана по координатам декарта, инициализируем новый вектор
    Console::WriteLine(f.data()[0]); //выводим координаты точки
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::WriteLine("\n"); // выводим пустую строку для удобства восприятия 

    r = MY(v*2, -0.785398); //поворачиваем по Y точку с координатами (2,2,2)
    r = MX(r, 0.61540309); //поворачиваем по X точку с неизвестными координатами

    f = toScreen(r); // ищем координаты на экране
    Console::WriteLine(f.data()[0]); //выводим координаты точки
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::WriteLine("\n");

    r = MY(v * 3, -0.785398); //поворачиваем по Y точку с координатами (3,3,3)
    r = MX(r, 0.61540309); //поворачиваем по X точку с неизвестными координатами

    f = toScreen(r);// ищем координаты на экране
    Console::WriteLine(f.data()[0]);//выводим координаты точки
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::ReadKey();

    return 0;
}

ref class BigCircle //класс большого круга
{
public:
    int getCountOfPoint(); //получить количество точек отрисовки круга 
    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis); //конструкторы
    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots);
    BigCircle(Graphics^ g, Pen^ pen);
    BigCircle(Graphics^ graphics);
    ~BigCircle(); //дискриптор 
    void rotationX(double u); //ворочаем по осям на угол u (в радианах)
    void rotationY(double u);
    void rotationZ(double u);
    void onDrowAll(); //нарисовать всё
    array<Vector3d*, 1>^ getCircleCoords(); //получить массив векторов, соответствующих точкам круга
    array<Vector3d*, 1>^ getAxis(); //получить массив векторов, соответствующих точкам оси круга
    Vector3d getNormal(); //получить вектор нормали 
    Vector3d getZenit();  //получить точку пересечения оси с небесной сферой 1
    Vector3d getNadir();  //получить точку пересечения оси с небесной сферой 2
    bool addDot(Vector3d v); //добавить на круг точку 
    bool isDotOnCircle(Vector3d v); //проверить находится ли точка на круге
    int getCountOfDots(); //получить счётчик пользовательских точек 
    array<Vector3d*, 1>^ getDots(); //получить массив пользовательских точек 
    void drowCircle(); // нарисовать по отдельности круг
    void drowAxis(); // ось
    void drowDots();// точки
    void drowZenit(); //точку пересечения 1
    void drowNadir(); //точку пересечения 2
    double dihedralAngle_rad(BigCircle b); //двугранный угол между двумя плоскостями больших кругов, в радианах
    double dihedralAngle_deg(BigCircle b);//двугранный угол между двумя плоскостями больших кругов, в градусах
    Vector3d intersectionDots(BigCircle b); //получить точку пересечения двух больших кругов, вторая точка является противоположной этой, то есть просто можно домножить на -1

private:

    Vector3d MX(Vector3d v, double u); //оперрации поворота вокруг соотвествующей оси
    Vector3d MY(Vector3d v, double u); //"здесь поворот совершается против часовой стрелки
    Vector3d MZ(Vector3d v, double u); //а объект по сути поворачивается в другую сторону, по часовой стрелке" (с) Аналитик
    //на вход идёт вектор и угол в радианах, на который нужно повернуть
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Graphics^ graphics; // объект для отрисовки
    Vector3d* normal; //нормаль к поверхности
    array<Vector3d*, 1>^ coords = gcnew array<Vector3d*>(700); // массив координат точек круга 
    array<Vector3d*, 1>^ axis = gcnew array<Vector3d*>(700); //массив координат точек оси
    array<Vector3d*, 1>^ dots = gcnew array<Vector3d*>(100); //массив координат пользовательских точек
    List<PointF>^ points = gcnew List<PointF>(); //коллекция точек большого круга, для отрисовки
    List<PointF>^ pointsAxis = gcnew List<PointF>(); //коллекция точек оси, для отрисовки
    double s = 0.01; //шаг по углу построения точек круга 
    double alpha = 0; // угол поворота, нужен для построения круга
    int tochki = 0; // колличество точек круга
    int unit = 25; //массштаб
    Pen^ pen; //ручка для круга
    Pen^ penForDots; //ручка для точек
    Pen^ penForAxis; //ручка для оси
    Vector3d* zenit; //точка пересечения оси с небсферой 1
    Vector3d* nadir;//точка пересечения оси с небсферой 2
    int countOfDots = 0; //счётчик пользовательских точек
    void init(); // построение круга
    const double EPS = 0.1; //бесконечно малое эпсилон

};
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
bool BigCircle::addDot(Vector3d v)
{
    if (isDotOnCircle(v)) {
        *dots[countOfDots] = v;
        countOfDots++;
        return true;
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
    graphics->DrawLines(pen, points->ToArray());

}
void BigCircle::drowAxis()
{
    graphics->DrawLines(pen, pointsAxis->ToArray());
}
void BigCircle::drowDots()
{
    if (countOfDots > 0) {
        for each (Vector3d * v in dots) {
                graphics->DrawEllipse(penForDots, v->data()[0], v->data()[1], 5, 5);
        }
    }
    
}
void BigCircle::drowZenit()
{
    graphics->DrawEllipse(penForDots, zenit->data()[0], zenit->data()[1], 5, 5);
}
void BigCircle::drowNadir()
{
    graphics->DrawEllipse(penForDots, nadir->data()[0], nadir->data()[1], 5, 5);
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
void BigCircle::rotationX(double u) {
    for each (Vector3d* v in coords)
    {
       
        (*v) = MX(*v, u); 
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    for each (Vector3d * v in axis)
    {
        (*v) = MX(*v, u);
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    *normal = MX(*normal, u);
    *zenit = (*normal * (double)unit);
    *nadir = -(*normal * (double)unit);
}
void BigCircle::rotationY(double u) {
    for each (Vector3d * v in coords)
    {
        (*v) = MY(*v, u);
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    for each (Vector3d * v in axis)
    {
        (*v) = MY(*v, u);
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    *normal = MY(*normal, u);
    *zenit = (*normal * (double)unit);
    *nadir = -(*normal * (double)unit);
}
void BigCircle::rotationZ(double u) {
    for each (Vector3d * v in coords)
    {
        (*v) = MZ(*v, u);
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    for each (Vector3d * v in axis)
    {
        (*v) = MZ(*v, u);
        points->Clear();
        points->Add(PointF(v->data()[0], v->data()[1]));
    }
    *normal = MZ(*normal, u);
    *zenit = (*normal * (double)unit);
    *nadir = -(*normal * (double)unit);
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis)
{
    this->graphics = graphics;
    this->pen = pen;
    this->penForDots = penForDots;

    init();
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForAxis)
{
    this->graphics = graphics;
    this->pen = pen;
    this->penForAxis = penForAxis;

    penForDots = gcnew Pen(Color::Red);
    penForDots->Width = 10;
    penForDots->LineJoin = LineJoin::Round;

    init();
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots)
{
    this->graphics = graphics;
    this->pen = pen;
    this->penForDots = penForDots;

    penForAxis = gcnew Pen(Color::Gray);
    penForAxis->Width = 5;
    penForAxis->LineJoin = LineJoin::Round;
    
    init();
}
BigCircle::BigCircle(Graphics^ graphics, Pen^ pen)
{
    this->graphics = graphics;
    this->pen = pen;

    penForDots = gcnew Pen(Color::Red);
    penForDots->Width = 10;
    penForDots->LineJoin = LineJoin::Round;

    penForAxis = gcnew Pen(Color::Gray);
    penForDots->Width = 5;
    penForDots->LineJoin = LineJoin::Round;
    init();
}
BigCircle::BigCircle(Graphics^ graphics)
{
    this->graphics = graphics;

    pen = gcnew Pen(Color::Black);
    pen->Width = 8;
    pen->LineJoin = LineJoin::Round;

    penForDots = gcnew Pen(Color::Red);
    penForDots->Width = 10;
    penForDots->LineJoin = LineJoin::Round;

    penForAxis = gcnew Pen(Color::Gray);
    penForDots->Width = 5;
    penForDots->LineJoin = LineJoin::Round;

    init();
}
Vector3d BigCircle::MX(Vector3d v, double u) {

    Matrix3d a;
    a << 1, 0, 0,
        0, cos(u), sin(u),
        0, -sin(u), cos(u);
    return a * v;
}
Vector3d BigCircle::MY(Vector3d v, double u) {

    Matrix3d a;
    a << cos(u), 0, -sin(u),
        0, 1, 0,
        sin(u), 0, cos(u);
    return a * v;
}
Vector3d BigCircle::MZ(Vector3d v, double u) {

    Matrix3d a;
    a << cos(u), sin(u), 0,
        -sin(u), cos(u), 0,
        0, 0, 1;
    return a * v;
}
void BigCircle::init()
{
    for (int i = 0; alpha < 2 * M_PI; i++, alpha += s, tochki++) {

        coords[i]->data()[0] = cos(alpha) * unit;
        coords[i]->data()[1] = sin(alpha) * unit;
        coords[i]->data()[2] = 0;
        points->Add(PointF(coords[i]->data()[0], coords[i]->data()[1]));
        //добавляем точку в коллекцию. полученные координаты сразу переводим в экранные единицы}

        axis[i]->data()[0] = 0;
        axis[i]->data()[1] = 0;
        axis[i]->data()[2] = pow(-1.0, i) * i / 2;
        pointsAxis->Add(PointF(axis[i]->data()[0], axis[i]->data()[1]));
    }
    (*normal) << 0, 0, 1;
    *zenit = (*normal * (double)unit);
    *nadir = -(*normal * (double)unit);
}
BigCircle::~BigCircle()
{
}

Vector3d toScreen(Vector3d v) {

    Matrix3d C;
    C << sqrt(3), 0, -sqrt(3),
        1, 2, 1,
        sqrt(2), -sqrt(2), sqrt(2);

    return (1 / sqrt(6)) * C * v;
}
Vector3d MX(Vector3d v, double u) {
        
    Matrix3d a;
    a << 1,   0,      0,
         0,  cos(u), sin(u),
         0, -sin(u), cos(u);
    return a * v;
}
Vector3d MY(Vector3d v, double u) {

    Matrix3d a;
    a << cos(u), 0, -sin(u),
           0,    1,    0,
         sin(u), 0, cos(u);
    return a * v;
}
Vector3d MZ(Vector3d v, double u) {

    Matrix3d a;
    a <<  cos(u), sin(u), 0,
         -sin(u), cos(u), 0,
            0,      0,    1;
    return a * v;
}
array<double>^ rotation(array<double>^ A) {

    double b = acos(A[0] / sqrt(pow(A[0], 2) + pow(A[1], 2)));
    double c = acos(A[2]/sqrt(pow(A[0], 2)+ pow(A[1], 2)+ pow(A[2], 2)));

    array<double>^ a = { b, c };
    array<double>^ d = { -b, c };
    return (A[1] > 0) ? a : d;
}
array<double>^ plane(array<double>^ A, array<double>^ B) {

    array<double>^ a = {
        A[1]*B[2]-A[2]*B[1], //коэф при х
        A[2]*B[0]-A[0]*B[2], //коэф при y
        A[0]*B[1]-A[1]*B[0]  //коэф при z
    };
    return a;
}
double reneY(array<double>^ Array) {
    Vector3d a = reneXYZ(Array[0], Array[1]);
    return a[1];
}
double reneX(array<double>^ Array) {
    Vector3d a = reneXYZ(Array[0], Array[1]);
    return a[0];
}
double reneY(double x, double y) {
    Vector3d a = reneXYZ(x, y);
    return a[1];
}
double reneX(double x, double y) {
    Vector3d a = reneXYZ(x, y);
    return a[0];
}
double reneZ(double x, double y) {
    Vector3d a = reneXYZ(x, y);
    return a[2];
}
double reneZ(array<double>^ Array) {
    Vector3d a = reneXYZ(Array[0], Array[1]);
    return a[2];
}
Vector3d reneXYZ(double x, double y) {

    Matrix3d C;
    Vector3d A;
    C << sqrt(3),   0,    -sqrt(3),
            1,      2,        1,
        sqrt(2), -sqrt(2), sqrt(2);

    C.inverse();
    A << x, y, 0; //TODO доделать, траблы с зет штрих

    A = sqrt(6) * C * A;

    return A;
}
double cordX(double fi, double teta) {

    double a = sqrt(0.5) * (cos(fi) * sin(teta) - cos(teta));
    return a;
}
double cordX(array<double>^ Array) {
    return cordX(Array[0], Array[1]);
}
double cordY(double fi, double teta) {
    return sqrt(1.0 / 6.0) * (cos(fi) * sin(teta) + 2 * sin(fi) * sin(teta) + cos(teta));
}
double cordY(array<double>^ Array) {
    return cordY(Array[0], Array[1]);
}
array<double>^ cordXY(double fi, double teta) {
    array<double>^ a = { cordX(fi, teta), cordY(fi, teta) };

    return a;
}
array<double>^ cordXY(array<double>^ Array) {
    array<double>^ a = { cordX(Array[0], Array[1]), cordY(Array[0], Array[1]) };

    return a;
}
double cordFI(double x, double y) {
    return y/x;
}
double cordTETA(double x, double y) {
    return atan(sqrt((pow(x, 2.0) + pow(y, 2)) / (sqrt(1 - pow(x, 2.0) - pow(y, 2)))));
}
double cordFI(array<double>^ Array) {
    return cordFI(Array[0], Array[1]);
}
double cordTETA(array<double>^ Array) {
    return cordTETA(Array[0], Array[1]);
}
array<double>^ cordFI_TETA(double x, double y) {
    array<double>^ a = { cordFI(x, y), cordTETA(x, y) };
    return a;
}
array<double>^ cordFI_TETA(array<double>^ Array) {
    array<double>^ a = { cordFI(Array[0], Array[1]), cordTETA(Array[0], Array[1]) };

    return a;
}