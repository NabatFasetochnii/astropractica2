#include "MyForm.h"
#include <Math.h>
#include <iostream>
#include <Eigen/Dense>

using namespace System;
using namespace Eigen;

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