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

//������� �������� �� ������ �� ����������� ����� � ����� �� ������
double cordX(double fi, double teta);
double cordX(array<double>^ Array); // � ������� �� ���� ������ ������� ��, ������ - ����
double cordY(double fi, double teta);
double cordY(array<double>^ Array);
array<double>^ cordXY(double fi, double teta); // � ������� �� ����� ������ ������� �, ������ - �
array<double>^ cordXY(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//���������� ����������� �� ������������ ���������� (������� � ���� ����� ������ ���� ���) - ����� ���������
double cordFI(double x, double y); //�� ���� �������� ��������� ���������� �����    
double cordFI(array<double>^ Array); 
double cordTETA(double x, double y);
double cordTETA(array<double>^ Array);
array<double>^ cordFI_TETA(double x, double y);
array<double>^ cordFI_TETA(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//���������� ����� � ���������� ����� ���������� �� ������
double reneX(double x, double y); // �� ���� �������� ���������� �� ������, �� ����� ��������� ���������� ����� �� �����
double reneX(array<double>^ Array);
double reneY(double x, double y);
double reneY(array<double>^ Array);
Vector3d reneXYZ(double x, double y);
double reneZ(double x, double y);
double reneZ(array<double>^ Array);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//������������ ��� ��������� ��������� �� ���� ������ (+ ������ ���������) � ����������
array<double>^ plane(array<double>^ A, array<double>^ B); //�� ���� �������� ��� �������,
//������� ������������� ����� �������� �����, �� ���� ��������� ���������� x, y,z
//�� ����� ������� ������ � ��������������
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//���� �������� �������� ����� (�� � ����)
array<double>^ rotation(array<double>^ A); //
//���������� ������� N - �� �����
//���� ��������(�� - ������ ��� ���, ���� - ������ ��� �����) ��� ���� ����� ��������� �������(0.0.1) 
//��� ����� ��� ����� ����������� ������� N - �� ������
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Vector3d MX(Vector3d v, double u); //��������� �������� ������ �������������� ���
Vector3d MY(Vector3d v, double u); //"����� ������� ����������� ������ ������� �������
Vector3d MZ(Vector3d v, double u); //� ������ �� ���� �������������� � ������ �������, �� ������� �������"
//�� ���� ��� ������ � ����, �� ������� ����� ���������
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Vector3d toScreen(Vector3d v);

int main() {

    array<Vector2cd*>^ a = gcnew array<Vector2cd*>(20);

    Vector3d v, r; //��������� ��� ������ ������, ������� ����� �������� ��������� ���� ����
    v << 1, 1, 1; //��������� ��� ���������

    //���� ��������� �� -0.785398 ������ �� Y � �� 0.61540309 ������ �� X

    r = MY(v, -0.785398); //������������ �� Y ����� � ������������ (1,1,1)
    r = MX(r, 0.61540309); //������������ �� X

    Vector3d f = toScreen(r); //��������� ������� ��������� ��������� ������ �� ����������� �������, �������������� ����� ������
    Console::WriteLine(f.data()[0]); //������� ���������� �����
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::WriteLine("\n"); // ������� ������ ������ ��� �������� ���������� 

    r = MY(v*2, -0.785398); //������������ �� Y ����� � ������������ (2,2,2)
    r = MX(r, 0.61540309); //������������ �� X ����� � ������������ ������������

    f = toScreen(r); // ���� ���������� �� ������
    Console::WriteLine(f.data()[0]); //������� ���������� �����
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::WriteLine("\n");

    r = MY(v * 3, -0.785398); //������������ �� Y ����� � ������������ (3,3,3)
    r = MX(r, 0.61540309); //������������ �� X ����� � ������������ ������������

    f = toScreen(r);// ���� ���������� �� ������
    Console::WriteLine(f.data()[0]);//������� ���������� �����
    Console::WriteLine(f.data()[1]);
    Console::WriteLine(f.data()[2]);

    Console::ReadKey();

    return 0;
}

ref class BigCircle //����� �������� �����
{
public:
    int getCountOfPoint(); //�������� ���������� ����� ��������� ����� 
    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis); //������������
    BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots);
    BigCircle(Graphics^ g, Pen^ pen);
    BigCircle(Graphics^ graphics);
    ~BigCircle(); //���������� 
    void rotationX(double u); //�������� �� ���� �� ���� u (� ��������)
    void rotationY(double u);
    void rotationZ(double u);
    void onDrowAll(); //���������� ��
    array<Vector3d*, 1>^ getCircleCoords(); //�������� ������ ��������, ��������������� ������ �����
    array<Vector3d*, 1>^ getAxis(); //�������� ������ ��������, ��������������� ������ ��� �����
    Vector3d getNormal(); //�������� ������ ������� 
    Vector3d getZenit();  //�������� ����� ����������� ��� � �������� ������ 1
    Vector3d getNadir();  //�������� ����� ����������� ��� � �������� ������ 2
    bool addDot(Vector3d v); //�������� �� ���� ����� 
    bool isDotOnCircle(Vector3d v); //��������� ��������� �� ����� �� �����
    int getCountOfDots(); //�������� ������� ���������������� ����� 
    array<Vector3d*, 1>^ getDots(); //�������� ������ ���������������� ����� 
    void drowCircle(); // ���������� �� ����������� ����
    void drowAxis(); // ���
    void drowDots();// �����
    void drowZenit(); //����� ����������� 1
    void drowNadir(); //����� ����������� 2
    double dihedralAngle_rad(BigCircle b); //���������� ���� ����� ����� ����������� ������� ������, � ��������
    double dihedralAngle_deg(BigCircle b);//���������� ���� ����� ����� ����������� ������� ������, � ��������
    Vector3d intersectionDots(BigCircle b); //�������� ����� ����������� ���� ������� ������, ������ ����� �������� ��������������� ����, �� ���� ������ ����� ��������� �� -1

private:

    Vector3d MX(Vector3d v, double u); //��������� �������� ������ �������������� ���
    Vector3d MY(Vector3d v, double u); //"����� ������� ����������� ������ ������� �������
    Vector3d MZ(Vector3d v, double u); //� ������ �� ���� �������������� � ������ �������, �� ������� �������" (�) ��������
    //�� ���� ��� ������ � ���� � ��������, �� ������� ����� ���������
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Graphics^ graphics; // ������ ��� ���������
    Vector3d* normal; //������� � �����������
    array<Vector3d*, 1>^ coords = gcnew array<Vector3d*>(700); // ������ ��������� ����� ����� 
    array<Vector3d*, 1>^ axis = gcnew array<Vector3d*>(700); //������ ��������� ����� ���
    array<Vector3d*, 1>^ dots = gcnew array<Vector3d*>(100); //������ ��������� ���������������� �����
    List<PointF>^ points = gcnew List<PointF>(); //��������� ����� �������� �����, ��� ���������
    List<PointF>^ pointsAxis = gcnew List<PointF>(); //��������� ����� ���, ��� ���������
    double s = 0.01; //��� �� ���� ���������� ����� ����� 
    double alpha = 0; // ���� ��������, ����� ��� ���������� �����
    int tochki = 0; // ����������� ����� �����
    int unit = 25; //��������
    Pen^ pen; //����� ��� �����
    Pen^ penForDots; //����� ��� �����
    Pen^ penForAxis; //����� ��� ���
    Vector3d* zenit; //����� ����������� ��� � ��������� 1
    Vector3d* nadir;//����� ����������� ��� � ��������� 2
    int countOfDots = 0; //������� ���������������� �����
    void init(); // ���������� �����
    const double EPS = 0.1; //���������� ����� �������

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
        //��������� ����� � ���������. ���������� ���������� ����� ��������� � �������� �������}

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
        A[1]*B[2]-A[2]*B[1], //���� ��� �
        A[2]*B[0]-A[0]*B[2], //���� ��� y
        A[0]*B[1]-A[1]*B[0]  //���� ��� z
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
    A << x, y, 0; //TODO ��������, ������ � ��� �����

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