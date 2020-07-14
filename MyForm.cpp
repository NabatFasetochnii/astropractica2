#include "MyForm.h"
#include <Math.h>
#include <iostream>
#include <Eigen/Dense>

using namespace System;
using namespace Eigen;

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
array<double>^ reneXY(double x, double y);
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


int main() {

    double x = 10, y = 40.5;

    Console::WriteLine(reneX(x, y));
    Console::WriteLine(reneY(x, y));

    /*double fi = 30;
    double teta = 60;

    Console::WriteLine("cordX1: ");
    Console::WriteLine(cordX(fi, teta));

    array<double>^ Array =  {fi, teta};

    Console::WriteLine("cordX2: ");
    Console::WriteLine(cordX(Array));

    Console::WriteLine("cordY1: ");
    Console::WriteLine(cordY(fi, teta));

    Console::WriteLine("cordY2: ");
    Console::WriteLine(cordY(Array));

    Console::WriteLine("cordXY1: ");
    Console::WriteLine(cordXY(fi, teta)[0]);
    Console::WriteLine(cordXY(fi, teta)[1]);

    Console::WriteLine("cordXY2: ");
    Console::WriteLine(cordXY(Array)[0]);
    Console::WriteLine(cordXY(Array)[1]);*/

    Console::ReadKey();

    return 0;
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
    array<double>^ a = reneXY(Array[0], Array[1]);
    return a[1];
}
double reneX(array<double>^ Array) {
    array<double>^ a = reneXY(Array[0], Array[1]);
    return a[0];
}
double reneY(double x, double y) {
    array<double>^ a = reneXY(x, y);
    return a[1];
}
double reneX(double x, double y) {
    array<double>^ a = reneXY(x, y);
    return a[0];
}
array<double>^ reneXY(double x, double y) {

    Matrix3d C;
    Vector3d A;
    C << sqrt(3), 0, -sqrt(3),
        1, 2, 1,
        sqrt(2), -sqrt(2), sqrt(2);

    C.inverse();
    A << x, y, 0;

    A = sqrt(6) * C * A;

    array<double>^ a = { (A[0]), (A[1]) };
    return a;
   
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