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