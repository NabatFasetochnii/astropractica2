#include "MyForm.h"
#include <Math.h>

using namespace System;

//формулы перехода от точкин на поверхности сферы к точке на экране
double cordX(double fi, double teta);
double cordX(array<double>^ Array); // в массиве на вход первый элемент фи, второй - тета
double cordY(double fi, double teta);
double cordY(array<double>^ Array);
array<double>^ cordXY(double fi, double teta); // в массиве на выход первый элемент х, второй - у
array<double>^ cordXY(array<double>^ Array);
//////////////////////////////////////////////////////////////////////
//координаты сферические по кооридинатам на экране (решения у этой штуки должно быть два) - слова аналитика
double cordFI(double x, double y);
double cordFI(array<double>^ Array); //аналогичная логика
double cordTETA(double x, double y);
double cordTETA(array<double>^ Array);
array<double>^ cordFI_TETA(double x, double y);
array<double>^ cordFI_TETA(array<double>^ Array);
//////////////////////////////////////////////////////////////////////


int main() {

    double fi = 30;
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
    Console::WriteLine(cordXY(Array)[1]);

    Console::ReadKey();

    return 0;
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