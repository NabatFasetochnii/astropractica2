#pragma once
#include "MyForm.h"
#include <Math.h>
#include <iostream>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <Math.h>
#include <string>

namespace astropractica {

	using namespace System::Drawing;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Collections::Generic;
	using namespace System;
	using namespace Eigen;
	using namespace System::Drawing::Drawing2D;
	using namespace System::Windows::Forms;

	/// <summary>
	/// Сводка для MyForm
	/// </summary>
	

	public ref class BigCircle //класс большого круга
	{
	public:
		int getCountOfPoint(); //получить количество точек отрисовки круга 
		BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis, int pW, int pH); //конструкторы
		BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, int pW, int pH);
		BigCircle(Graphics^ graphics, Pen^ pen, int pW, int pH);
		BigCircle(Graphics^ graphics, int pW, int pH);
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
		bool addDot(int x, int y); //добавить на круг точку 
		bool isDotOnCircle(Vector3d v); //проверить находится ли точка на круге
		bool isDotOnCircle(int x, int y); //проверить находится ли точка на круге
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
		const int CONT_OF_ARR_POINTS = 629;
		Vector3d MX(Vector3d v, double u); //оперрации поворота вокруг соотвествующей оси
		Vector3d MY(Vector3d v, double u); //"здесь поворот совершается против часовой стрелки
		Vector3d MZ(Vector3d v, double u); //а объект по сути поворачивается в другую сторону, по часовой стрелке" (с) Аналитик
		//на вход идёт вектор и угол в радианах, на который нужно повернуть
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Graphics^ graphics; // объект для отрисовки
		Vector3d* normal; //нормаль к поверхности
		array<Vector3d*, 1>^ coords = gcnew array<Vector3d*>(CONT_OF_ARR_POINTS); // массив координат точек круга //исправил
		array<Vector3d*, 1>^ axis = gcnew array<Vector3d*>(CONT_OF_ARR_POINTS); //массив координат точек оси //исправил
		array<Vector3d*, 1>^ dots = gcnew array<Vector3d*>(CONT_OF_ARR_POINTS); //массив координат пользовательских точек // исправил 
		List<PointF>^ points = gcnew List<PointF>(); //коллекция точек большого круга, для отрисовки
		List<PointF>^ pointsAxis = gcnew List<PointF>(); //коллекция точек оси, для отрисовки
		double s = 0.01; //шаг по углу построения точек круга 
		double alpha = 0; // угол поворота, нужен для построения круга
		int tochki = 0; // колличество точек круга
		int unit = 200; //массштаб
		Pen^ pen; //ручка для круга
		Pen^ penForDots; //ручка для точек
		Pen^ penForAxis; //ручка для оси
		Vector3d* zenit; //точка пересечения оси с небсферой 1
		Vector3d* nadir;//точка пересечения оси с небсферой 2
		int countOfDots = 0; //счётчик пользовательских точек
		void init(); // построение круга
		const double EPS = 8; //бесконечно малое эпсилон
		int pW, pH;
		Vector3d *buf = new Vector3d;
	};

	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		int pW;
		int pH;
		Bitmap^ img;
		Point p;
		int x, y;
		BigCircle^ big;
		String^ buf;
		bool isPicStart = false;

		MyForm(void)
		{
			
			InitializeComponent();
			//
			//TODO: добавьте код конструктора
			//
			pW = pic->Width;
			pH = pic->Height;
			img = gcnew Bitmap(pW, pH);
		}

	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}

	
	private: PictureBox^ pic;
	private: System::Windows::Forms::Button^ button1;
	
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
	
		void InitializeComponent(void)
		{
			this->pic = (gcnew System::Windows::Forms::PictureBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pic))->BeginInit();
			this->SuspendLayout();
			// 
			// pic
			// 
			this->pic->Location = System::Drawing::Point(12, 12);
			this->pic->Name = L"pic";
			this->pic->Size = System::Drawing::Size(819, 571);
			this->pic->TabIndex = 0;
			this->pic->TabStop = false;
			this->pic->Click += gcnew System::EventHandler(this, &MyForm::pictureBox1_Click);
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(847, 113);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(255, 204);
			this->button1->TabIndex = 1;
			this->button1->Text = L"button1";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1114, 718);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->pic);
			this->Name = L"MyForm";
			this->Text = L"astropractica";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pic))->EndInit();
			this->ResumeLayout(false);

		}
#pragma endregion
	private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
	}
	private: System::Void pictureBox1_Click(System::Object^ sender, System::EventArgs^ e) {
			
		if (isPicStart) {
		p = pic->PointToClient(Cursor->Position);
		big->addDot(p.X, p.Y);
		big->drowDots();
		this->pic->Image = img;
		}
		
	}
	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {

		isPicStart = true;
		Graphics^ g = Graphics::FromImage(img);
		big = gcnew BigCircle(g, pW, pH); //TODO защита от идиота
		big->rotationY(0.5);
		big->rotationX(0.5);
		big->rotationZ(1);
		//big->drowCircle();
		big->onDrowAll();
		/*big->drowAxis();
		big->drowDots();
		big->drowZenit();
		big->drowNadir();*/
		this->pic->Image = img;
	}
	
	};
	
	

}
