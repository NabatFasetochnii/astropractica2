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
	/// ������ ��� MyForm
	/// </summary>
	

	public ref class BigCircle //����� �������� �����
	{
	public:
		int getCountOfPoint(); //�������� ���������� ����� ��������� ����� 
		BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, Pen^ penForAxis, int pW, int pH); //������������
		BigCircle(Graphics^ graphics, Pen^ pen, Pen^ penForDots, int pW, int pH);
		BigCircle(Graphics^ graphics, Pen^ pen, int pW, int pH);
		BigCircle(Graphics^ graphics, int pW, int pH);
		BigCircle(Graphics^ graphics, int pW, int pH, int x1, int y1, int x2, int y2);
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
		bool addDot(int x, int y); //�������� �� ���� �����. �� ��������
		bool isDotOnCircle(Vector3d v); //��������� ��������� �� ����� �� �����
		bool isDotOnCircle(int x, int y); //��������� ��������� �� ����� �� �����
		int getCountOfDots(); //�������� ������� ���������������� ����� 
		array<Vector3d*, 1>^ getDots(); //�������� ������ ���������������� ����� 
		void drowCircle(); // ���������� �� ����������� ����
		void drowAxis(); // ���
		void drowDots();// �����
		void drowZenit(); //����� ����������� 1
		void drowNadir(); //����� ����������� 2
		void drowArc();//�� �������� 
		double dihedralAngle_rad(BigCircle b); //���������� ���� ����� ����� ����������� ������� ������, � ��������
		double dihedralAngle_deg(BigCircle b);//���������� ���� ����� ����� ����������� ������� ������, � ��������
		Vector3d intersectionDots(BigCircle b); //�������� ����� ����������� ���� ������� ������, ������ ����� �������� ��������������� ����, �� ���� ������ ����� ��������� �� -1
		void metod1(int X1,int Y1, int X2, int Y2);
		void setArc(int X1, int Y1, int X2, int Y2); //�� ��������
	private:
		const int COUNT_OF_ARR_POINTS = 629;
		Vector3d MX(Vector3d v, double u); //��������� �������� ������ �������������� ���
		Vector3d MY(Vector3d v, double u); //"����� ������� ����������� ������ ������� �������
		Vector3d MZ(Vector3d v, double u); //� ������ �� ���� �������������� � ������ �������, �� ������� �������" (�) ��������
		//�� ���� ��� ������ � ���� � ��������, �� ������� ����� ���������
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Graphics^ graphics; // ������ ��� ���������
		Vector3d* normal; //������� � �����������
		array<Vector3d*, 1>^ coords = gcnew array<Vector3d*>(COUNT_OF_ARR_POINTS); // ������ ��������� ����� ����� //��������
		array<Vector3d*, 1>^ axis = gcnew array<Vector3d*>(COUNT_OF_ARR_POINTS); //������ ��������� ����� ��� //��������
		array<Vector3d*, 1>^ dots = gcnew array<Vector3d*>(COUNT_OF_ARR_POINTS); //������ ��������� ���������������� ����� // �������� 
		//List<PointF>^ points = gcnew List<PointF>(); //��������� ����� �������� �����, ��� ���������
		List<PointF>^ pointsAxis = gcnew List<PointF>(); //��������� ����� ���, ��� ���������
		double s = 0.01; //��� �� ���� ���������� ����� ����� 
		double alpha = 0; // ���� ��������, ����� ��� ���������� �����
		int tochki = 0; // ����������� ����� �����
		int unit = 200; //�������
		Pen^ pen; //����� ��� �����
		Pen^ penForDots; //����� ��� �����
		Pen^ penForAxis; //����� ��� ���
		Vector3d* zenit; //����� ����������� ��� � ��������� 1
		Vector3d* nadir;//����� ����������� ��� � ��������� 2
		int countOfDots = 0; //������� ���������������� �����
		void init(); // ���������� �����
		const double EPS = 8; //���������� ����� �������
		const double EPS2 = 1; //���������� ����� �������2
		int pW, pH; // ������ �����������
		Vector3d *buf = new Vector3d;
		double cordFI(Vector3d v); //������� �� oz
		double cordTETA(Vector3d v); //������� �� oy
		//Vector3d plane(Vector3d v, Vector3d b); // ���� �����������
		int arcStatPoint;
		int arcEndPoint;
		bool isArcSet = false;
		Pen^ penForBack;

	};

	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		int pW;
		int pH;
		Bitmap^ img;
		Graphics^ g;
		Point p;
		int x, y;
		BigCircle^ big;
		BigCircle^ big2;
		String^ buf;
		bool isPicStart = false;
		bool firstClick = true;
		int x1;
		int y1;
		bool b = true;

		MyForm(void)
		{
			
			InitializeComponent();
			//
			//TODO: �������� ��� ������������
			//
			pW = pic->Width;
			pH = pic->Height;
			img = gcnew Bitmap(pW, pH);
		}

	protected:
		/// <summary>
		/// ���������� ��� ������������ �������.
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
			//big->rotationX(0.1);
			if (firstClick) {
			
				x1 = p.X;
				y1 = p.Y;
				firstClick = false;
			}
			else
			{
				//g->Clear(Color::White);
				/*big2 = gcnew BigCircle(g, pW, pH, x1, y1, p.X, p.Y);
				big2->onDrowAll();*/
				b = false;
				big->setArc(x1, y1, p.X, p.Y);
				big->drowArc();
				
				x1 = p.X;
				y1 = p.Y;
			}

		this->pic->Image = img;
		}
		
	}
	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {

		isPicStart = true;
		g = Graphics::FromImage(img);
		big = gcnew BigCircle(g, pW, pH); 
		big->rotationY(1);
		//big->rotationX(0.5);
		//big->rotationZ(1);
		//big->drowCircle();
		if (b) {
			big->onDrowAll();
			/*big->drowAxis();
			big->drowDots();
			big->drowZenit();
			big->drowNadir();*/
			this->pic->Image = img;
		}
		else
		{
			//big->rotationY(0.1);
			//big->onDrowAll();
		}

		
	}
	
	};
	
	

}
