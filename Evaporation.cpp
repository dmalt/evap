#include "iostream"
#include "Evaporation.h"
#include "math.h"

using namespace std;

enum{H2,O2,N2,H2O,OH,H,O,HO2,H2O2,	EC};     //


Evaporation::Evaporation(double * ExtConc, double ExtTemp, double ExtPressure){
	A[0]	=	-1649810.;
	A[1]	=	79239.;	
	A[2]	=	-1252.;
	A[3]	=	6.5;
// ****************Число компонент*******************
	CN=9;							
// ********Критическа температура для кислорода******
	T_c = 154.77; 						
// *********Выделение памяти под массивы*************
	Y_ex	= new double[CN];			
	Cp 		= new double[CN]; 			 
	Y_w		= new double[CN];			
// *****Инициализация внешних концентраций***********
	for (int i = 0; i < CN; ++i)
		Y_ex[i] = ExtConc[i];

	T_ex 	= ExtTemp;
	P_ex 	= ExtPressure;
// ****************Теплоемкости компонент**************
	Cp[H2]	=	16000.;	Cp[O2]	=	1670.;	Cp[N2]	=	1300;
	Cp[H2O]	=	0.; 	Cp[OH]	=	0.;		Cp[H]	=	0.;
	Cp[O]	=	0.; 	Cp[HO2]	=	0.; 	Cp[H2O2]=	0.;

}

Evaporation::~Evaporation(){
	delete[] Y_ex;
	delete[] Cp; 				// Теплоемкости компонент
	delete[] Y_w;
}


// ***Расчет теплоемкости смеси по теплоемкостям и концентрациям компонент на бесконечности***
double Evaporation::GetCp_mixt_ex(){
	double Cp_mixt_ex = 0.;
	for (int i = 0; i < CN; ++i)
		Cp_mixt_ex+=Y_ex[i]*Cp[i];
	return Cp_mixt_ex;
}

// ***Расчет теплоемкости смеси по теплоемкостям и концентрациям компонент на границе капли***
double Evaporation::GetCp_mixt_w(){
	double GetCp_mixt_w = 0.;
	for (int i = 0; i < CN; ++i)
		GetCp_mixt_w+=Y_w[i]*Cp[i];
	return GetCp_mixt_w;
}

// ***Расчет тепла парообразования для заданной температуры на поверхности капли***
double Evaporation::GetVapHeat(){
	double HL=0.;
	if(T_w<=T_c )
		HL=1000.*(1-0.006468*T_w)/(0.003909-0.000022*T_w);
	return HL;
}
// ******************Расчет кси************************
double Evaporation::GetXi(){
	double xi;
	double Cpe 	=	GetCp_mixt_ex();
	double Cpw 	=	GetCp_mixt_w();
	double HL	=	GetVapHeat();
	xi=log(1.+(Cpe*T_ex-Cpw*T_w)/(HL+Cpw*(T_w-T_av)));
	return xi;
}


int main(int argc, char  *argv[])
{	
	int CN = 9;
	double * Ye=new double[CN];
		for (int i = 0; i < CN; ++i)
		{
			if(i == 0) Ye[i] = 0.2;
			else if(i == 1) Ye[i] = 0.8;
			else Ye[i] = 0;
		}
	double Te = 150., Pe = 1e5;
	Evaporation dropplet(Ye, Te, Pe);
	return 0;
}