#include "iostream"
#include "Evaporation.h"
#include "math.h"

using namespace std;

enum{H2,O2,N2,H2O,OH,H,O,HO2,H2O2,	EC};     //


double pow(double val, int pow){
	double power = val;
	for (int i = 1; i < pow; ++i)
		power*=val;
	return power;
}


Evaporation::Evaporation(double * ExtConc, double ExtTemp, double ExtPressure){
/*****************Число компонент*******************/
	CN=9;							
/*********Критическа температура для кислорода******/
	T_c = 154.77; 						
/*********Выделение памяти под массивы**************/
	Y_ex	= new double[CN];			
	Cp 		= new double[CN]; 			 
	Y_w		= new double[CN];			
/*******Инициализация внешних концентраций**********/
	for (int i = 0; i < CN; ++i)
		Y_ex[i] = ExtConc[i];

	T_ex 	= ExtTemp;
	P_ex 	= ExtPressure;
/*******************Теплоемкости компонент******************/
	Cp[H2]	=	16000.;	Cp[O2]	=	1670.;	Cp[N2]	=	1300;
	Cp[H2O]	=	0.; 	Cp[OH]	=	0.;		Cp[H]	=	0.;
	Cp[O]	=	0.; 	Cp[HO2]	=	0.; 	Cp[H2O2]=	0.;
/******************Молярные массы компонент*****************/
	mu[H2]	=	2.e-3;		mu[O2]	=	32.e-3;	mu[N2]	=	28.e-3;
	mu[H2O]	=	18.e-3; 	mu[OH]	=	17.e-3;	mu[H]	=	1.e-3;
	mu[O]	=	16.e-3; 	mu[HO2]	=	33.e-3; mu[H2O2]=	34.e-3;

}

Evaporation::~Evaporation(){
	delete[] Y_ex;
	delete[] Cp; 				// Теплоемкости компонент
	delete[] Y_w;
}


/***Расчет теплоемкости смеси по теплоемкостям и концентрациям компонент на бесконечности***/
double Evaporation::GetCp_mixt_ex(){
	double Cp_mixt_ex = 0.;
	for (int i = 0; i < CN; ++i)
		Cp_mixt_ex+=Y_ex[i]*Cp[i];
	return Cp_mixt_ex;
}

/***Расчет теплоемкости смеси по теплоемкостям и концентрациям компонент на границе капли***/
double Evaporation::GetCp_mixt_w(){
	double GetCp_mixt_w = 0.;
	for (int i = 0; i < CN; ++i)
		GetCp_mixt_w+=Y_w[i]*Cp[i];
	return GetCp_mixt_w;
}

/***Расчет тепла парообразования для заданной температуры на поверхности капли***/
double Evaporation::GetVapHeat(double T){
	double HL=0.;
	if(T_w<=T_c )
		HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);
	return HL;
}
/******************Расчет кси************************/
double Evaporation::GetXi(double T){
	double xi;
	double Cpe 	=	GetCp_mixt_ex();
	double Cpw 	=	GetCp_mixt_w();
	double HL	=	GetVapHeat(T);
	xi=log(1.+(Cpe*T_ex-Cpw*T)/(HL+Cpw*(T-T_av)));
	return xi;
}

/*********Расчет дисбаланса уравнения для температуры*******/
double Evaporation::GetDelta(double T){
	double delta;
	double xi     = GetXi(T);
	double exi    = exp(-xi);
	double Y_O2_w = 1 - (1 - Y_ex[O2])*exp(-xi);
	double mu_e   = GetMolarMassEx();
	double IN     = 0.1;											// IN - внешний параметр и его надо бы посчитать.
	double Cpe    = GetCp_mixt_ex();
	double Cpw    = GetCp_mixt_w();
	double Gamma  = 1.;												// Эту гамму нужно посчитать, хотя вроде бы она примерно равна 1
	double XN = GetO2PartPres(T) / P_ex;							// Давление насыщенных паров кислорода, обезразмеренное на внешнее давление

	double mu_w_divby_mu_O2 = 1/(1 - (1 - mu_e/mu[O2]) * exi);

	delta = Y_O2_w * mu_w_divby_mu_O2;									// Добавили к функции дисбаланса левую часть уравнения (15)
	delta+= xi * IN * sqrt(T * Cpw * mu_w_divby_mu_O2 * Gamma / (Cpe * T_ex));
	delta-= XN;
	return delta;
}

double Evaporation::GetO2PartPres(double T){
	double A[4] = {-1649810., 79239., -1252., 6.5};
	double P_O2 = 0;
	for (int i = 0; i < 4; ++i)
		P_O2+=A[i] * pow(T,i);
	return P_O2;
}

int Evaporation::SolveNewton(){
	double T_next, T_prev, f_prime, f_next, f_prev;
	double dT = 0.1, df;
	T_prev = 150.;
	T_next = T_prev - dT;     // Начальное приближение
	int iter = 0;
	const double a_tol = 1.e-7;
	do{
		/******Считаем производную********/
		f_next = GetDelta(T_next);
		f_prev = GetDelta(T_prev);
		df = f_next - f_prev;
		dT = T_next - T_prev;
		f_prime = df/dT;
		/********************************/
		T_prev = T_next;
		f_prev = f_next;
		T_next = T_prev - f_prev/f_prime;
		iter++;
		cout<<".";
	}while(fabs(T_next-T_prev) > a_tol && iter < 500);
	cout<<endl;
	T_w = T_next;
	return iter;
}

double Evaporation::GetMolarMassEx(){
	double mu_e=0;
	for (int i = 0; i < CN; ++i)
		mu_e+=Y_ex[i]/mu[i];
	mu_e = 1./mu_e;
	return mu_e;
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
	int iter = dropplet.SolveNewton();
	cout<<"In "<<iter <<" iterations we`ve got T equal to "<<dropplet.T_w<<endl;
	return 0;
}