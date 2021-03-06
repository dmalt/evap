#include "iostream"
#include "Evaporation.h"
#include "math.h"
#define DEBUG 1
using namespace std;

enum{H2,O2,N2,H2O,OH,H,O,HO2,H2O2,	EC};     //


double pow(double val, int pow){
	double result;
	if(pow){
		result=val;
		for (int i = 1; i < pow; ++i)
			result*=val;
	}
	else if(!pow) result = 1.;
	return result;
}


Evaporation::Evaporation(double * ExtConc, double ExtTemp, double ExtPressure, double T_dropplet, double In){
/*****************Число компонент*******************/
	CN=9;							
	IN = In;
/*********Критическа температура для кислорода******/
	T_c = 154.77; 		

	T_av = T_dropplet;				

/*********Выделение памяти под массивы**************/
	Y_ex	= new double[CN];			
	Cp 		= new double[CN]; 			 
	Y_w		= new double[CN];
	mu 		= new double[CN];			
/*******Инициализация внешних концентраций**********/
	for (int i = 0; i < CN; ++i){
		Y_ex[i] = ExtConc[i];
		Y_w[i] 	= 0.;	
	}
	Y_w[O2] = 1.;

	T_ex 	= ExtTemp;
	P_ex 	= ExtPressure;
/*******************Теплоемкости компонент******************/
	Cp[H2]	=	14230.;	Cp[O2]	=	1670.;	Cp[N2]	=	1300;
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
	// return 14230.;
}





/***Расчет тепла парообразования для заданной температуры на поверхности капли***/
double Evaporation::GetVapHeat(double T){
	double HL=0.1;
	if(T<=T_c )
		HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);
	return HL;
}
/******************Расчет кси************************/
double Evaporation::GetT(double xi, double T_prev){
	double exi = exp(-xi);
	double Cpe 	=	GetCp_mixt_ex();
	double HL	= 	GetVapHeat(T_prev);
	if (DEBUG) cout <<"HL = "<< HL <<endl;
	double T = (Cpe * T_ex - HL * (1./exi - 1)) / (Cpe * exi + Cp[O2] * (1 - exi));
	return T;
}

/*********Расчет дисбаланса уравнения для температуры*******/
double Evaporation::GetDelta(double xi, double T){
	// if(T<0.) T=0.;
	// if (T>1000.) T=1000.;
	double delta;
	double exi    = exp(-xi);
	for (int i = 0; i < CN; ++i){
		if(i!=O2) Y_w[i] = Y_ex[i] * exi;
		else Y_w[i] = 1 - (1 - Y_ex[i]) * exi;
		if(DEBUG){
			cout<<Y_w[i]<<endl;
			cout<<Y_ex[i]<<endl;
		}
	}
	double Y_O2_w = 1 - (1 - Y_ex[O2]) * exi;
	double mu_e   = GetMolarMassEx();
	// double IN     = 0.1;											// IN - внешний параметр и его надо бы посчитать.
	double Cpe    = GetCp_mixt_ex();
	double XN = GetO2PartPres(T) / P_ex;						// Давление насыщенных паров кислорода, обезразмеренное на внешнее давление


	double mu_w_divby_mu_O2 = 1/(1 - (1 - mu[O2]/mu_e) * exi);
	// double mu_w_divby_mu_O2 = 1/(1 + 15*(1 - Y_ex[O2]) * exi);

	delta = Y_O2_w * mu_w_divby_mu_O2;									// Добавили к функции дисбаланса левую часть уравнения (15)
	double add = xi * IN * sqrt(T * mu_e/ (T_ex * mu[O2]));
	// double add = xi * IN *sqrt((0.117358/T_ex)*T/(1.+15*(1-Y_ex[O2])*exi));
	delta+= add;
	delta-= XN;
	return delta;
}

double Evaporation::GetO2PartPres(double T){
	double A[4] = {-3630740., 150530., -2048.72, 9.31122};
	// double A[4] = {-1649810., 79239., -1252., 6.5};
	double P_O2 = 0;
	for (int i = 0; i < 4; ++i)
		P_O2+=A[i] * pow(T,i);
	return P_O2;
}

int Evaporation::SolveNewton(){
	double xi_next, xi_prev, f_prime, f_next, f_prev;
	double dxi = 0.01, df;
	double T_prev =130., T_next = 130., T_const = 70.;
	xi_prev = 1.;				// Начальное приближение.
	xi_next = xi_prev - dxi;     	// 
	int iter = 0;
	const double a_tol = 1.e-5;
	do{
		/******Считаем производную********/

		
		T_prev     = GetT(xi_prev, T_const);
		T_next     = GetT(xi_next, T_const);
		if (DEBUG) cout<< "T_prev="<< T_prev<<" T_next = "<<T_next<< endl;
		f_next = GetDelta(xi_next, T_next);
		f_prev = GetDelta(xi_prev, T_prev);
		df = f_next - f_prev;
		dxi = xi_next - xi_prev;
		f_prime = df/dxi;
		/********************************/
		xi_prev = xi_next;
		f_prev = f_next;
		if(DEBUG) cout<<"f_prev = "<<f_prev<<" f_prime = "<<f_prime<<endl<<endl;
		xi_next = xi_prev - f_prev/f_prime;
		// if(T_next<70.) T_next=70.;
		iter++;
		cout<<".";
	}while(fabs(xi_next-xi_prev) > a_tol );
	cout<<"xi_next = "<<xi_next<<endl;
	T_w = GetT(xi_next, T_const);
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
			if(i == O2) Ye[i] = 0.2;
			else if(i == H2) Ye[i] = 0.8;
			else Ye[i] = 0;
		}

	double Te = 100., Pe = 5066250, Tav = 80.;
	// double In = 0.1;
	for (double In = 0.1; In < 100.; In += 0.1){
		Evaporation dropplet(Ye, Te, Pe, Tav, In);
		int iter = dropplet.SolveNewton();
 		cout<<"In "<<iter <<" iterations we`ve got T equal to "<<dropplet.T_w<<endl;
	}
	return 0;
}