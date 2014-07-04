#include "fstream"
#include "iostream"
#include "Evaporation.h"
#include "math.h"
#include <iomanip>
#define DEBUG 0
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


Evaporation::Evaporation(double * ExtConc, double ExtTemp, double ExtPressure, double T_droplet, double In){
/*****************Число компонент*******************/
	CN=9;		
	IN = In;					
/*********Критическа температура для кислорода******/
	T_c = 154.77; 		

	T_av = T_droplet;				

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


/***Расчет теплоемкости смеси по теплоемкостям и концентрациям компонент на границе капли***/
double Evaporation::GetCp_mixt_w(){

	double Cp_mixt_w = 0.;
	for (int i = 0; i < CN; ++i)
		Cp_mixt_w+=Y_w[i]*Cp[i];
	return Cp_mixt_w;
	 // return 1670.;

}

/***Расчет тепла парообразования для заданной температуры на поверхности капли***/
double Evaporation::GetVapHeat(double T){
	double HL=0.1;
	if(T<=T_c )
		HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);
	// cout<<" HL = "<<HL<<endl;
	return HL;
}

double GetVapHeat(double T){
	double HL=0.1;
	double T_c = 154.77;
	if(T<=T_c )
		HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);
	// cout<<" HL = "<<HL<<endl;
	return HL;
}
/******************Расчет кси************************/
double Evaporation::GetXi(double T){
	double xi;
	double alpha = 0.6;   // Параметр альфа нужен для расчета теплоемкости 
						  // на границе как взвешенной суммы теплоекости кислорода и водорода 
	// double Cpe 	=	GetCp_mixt_ex();
	// double Cpw 	=	GetCp_mixt_w();
	double Cpe;
	double Cpw;
	Cpe = Cpw = Cp[H2] * alpha + Cp[O2] * (1 - alpha);
	double HL	=	GetVapHeat(T);
	if(DEBUG) cout<<"GetXi:"<<(Cpe*T_ex-Cpw*T)<<" "<<(HL+Cpw*(T-T_av))<<" HL="<<HL<<endl;
	xi=log(1.+(Cpe*T_ex-Cpw*T)/(HL+Cpw*(T-T_av)));
	return xi;
}


void Evaporation::GetConc(double xi){
	double exi = exp(-xi);
	for (int i = 0; i < CN; ++i){
		if(i!=O2) Y_w[i] = Y_ex[i] * exi;
		else Y_w[i] = 1 - (1 - Y_ex[i]) * exi;
		if(DEBUG){
			cout<<Y_w[i]<<endl;
			cout<<Y_ex[i]<<endl;
		}
	}
}

/*********Расчет дисбаланса уравнения для температуры*******/
double Evaporation::GetDelta(double T){
	// if(T<0.) T=0.;
	// if (T>T_c) T=T_c;
	double delta;
	double xi     = GetXi(T);
	double exi    = exp(-xi);
	GetConc(xi);
	double Y_O2_w = 1 - (1 - Y_ex[O2])*exi;
	double mu_e   = GetMolarMassEx();
	// double IN     = 5.1;										// IN - внешний параметр и его надо бы посчитать.
	double Cpe    = GetCp_mixt_ex();
	double Cpw    = GetCp_mixt_w();
	double Gamma  = 1.;											// Эту гамму нужно посчитать, хотя вроде бы она примерно равна 1
	double XN = GetO2PartPres(T) / P_ex;						// Давление насыщенных паров кислорода, обезразмеренное на внешнее давление


	double mu_w_divby_mu_O2 = 1/(1 - (1 - mu[O2]/mu_e) * exi);
	// double mu_w_divby_mu_O2 = 1/(1 + 15*(1 - Y_ex[O2]) * exi);

	delta = Y_O2_w * mu_w_divby_mu_O2;									// Добавили к функции дисбаланса левую часть уравнения (15)
	double add = xi * IN * sqrt(T * Cpw * mu_w_divby_mu_O2 * Gamma / (Cpe * T_ex));
	// double add = xi * IN *sqrt((0.117358/T_ex)*T/(1.+15*(1-Y_ex[O2])*exi));

	if(DEBUG == 1){
		cout<<"add = "<<add<<endl<<"T= "<<T<<endl<<"Cpw = "<<Cpw<<endl<<"Cpe="<<Cpe<<endl<<"mu_div = "<<mu_w_divby_mu_O2<<endl; 
		cout<<"exi = "<<exi<<endl<<"xi = "<<xi<<endl<<"delta="<<delta<<endl<<"XN="<<XN<<endl<<endl;
	}

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
	double T_next, T_prev, f_prime, f_next, f_prev;
	double dT = 0.1, df;
	T_prev = 130.;				// Начальное приближение.
	T_next = T_prev - dT;     	// 
	int iter = 0;
	const double a_tol = 1.e-5;
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
		if(DEBUG) cout<<"f_prev = "<<f_prev<<" f_prime = "<<f_prime<<endl<<endl;
		T_next = T_prev - f_prev/f_prime;
		// if(T_next<70.) T_next=70.;
		iter++;
		cout<<".";
	}while(fabs(T_next-T_prev) > a_tol );
	cout<<endl;

	if(T_next<T_c)
		T_w = T_next;
	else if(T_next>=T_c){
		T_w = T_c;
		GetConc(GetXi(T_w));
	}

	Peclet = GetXi(T_w);
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
	ofstream out;
	out.open("out.txt");
	// out.precision(5);
	 	
	out<<endl;
	int CN = 9;
	double * Ye=new double[CN];
		for (int i = 0; i < CN; ++i)
		{
			if(i == O2) Ye[i] = 0.3;
			else if(i == H2) Ye[i] = 0.7;
			else Ye[i] = 0;
		}

	double Te = 170., Pe = 15.e5, Tav = 60.;
	double T_boil = 91.;
	double Cp = 4046.;
	double HL = GetVapHeat(T_boil);
	double exi_boil = 1./(1.+(Cp*Te-Cp*T_boil)/(HL+Cp*(T_boil-Tav)));
	//**********************************//
	double Y_H2_boil = Ye[H2] * exi_boil;
	double Y_O2_boil = 1 - (1 - Ye[O2]) * exi_boil; 
	out<<"Y_H2_boil = " << Y_H2_boil<<"\t Y_O2_boil = "<< Y_O2_boil << endl;
	// ******************************** //
	// ********* Шапка вывода *******  //
	out<<""<<setw(3)<<"In";
		// out.width(9);
		out<<"\t"<<setw(5)<<"T_w";
		out<<"\t"<<setw(5)<<"Te";
		out<<"\t"<<setw(6)<<"Peclet";
		 // out.width(7);
		string ConcName[EC] = {"H2","O2","N2","H2O","OH","H","O","HO2","H2O2"};

		for (int i = 0; i < CN; ++i){
			 // out.precision(3);
			 // out.width(4);	
			// enum{H2,O2,N2,H2O,OH,H,O,HO2,H2O2,	EC}; 

			out<<"\t"<<setw(3)<< "["<<ConcName[i]<<"]";
		}
		out<<endl<<endl;
	// ***************************** //

	for (double In = 0.; In < 100.; In += 1.){
		Evaporation droplet(Ye, Te, Pe, Tav,In);
		int iter = droplet.SolveNewton();
		cout<<"In "<<iter <<" iterations we`ve got T equal to "<<droplet.T_w<<endl;
		// out.width(4);
		 out.precision(4);
		out<<""<<setw(3)<<In;
		// out.width(9);
		out<<"\t"<<setw(5)<<droplet.T_w;
		out<<"\t"<<setw(5)<<Te;
		out<<"\t"<<setw(6)<<droplet.Peclet;
		 // out.width(7);
		for (int i = 0; i < CN; ++i){
			 // out.precision(3);
			  // out.width(6);	

			out<<"\t"<<setw(6)<< droplet.Y_w[i];
		}
		out<<endl;
	}
	out.close();
	return 0;
}