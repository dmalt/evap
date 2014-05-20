#include "iostream"
#include "math.h"
using namespace std;

double T_ex=120,T_av=70;

double GetVapHeat(double T){
	double HL=0.1;
	if(T<=154.6 )
		HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);
	return HL;
}
double GetCp_mixt_ex(){
	// double Cp_mixt_ex = 0.;
	// for (int i = 0; i < CN; ++i)
	// 	Cp_mixt_ex+=Y_ex[i]*Cp[i];
	// return Cp_mixt_ex;
	return 14000.;
}
double GetCp_mixt_w(){
	// double GetCp_mixt_w = 0.;
	// for (int i = 0; i < CN; ++i)
	// 	GetCp_mixt_w+=Y_w[i]*Cp[i];
	// return GetCp_mixt_w;
	return 1670.;
}
double GetXi(double T){
	double xi;
	double Cpe 	=	GetCp_mixt_ex();
	double Cpw 	=	GetCp_mixt_w();
	double HL	=	GetVapHeat(T);
	cout<<"GetXi:"<<(Cpe*T_ex-Cpw*T)<<" "<<(HL+Cpw*(T-T_av))<<endl;
	xi=log(1.+(Cpe*T_ex-Cpw*T)/(HL+Cpw*(T-T_av)));
	return xi;
}

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


double GetO2PartPres(double T){
	double A[4] = {-3630740., 150530., -2048.72, 9.31122};
	double P_O2 = 0;
	for (int i = 0; i < 4; ++i)
		P_O2+=A[i] * pow(T,i);
	return P_O2;
}
int main(int argc, char const *argv[])
{
	for (double T = 0.; T < 200; ++T)
	{
		cout<<"T = "<< T<< " P = "<<GetO2PartPres(T)<<" HL = "<<GetVapHeat(T)<<"  xi = "<<GetXi(T)<<endl;
	}
	return 0;
}