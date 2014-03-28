/* get the thermal data for LOGOS */


#include <stdio.h>
#include <math.h>
#include <string.h>

#define A0	(-1649810.)
#define A1	(79239.)
#define A2	(-1252.)
#define A3	6.5

#define Ye 0.2
#define Te 160
#define pe 5066250
#define Tcp 130




#define Tc	154.6

double Fi1(double xi, double T, double In) {
 double p,ex=exp(-xi);
 p=(1.-(1-Ye)*ex)/(1.+15*(1-Ye)*ex);
 p+=xi*In*sqrt((0.117358/Te)*T/(1.+15*(1-Ye)*ex));
 p-=(A0+T*(A1+T*(A2+T*A3)))/pe;
 return p;
}

double Fi2(double T) {   //Функция считает кси
  //Рассчитываем теплоту парообразования для заданной температуры
 double HL=1000.*(1-0.006468*T)/(0.003909-0.000022*T);  //Тут наверное нужно поправить, чтобы лишний раз не считалось зря
 if(T>Tc || HL<=0.1) HL=0.1;   
 double cpe=14230., cpw=1670.; //Здесь cpe вообще-то нужно задавать извне;cpw нужно рассчитывать
 double xi=1.+(cpe*Te-cpw*T)/(HL+cpw*(T-Tcp));                       
 return xi;
}

double CritRes(double xi, double In) {
 double p,ex=exp(-xi);
 p=(1.-(1-Ye)*ex)/(1.+15*(1-Ye)*ex);
 p+=xi*In*sqrt((0.117358*(Tcp/Te)+(1-0.117358*(Tcp/Te))*ex)/(1.+15*(1-Ye)*ex))-4695670./pe;
 return p;
}

double GetRes(double T, double *ex, double In) {
 *ex=log(Fi2(T));  // Посчитали кси
 return Fi1(*ex,T,In);
}

int GetResultSub(double *xi, double *T, double In) {
 int it=0;
 double r0,r1;
 double x0,x1,dx0,dx1;
 double dT=0.1;
 r0=GetRes(x0=*T,xi,In);
 x1=x0-(dx0=dT);                    //Зачем?
 do {
  r1=GetRes(x1,xi,In);  //Здесь x1 выступает в роли температуры
  dx1=r1*(x1-x0)/(r1-r0);
  if(x1-dx1<0.5*x1) dx1=0.5*x1;     //И это что за костыль?
  ++it;
//  printf("it=%d x0=%e x1=%e r0=%e r1=%e dx=%e\n",it,x0,x1,r0,r1,dx1);
  x0=x1;r0=r1;x1-=dx1;dx0=dx1;
 } while(fabs(dx1)>1.e-5);
 *T=x1;
 return it;
}

int GetResultCrit(double *xi, double *T, double In) {
 int it=0;
 double r0,r1;
 double x0,x1,dx0=0.01,dx1;
 *T=Tc;x0=*xi;
 r0=CritRes(x0,In);
 x1=x0-dx0;
 do {
  r1=CritRes(x1,In);
  dx1=r1*(x1-x0)/(r1-r0);
//  if(x1-dx1<0.5*x1) dx1=0.5*x1;
  ++it;
//  printf("it=%d x0=%e x1=%e r0=%e r1=%e dx=%e\n",it,x0,x1,r0,r1,dx1);
  x0=x1;r0=r1;x1-=dx1;dx0=dx1;
 } while(fabs(dx1)>1.e-5);
 *xi=x1;
 return it;
}

#define IN	101

#define Inmax	1.
double Ain[IN],Xis[IN],Tes[IN];
int IT[IN];

int main(void) {
 int k,it,ince=0;
 double xi=2.,T=150.,In=10. ,di=Inmax/(IN-1.);  //xi - это переменная кси, равная интегралу от Пекле. В нашей системе уравнений выступает в роли неизвестной
 for(k=0;k<IN;k++) {                   //В цикле меняем параметр In
  In=k*di;Ain[k]=In;                                                        //Какая-то ерунда
  if(!ince) {                                                               //Реализуется две ветки расчета - до- и сверхкритическая 
    //Функция GetResultSub меняет значения 
   IT[k]=GetResultSub(&xi,&T,In);Xis[k]=xi;Tes[k]=T;
   if(T>=Tc) ince=1;                                                        // Переход в сверхкритическую ветку расчета 
  }
  if(ince) {                                                                // Сверхкритическая ветка
   IT[k]=GetResultCrit(&xi,&T,In);Xis[k]=xi;Tes[k]=T;
  }
 }
 

 printf("Result= [\n");
 for(k=0;k<IN;k++) {
  printf(" %15.10f %15.10f %15.10f;\n",Ain[k],Xis[k],Tes[k]);
 }
 printf("];\n\n");

 return 0;
}



