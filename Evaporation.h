class Evaporation
{
private:
	/*************** VARIABLES **************/
	short CN; 			// Число компонент
	double * Cp; 		// Теплоемкости компонент

	double T_c;			// Критическая температура для кислорода
	double T_ex;		// Температура на бесконечности
	double T_av;		// Средняя температура капли

	double * Y_ex;		// Массив внешних концентраций
	double * Y_w;		// Массив концентраций на границе

	double P_ex;		// Внешнее давление

	
	double Peclet;		// Число Пекле

	/************** FUNCTIONS **************/
	double GetCp_mixt_ex(); 		// Расчет теплоемкости смеси на бесконечности
	double GetCp_mixt_w();			// Расчет теплоемкости смеси на границе
	double GetVapHeat(double T);	// Расчет теплоемкости парообразования в зависимости от температуры
	double GetXi(double T);			// Расчет кси
	double GetDelta(double T);		// Расчет дисбаланса уравнения, которое решаем методом Ньютона. В функцию передается температура на шаге метода Ньютона
	double GetO2PartPres(double T);	// Расчет парциального давления паров кислорода на границе капли в зависимости от температуры
	

public:
	double T_w;			// Температура на границе
	Evaporation(double * ExtConc, double ExtTemp, double ExtPres);
	~Evaporation();
	int	SolveNewton();		// Расчет температуры методом Ньютона; возвращает число итераций
};