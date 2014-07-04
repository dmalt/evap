



#include  <iostream>
#include  <string> 
#include "iomanip"

using namespace std; 

int main(void) { 
        string TheGood = "Jekyll", TheBad = "Hyde"; 
        cout << TheGood + " & " + TheBad << endl; 
        cout << TheBad + " & " + TheGood << endl; 
        cout<<scientific<<5.5<<endl;


        enum Conc {H2,O2,N2,H2O,OH,H,O,HO2,H2O2,	EC}; 
        Conc C;	
        cout << Conc(O2);
        return 0; 
} 