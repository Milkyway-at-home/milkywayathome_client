#include <iomanip>
#include <iostream>

using namespace std;

int main()
{

    double a = 3000000000000000.;

    for( int i = 0; i < 32; i++, a /= 10. )
        cout << setprecision(10) << a << endl;

    return 0;
    
}
