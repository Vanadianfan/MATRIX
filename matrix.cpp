#include <string>
#include <cstring>
#include <unistd.h>
#include <libgen.h>
#include <sstream>
// #include <matrix.h>  // optioanl, if already put matrix.h in your include path.
#include "matrix.h"
using namespace std;
int main(void) {
    #ifdef __APPLE__
    {
        char pathbuf[1024];
        strcpy(pathbuf, __FILE__);
        chdir(dirname(pathbuf));
    }
    #endif
    freopen("test_data.txt", "r", stdin);

    frac a1(12,25), a2(6,5), b1(-3,5), b2(1), c1(-16,25), c2(-8,5);
    cout << "12/25 * 6/5 + -3/5 * 1 + -16/25 * -8/5 = " << a1*a2+b1*b2+c1*c2 << endl;
    cout << "12/25 + -16/25 = " << a1 + c1 << endl << endl;

    Square_Matrix<frac> frac_identity = Square_Matrix<frac>::identity(3);
    istringstream frac_input("1 2 3 4 5 6 7 8 9");
    Square_Matrix<frac> frac_matrix(3);
    frac_input >> frac_matrix;

    cout << "Square_Matrix<frac> identity(2):" << endl << frac_identity;
    cout << "Square_Matrix<frac> sample matrix:" << endl << frac_matrix;
    cout << "det(sample) = " << frac_matrix.determinant() << endl;
    cout << "identity * sample:" << endl << frac_identity * frac_matrix;
    cout << "sample * identity:" << endl << frac_matrix * frac_identity;

    Square_Matrix<frac> frac_cofactor = frac_matrix.cofator_matrix();
    cout << "cofactor(sample):" << endl << frac_cofactor;
}
