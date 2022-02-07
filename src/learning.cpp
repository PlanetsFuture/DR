#include "dr.h"
#include "lu_decomposition.h"

int main() {
    int n = 3;
    dr_matrix_double a(n, n);
    dr_vector_double b(n), x(n);

    ifstream fin;
    ofstream fout;
    fin.open("data.txt");
    // read data and save it in the matrix a
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fin >> a[i][j];

    cout << "matrix = " << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j] << '\t';
        }
        cout << endl;
    }
    return 0;
}
