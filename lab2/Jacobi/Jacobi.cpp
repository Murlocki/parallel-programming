#include<iostream>
#include<string>
#include<sstream>
using namespace std;
double** input_matr(int N, int M) {
    double** matr = new double* [N];
    for (int i = 0; i < N; i++) {
        matr[i] = new double[M];
    }
    for (int i = 0; i < N + 1; i++) {
        string inp = "";
        getline(std::cin, inp);
        stringstream ss(inp);
        int j = 0;
        string word;
        while (ss >> word) {
            matr[i - 1][j] = stod(word);
            j = j + 1;
        }
    }
    return matr;
}
void print_matr(double** matr, int N, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << matr[i][j] << ' ';
        }
        cout << endl;
    }
}

void Jacobi(int N, double** A, double* F, double* X)
{
    double* TempX = new double[N];
    double norm;
    double eps = 0.0001;
    do {
        for (int i = 0; i < N; i++) {
            TempX[i] = F[i];
            for (int g = 0; g < N; g++) {
                if (i != g)
                    TempX[i] -= A[i][g] * X[g];
            }
            TempX[i] /= A[i][i];
        }
        norm = fabs(X[0] - TempX[0]);
        for (int h = 0; h < N; h++) {
            if (fabs(X[h] - TempX[h]) > norm)
                norm = fabs(X[h] - TempX[h]);
            X[h] = TempX[h];
        }
    } while (norm > eps);
    for (int i = 0; i < N; i++) {
        cout << X[i] << " ";
    }
    delete[] TempX;
}

int main()
{
    setlocale(LC_ALL, "russian");
    int N1;
    cin >> N1;
    double** matr1 = input_matr(N1, N1 + 1);
    print_matr(matr1, N1, N1 + 1);
    double* start = new double[N1];
    for (int i = 0; i < N1; i++) {
        start[i] = rand() % 100;
    }
    double* F = new double[N1];
    for (int i = 0; i < N1; i++) {
        F[i] = matr1[i][N1];
    }
    Jacobi(N1,matr1, F,start);
    


}