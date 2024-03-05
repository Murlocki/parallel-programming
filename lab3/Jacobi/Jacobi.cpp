#include<iostream>
#include<string>
#include<sstream>
#include<omp.h>
using namespace std;

double func_to_calc(double x) {
    return x*x;
}

double Simpson(double a, double b, double N) {
    double step = (double)(b - a) / N;
    double* X_int = new double[N + 1];
    X_int[0] = a;
    for (int i = 1; i < N + 1; i++) {
        X_int[i] = X_int[i - 1] + step;
    }
    double res = 0;
    for (int i = 0; i < N-1; i=i+2) {
        res += step / 3 * (func_to_calc(X_int[i]) + 4 * func_to_calc((double)(X_int[i + 1])) + func_to_calc(X_int[i + 2]));
        cout << func_to_calc(X_int[i]) << " " << 4 * func_to_calc((double)(X_int[i + 1])) << " " << func_to_calc(X_int[i + 2]) << endl;
    }
    return res;
}

double Simpson_parallel(double a, double b, int N) {
    double step = (double)(b - a) / N;
    double* X_int = new double[N + 1];
    X_int[0] = a;
    for (int i = 1; i < N + 1; i++) {
        X_int[i] = X_int[i - 1] + step;
    }
    double res = 0;
    #pragma omp parallel num_threads(3)
        {
            #pragma omp for 
            for (int i = 0; i < N - 1; i = i + 2) 
            {
                res += step / 3 * (func_to_calc(X_int[i]) + 4 * func_to_calc((double)(X_int[i + 1])) + func_to_calc(X_int[i + 2]));
                cout << func_to_calc(X_int[i]) << " " << 4 * func_to_calc((double)(X_int[i + 1])) << " " << func_to_calc(X_int[i + 2]) << endl;
            }
        };
    return res;
}

int main()
{
    setlocale(LC_ALL, "russian");
    double a,b;
    int N;
    cout << "Начало отрезка" << endl;
    cin >>a;
    cout << "Конец отрезка" << endl;
    cin >> b;
    cout << "Количество отрезков для разбиения" << endl;
    cin >> N;
    cout << Simpson_parallel(a, b, N);
    cout << Simpson(a, b, N);
}


