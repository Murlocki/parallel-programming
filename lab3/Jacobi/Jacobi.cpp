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
    }
    return res;
}

double Simpson_parallel(double a, double b, int N,int num_thr) {
    double step = (double)(b - a) / N;
    double* X_int = new double[N + 1];
    X_int[0] = a;
    for (int i = 1; i < N + 1; i++) {
        X_int[i] = X_int[i - 1] + step;
    }
    double res = 0;
    #pragma omp parallel num_threads(num_thr) reduction(+:res)
        {
            #pragma omp for collapse(1)
            for (int i = 0; i < N - 1; i = i + 2) 
            {
                res += step / 3 * (func_to_calc(X_int[i]) + 4 * func_to_calc((double)(X_int[i + 1])) + func_to_calc(X_int[i + 2]));
            }
        };
    return res;
}
double Simpson_parallel_2(double a, double b, int N,int num_thr) {
    double step = (double)(b - a) / N;
    double* X_int = new double[N + 1];
    X_int[0] = a;
    for (int i = 1; i < N + 1; i++) {
        X_int[i] = X_int[i - 1] + step;
    }
    double* res_f = new double[3];
    double res = 0;
    #pragma omp parallel num_threads(num_thr) reduction(+:res)
    {
        #pragma omp for collapse(2) 
        for (int i = 0; i < N - 1; i = i + 2)
        {
            for (int j = 0; j < 3; j++) {
                res_f[j] = func_to_calc(X_int[i + j]);
            }
            res += step / 3 * (res_f[0] + 4 * res_f[1] + res_f[2]);
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

    double start_time = clock();
    cout << Simpson_parallel_2(a, b, N, 2) << endl;
    double end_time = clock();
    cout <<"Второй вид 2 ядра:"<<end_time - start_time << endl;
    start_time = clock();
    cout << Simpson_parallel(a, b, N, 2) << endl;
    end_time = clock();
    cout << "Первый вид 2 ядра:" << end_time - start_time << endl;

    start_time = clock();
    cout << Simpson_parallel_2(a, b, N, 3) << endl;
    end_time = clock();
    cout << "Второй вид 3 ядра:" << end_time - start_time << endl;
    start_time = clock();
    cout << Simpson_parallel(a, b, N, 3) << endl;
    end_time = clock();
    cout << "Первый вид 3 ядра:" << end_time - start_time << endl;

    start_time = clock();
    cout << Simpson_parallel_2(a, b, N,4)<<endl;
    end_time = clock();
    cout << "Второй вид 4 ядра:" << end_time - start_time << endl;
    start_time = clock();
    cout << Simpson_parallel(a, b, N,4)<<endl;
    end_time = clock();
    cout << "Первый вид 4 ядра:" << end_time - start_time << endl;
    start_time = clock();
    cout << Simpson(a, b, N) << endl;
    end_time = clock();
    cout << "Последовательный вид:" << end_time - start_time << endl;
}


