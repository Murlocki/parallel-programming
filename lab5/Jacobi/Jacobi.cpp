#include<iostream>
#include<string>
#include<sstream>
#include<omp.h>
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

double** rand_matr(int N) {
    double** res = new double* [N];
    for (int i = 0; i < N; i++) {
        res[i] = new double[N];
        for (int j = 0; j < N; j++) {
            res[i][j] = rand() % 30 + 1;
        }
    }
    return res;
}

double** seq_fload_yor(double** matr, int N) {
    double** res = new double* [N];

    for (int i = 0; i < N; i++) {
        res[i] = new double[N];
        for (int j = 0; j < N; j++) {
            if (matr[i][j] == 0 && i!=j) {
                res[i][j]= numeric_limits<double>::infinity();
            }
            else { 
                res[i][j] = matr[i][j]; 
            }
        }
    }
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (res[i][j] > res[i][k] + res[k][j]) {
                    res[i][j] = res[i][k] + res[k][j];
                }
            }
        }
    }
    return res;
}

double** parallel_fload_yor(double** matr, int N,int numberthreads) {
    double** res = new double* [N];

    for (int i = 0; i < N; i++) {
        res[i] = new double[N];
        for (int j = 0; j < N; j++) {
            if (matr[i][j] == 0 && i != j) {
                res[i][j] = numeric_limits<double>::infinity();
            }
            else {
                res[i][j] = matr[i][j];
            }
        }
    }
    for (int k = 0; k < N; k++) {
#pragma omp parallel for collapse(2) num_threads(numberthreads)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (res[i][j] > res[i][k] + res[k][j]) {
                    res[i][j] = res[i][k] + res[k][j];
                }
            }
        }
    }
    return res;
}


int main()
{
    setlocale(LC_ALL, "russian");
    //double a, b;
    //int N;
    //cout << "Количество вершин графа массива" << endl;
    //cin >> N;
    //cout << "Введите матрицу смежности" << endl;
    //double** matr_sm = input_matr(N, N);
    for (int N = 100; N <= 1000; N = N + 50) {
        double** matr_sm = rand_matr(N);
        //print_matr(matr_sm, N, N);

        cout << N << endl;

        double start = clock();
        double end = clock();

        start = clock();
        double** res = seq_fload_yor(matr_sm, N);
        end = clock();
        cout << "Матрица кратчайших путей" << endl;
        //print_matr(res, N, N);
        cout << end - start << endl;

        cout << "Кратчашиие для пар 2" << endl;
        start = clock();
        res = parallel_fload_yor(matr_sm, N, 2);
        end = clock();
        //print_matr(res, N, N);
        cout << end - start << endl;

        cout << "Кратчашиие для пар 3" << endl;
        start = clock();
        res = parallel_fload_yor(matr_sm, N, 3);
        end = clock();
        //print_matr(res, N, N);
        cout << end - start << endl;

        cout << "Кратчашиие для пар 4" << endl;
        start = clock();
        res = parallel_fload_yor(matr_sm, N, 4);
        end = clock();
        //print_matr(res, N, N);
        cout << end - start << endl;
        cout << endl;
    }

}


