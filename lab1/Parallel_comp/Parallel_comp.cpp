
#include <iostream>
using namespace std;
#include<string>
#include<sstream>
#include<vector>
#include<omp.h>
#include<cmath>
#include<ctime>

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

double** seq_matr_multip(double** matr1, int N1, int M1, double** matr2, int N2, int M2) {
    double** result = new double* [N1];
    for (int i = 0; i < N1; i++) {
        result[i] = new double[M2];
    }

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < M2; j++) {
            double s = 0;
            for (int m = 0; m < M1; m++) {
                s = s + matr1[i][m] * matr2[m][j];
            }
            result[i][j] = s;
        }
    }
    return result;
}

double** lent_mult(double** matr1, double** matr2, int N1, int M1, int N2, int M2, int threads) {
    double** result = new double* [N1];
    for (int i = 0; i < N1; i++) {
        result[i] = new double[M2];
        for (int j = 0; j < M2; j++) {
            result[i][j] = 0;
        }
    }
    omp_set_num_threads(threads);
#pragma omp parallel for collapse(3) num_threads(threads)
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < M2; j++) {
            for (int k = 0; k < M1; k++) {
                result[i][j] = result[i][j] + matr1[i][k] * matr2[k][j];
            }
        }
    }

    return result;
}


double** block_mult(double** matr1, double** matr2, int N1, int M1, int N2, int M2,int threads) {
    double** result = new double* [N1];
    for (int i = 0; i < N1; i++) {
        result[i] = new double[M2];
        for (int j = 0; j < M2; j++) {
            result[i][j] = 0;
        }
    }
    omp_set_num_threads(threads);
    int blockSize = max(M1 / threads, 1);
    cout << "Размер блока: " << blockSize << endl;
#pragma omp parallel for collapse(3) num_threads(threads)
    for (int i = 0; i < N1; i += blockSize) {
        for (int j = 0; j < M2; j += blockSize) {
            for (int k = 0; k < M1; k += blockSize) {
                // Умножение блоков матриц
                for (int ii = i; ii < min(i + blockSize, N1); ii++) {
                    for (int jj = j; jj < min(j + blockSize, M2); jj++) {
                        for (int kk = k; kk < min(k + blockSize, M1); kk++) {
                            result[ii][jj] += matr1[ii][kk] * matr2[kk][jj];
                        }
                    }
                }
                
            }
        }
    }
    return result;
}

int main()
{
    setlocale(LC_ALL, "russian");
    int N1, M1;
    cin >> N1 >> M1;
    double** matr1 = input_matr(N1, M1);

    int N2, M2;
    cin >> N2 >> M2;
    double** matr2 = input_matr(N2, M2);

    double time = clock();
    double**result=seq_matr_multip(matr1, N1, M1, matr2, N2, M2);
    double end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);

    cout << "Блок 2 нить" << endl;
    time = clock();
    result = block_mult(matr1, matr2, N1, M1, N2, M2,2);
    end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);
    cout << "Блок 3 нить" << endl;
    time = clock();
    result = block_mult(matr1, matr2, N1, M1, N2, M2,3);
    end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);
    cout << "Блок 4 нить" << endl;
    time = clock();
    result = block_mult(matr1, matr2, N1, M1, N2, M2,4);
    end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);

    cout << "лента 2 нить" << endl;
    time = clock();
    result = lent_mult(matr1, matr2, N1, M1, N2, M2,2);
    end_time = clock();
    cout << "Result time: " << double((end_time-time))<< endl;
    print_matr(result, N1, M2);
    cout << "лента 3 нить" << endl;
    time = clock();
    result = lent_mult(matr1, matr2, N1, M1, N2, M2,3);
    end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);
    cout << "лента 4 нить" << endl;
    time = clock();
    result = lent_mult(matr1, matr2, N1, M1, N2, M2,4);
    end_time = clock();
    cout << "Result time: " << double((end_time - time)) << endl;
    print_matr(result, N1, M2);


}

