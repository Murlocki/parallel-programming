#include<iostream>
#include<string>
#include<sstream>
#include<omp.h>
using namespace std;

double* createArray(int N) {
    double* res = new double[N];
    for (int i = 0; i < N; i++) {
        res[i] = rand() % 10 - 10;
    }
    return res;
}
void printArray(double* array, int N) {
    for (int i = 0; i < N; i++) {
        cout << array[i] << ' ';
    }
    cout << endl;
}

double* bubble(double* ar, int N) {
    double* res = new double[N];
    for (int i = 0; i < N; i++) {
        res[i] = ar[i];
    }
    for (int i = 0; i < N - 1; i++)
    {
        int first = i % 2;
        for (int j = first; j < N - 1; j += 2)
        {
            if (res[j] > res[j + 1])
            {
                double prom = res[j];
                res[j] = res[j + 1];
                res[j + 1] = prom;
            }
        }
    }
    return res;
}

double* bubble_parallel(double* ar, int N,int numthreads) {
    double* res = new double[N];
    for (int i = 0; i < N; i++) {
        res[i] = ar[i];
    }
    for (int i = 0; i < N - 1; i++)
    {
        int first = i % 2;
        #pragma omp parallel for collapse(1) num_threads(numthreads) shared(res,first,N)
        for (int j = first; j < N - 1; j += 2)
        {
            if (res[j] > res[j + 1])
            {
                double prom = res[j];
                res[j] = res[j + 1];
                res[j + 1] = prom;
            }
        }
    }
    return res;
}

double*ShellSort(double*arr,int N)
{
    int i, j, step;
    int tmp;
    double* res = new double[N];
    for (int i = 0; i < N; i++) {
        res[i] = arr[i];
    }
    for (step = N / 2; step > 0; step /= 2)
        for (i = step; i < N; i++)
        {
            tmp = res[i];
            for (j = i; j >= step; j -= step)
            {
                if (tmp < res[j - step])
                    res[j] = res[j - step];
                else
                    break;
            }
            res[j] = tmp;
        }
    return res;
}

double* ShellSort_parallel(double* arr, int N, int numthreads)
{
    int i, j, step;
    int tmp;
    double* res = new double[N];
    for (int i = 0; i < N; i++) {
        res[i] = arr[i];
    }
    for (step = N /2; step > 0; step /= 2)

        #pragma omp parallel for collapse(1) num_threads(numthreads),shared(res,step,N)
        for (i = 0; i < step; i++)
        {
            for (int pr = i; pr < N; pr += step) {
                tmp = res[pr];
                for (j = pr; j >= step; j -= step)
                {
                    if (tmp < res[j - step])
                        res[j] = res[j - step];
                    else
                        break;
                }
                res[j] = tmp;
            }
        }
    return res;
}

void qsortRecursive(double*arr, int size){
   
    long i = 0;
    int j = size - 1;

    int mid = arr[size / 2];

    do
    {
        while (arr[i] < mid)
            i++;
        while (arr[j] > mid)
            j--;
        if (i <= j)
        {
            swap(arr[i], arr[j]);
            i++;
            j--;
        }
    } while (i <= j);
        if (j > 0)
            qsortRecursive(arr, j + 1);
        if (size > i)
            qsortRecursive(arr + i, size - i);
}



// Function to swap two numbers a and b
void swap(double* a, double* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

// Function to perform the partitioning
// of array arr[]
int partition(double*arr, int start, int end)
{
    // Declaration
    double pivot = arr[end];
    int i = (start - 1);

    // Rearranging the array
    for (int j = start; j <= end - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[end]);

    // Returning the respective index
    return (i + 1);
}

// Function to perform QuickSort Algorithm
// using openmp
void quicksort(double* arr, int start, int end)
{
    // Declaration
    int index;

    if (start < end) {

        // Getting the index of pivot
        // by partitioning
        index = partition(arr, start, end);

        // Parallel sections
#pragma omp parallel sections
        {
#pragma omp section
            {
                // Evaluating the left half
                quicksort(arr, start, index - 1);
            }
#pragma omp section
            {
                // Evaluating the right half
                quicksort(arr, index + 1, end);
            }
        }
    }
}
int main()
{
    setlocale(LC_ALL, "russian");
    double a,b;
    int N;
    cout << "Количество элементов массива" << endl;
    //cin >> N;

    double start_time = clock();
    double end_time = clock();
    for (int N = 100; N <= 1100; N = N + 100) {
        double* newAr = createArray(N);
        //printArray(newAr, N);
        cout << endl;
        start_time = clock();
        double* res = bubble(newAr, N);
        end_time = clock();
        cout << "Пузырек послед" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = bubble_parallel(newAr, N, 2);
        end_time = clock();
        cout << "Пузырек пар 2" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = bubble_parallel(newAr, N, 3);
        end_time = clock();
        cout << "Пузырек пар 3" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = bubble_parallel(newAr, N, 4);
        end_time = clock();
        cout << "Пузырек пар 4" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = ShellSort(newAr, N);
        end_time = clock();
        cout << "Шел послед" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = ShellSort_parallel(newAr, N, 2);
        end_time = clock();
        cout << "Шел пар 2" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = ShellSort_parallel(newAr, N, 3);
        end_time = clock();
        cout << "Шел пар 3" << end_time - start_time << endl;
        printArray(res, N);

        start_time = clock();
        res = ShellSort_parallel(newAr, N, 4);
        end_time = clock();
        cout << "Шел пар 4" << end_time - start_time << endl;
        printArray(res, N);

        res = new double[N];
        for (int i = 0; i < N; i++) {
            res[i] = newAr[i];
        }
        start_time = clock();
        qsortRecursive(res, N);
        end_time = clock();
        cout << "Быстрая посл" << end_time - start_time << endl;
        printArray(res, N);

        res = new double[N];
        for (int i = 0; i < N; i++) {
            res[i] = newAr[i];
        }
# pragma omp parallel num_threads(2)
        start_time = clock();
        quicksort(res, 0, N);
        end_time = clock();
        cout << "Быстрая пар 2" << end_time - start_time << endl;
        printArray(res, N);

        res = new double[N];
        for (int i = 0; i < N; i++) {
            res[i] = newAr[i];
        }
        start_time = clock();
# pragma omp parallel num_threads(3)
        quicksort(res, 0, N);
        end_time = clock();
        cout << "Быстрая пар 3" << end_time - start_time << endl;
        printArray(res, N);

        res = new double[N];
        for (int i = 0; i < N; i++) {
            res[i] = newAr[i];
        }
# pragma omp parallel num_threads(4)
        start_time = clock();
        quicksort(res, 0, N);
        end_time = clock();
        cout << "Быстрая пар 3" << end_time - start_time << endl;
        printArray(res, N);
        int u;
        cin >> u;
    }
}


