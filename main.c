#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

using namespace std;

const double n = 100;
const double Q = 1.001 - (double)12 / 1000;

double norm(double **a, int N) {
    double ans = 0;
    for (int i = 0; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += abs(a[i][j]);
        }
        if (ans < sum) {
            ans = sum;
        }
    }
    return ans;
}

int main(int argc, char *argv[]) {
    int mode = -1;
    int N;
    sscanf(argv[1], "%d", &mode);
    if (mode == 0) {
        printf("Введите n: ");
        scanf("%d", &N);
    } else if (mode == 1) {
        N = n;
    }
    double **A, **B, **E, **C, **MAS, **NEW, **T;
    double *x, *xs, *y, *b, *b1, *f, *f1, *nevyaz, det = 0.0;
    x = (double *)calloc(N, sizeof(double));
    nevyaz = (double *)calloc(N, sizeof(double));
    xs = (double *)calloc(N, sizeof(double));
    y = (double *)calloc(N, sizeof(double));
    b = (double *)calloc(N, sizeof(double));
    b1 = (double *)calloc(N, sizeof(double));
    f = (double *)calloc(N, sizeof(double));
    f1 = (double *)calloc(N, sizeof(double));
    A = (double **)calloc(N, sizeof(double *));
    B = (double **)calloc(N, sizeof(double *));
    E = (double **)calloc(N, sizeof(double *));
    C = (double **)calloc(N, sizeof(double *));
    MAS = (double **)calloc(N, sizeof(double *));
    NEW = (double **)calloc(N, sizeof(double *));
    T = (double **)calloc(N, sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)calloc(N, sizeof(double));
        B[i] = (double *)calloc(N, sizeof(double));
        E[i] = (double *)calloc(N, sizeof(double));
        C[i] = (double *)calloc(N, sizeof(double));
        MAS[i] = (double *)calloc(N, sizeof(double));
        NEW[i] = (double *)calloc(N, sizeof(double));
        T[i] = (double *)calloc(N, sizeof(double));
    }
    for (int i = 0; i < N; i++) {
        E[i][i] = 1;
    }
    if (mode == 0) {
        FILE *fin = fopen("input.txt", "r");
        // Считывание элементов матрицы
        if (fin) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    fscanf(fin, "%lf", &A[i][j]);
                    B[i][j] = A[i][j];
                    C[i][j] = A[i][j];
                    MAS[i][j] = A[i][j];
                    T[j][i] = A[i][j];
                }
            }
            for (int i = 0; i < N; i++) {
                fscanf(fin, "%lf", &b[i]);
                b1[i] = b[i];
                f[i] = b[i];
            }
        } else {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    scanf("%lf", &A[i][j]);
                    B[i][j] = A[i][j];
                    C[i][j] = A[i][j];
                    MAS[i][j] = A[i][j];
                    T[j][i] = A[i][j];
                }
            }
            for (int i = 0; i < N; i++) {
                scanf("%lf", &b[i]);
                b1[i] = b[i];
                f[i] = b[i];
            }
        }
    } else if (mode == 1) {
        printf("Задайте x для задания вектор столбца b: \n");
        double X;
        scanf("%lf", &X);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    A[i][j] = pow(Q - 1, i + j + 2);
                    B[i][j] = A[i][j];
                    C[i][j] = A[i][j];
                    MAS[i][j] = A[i][j];
                    T[j][i] = A[i][j];
                } else {
                    A[i][j] = pow(Q, i + j + 2) + 0.1 * (j - i);
                    B[i][j] = A[i][j];
                    C[i][j] = A[i][j];
                    MAS[i][j] = A[i][j];
                    T[j][i] = A[i][j];
                }
            }
            b[i] = X * exp(X / (i + 1)) * cos(X / (i + 1));
            b1[i] = b[i];
            f[i] = b[i];
        }
    }
    // Метод Гаусса - сведение к верхнему треуг виду
    for (int q = 0; q < N - 1; q++) {
        if (!B[q][q]) {
            for (int i = q + 1; i < N; i++) {
                if (B[i][q]) {
                    swap(B[q], B[i]);
                    swap(b[q], b[i]);
                    break;
                }
            }
        }
        if (!B[q][q]) {
            printf("\nError! det A = 0\n");
            return 0;
        }
        for (int i = q + 1; i < N; i++) {
            double con = B[i][q] / B[q][q];;
            for (int j = q; j < N; j++) {
                B[i][j] = B[i][j] - B[q][j] * con;
            }
            b[i] = b[i] - b[q] * con;
        }
    }
    int ff = 0;
    // Метод Гаусса с выбором главного элемента
    for (int q = 0; q < N - 1; q++) {
        double maxa = -1;
        double maxx = -1;
        int number = -1;
        for (int i = q; i < N; i++) {
            if (abs(A[i][q]) > maxa) {
                maxa = abs(A[i][q]);
                maxx = A[i][q];
                number = i;
            }
        }
        for (int i = q; i < N; i++) {
            double mi = A[i][q] / maxx;
            if (i != number) {
                for (int j = q; j < N; j++) {
                    A[i][j] = A[i][j] - A[number][j] * mi;
                }
                b1[i] = b1[i] - b1[number] * mi;
            }
        }
        swap(A[number], A[q]);
        ff = (ff + 1) % 2;
        swap(b1[number], b1[q]);
    }
    printf("Матрица после преобразования к верхней треугольной матрице методом Гаусса с выбором главного элемента:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    // Подсчет det A
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            det = A[i][i];
        } else {
            det *= A[i][i];
        }
    }
    if (ff == 1) {
        det *= -1;
    }
    // Подсчет обратной матрицы
    for (int q = 0; q < N - 1; q++) {
        if (!C[q][q]) {
            for (int i = q + 1; i < N; i++) {
                if (C[i][q]) {
                    swap(C[q], C[i]);
                    swap(E[q], E[i]);
                    break;
                }
            }
        }
        for (int i = q + 1; i < N; i++) {
            double con = C[i][q] / C[q][q];
            for (int j = 0; j < N; j++) {
                if (j >= q) {
                    C[i][j] = C[i][j] - C[q][j] * con;
                }
                E[i][j] = E[i][j] - E[q][j] * con;
            }
        }
    }
    double nor = 1;
    for (int i = 0; i < N; i++) {
        nor = C[i][i];
        for (int j = 0; j < N; j++) {
            if (C[i][j] != 0) {
                C[i][j] /= nor;
            }
            if (E[i][j] != 0) {
                E[i][j] /= nor;
            }
        }
    }
    for (int i = N - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double con = C[j][i];
            C[j][i] = 0;
            for (int q = 0; q < N; q++) {
                E[j][q] = E[j][q] - E[i][q] * con;
            }
        }
    }
/*
    // Вывод обратной матрицы
    printf("Обратная матрица:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", E[i][j]);
        }
        printf("\n");
    }

    // Вывод на экран матрицы
    printf("Матрица после преобразования к верхней треугольной матрице:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", B[i][j]);
        }
        printf("%lf\n", b[i]);
    }
    */
    printf("\n\n");
    // Вычисление решений
    for (int i = N - 1; i >= 0; i--) {
        if (i == N - 1) {
            x[i] = b[i] / B[i][i];
        } else {
            double u = b[i];
            for (int j = N - 1; j > i; j--) {
                u = u - x[j] * B[i][j];
            }
            x[i] = u / B[i][i];
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        if (i == N - 1) {
            y[i] = b1[i] / A[i][i];
        } else {
            double u = b1[i];
            for (int j = N - 1; j > i; j--) {
                u = u - y[j] * A[i][j];
            }
            y[i] = u / A[i][i];
        }
    }
    printf("\nОпределитель матрицы:\ndet(A) = %.100lf\n", det);
    int count = 0;
    while (abs(det) < 1) {
        det *= 10;
        count--;
    }
    printf("Определитель матрицы примерно:\ndet(A) = %lf * 10^%d\n", det, count);
    // Подсчет нормы и вывод
    double normaA = norm(A, N);
    printf("Норма матрицы:\n||A|| = %lf\n", normaA);
    double normaAI = norm(E, N);
    printf("Норма обратной матрицы:\n||A^(-1)|| = %lf\n", normaAI);
    double M = normaA * normaAI;
    printf("Число Обусловленности:\nM = %lf\n", M);
    // Вывод решений
    printf("Решение системы методом Гаусса:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %lf\n", i + 1, x[i]);
    }
    printf("Решение системы методом Гаусса с выбором элемента:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %lf\n", i + 1, y[i]);
    }

    //Релаксация
    printf("\n\n\nМетод верхней релаксации:\n\n");
    double w, eps;
    printf("Введите w(0 < w < 2):\n");
    scanf("%lf", &w);
    printf("Введите точность eps:\n");
    scanf("%lf", &eps);
    double discr = INT_MAX;
    int t = 0;
    //Приведение к симметричному виду
    int flag = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (MAS[i][j] != MAS[j][i]) {
                flag = 1;
                break;
            }
        }
    }
    if (flag) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int q = 0; q < N; q++) {
                    NEW[i][j] += T[i][q] * MAS[q][j];
                }
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                f1[i] += T[i][j] * f[j];
            }
        }
    }

    while (discr > eps) {
        discr = 0.0;
        for (int i = 0; i < N; i++) {
            if (NEW[i][i] == 0) {
                printf("Нельзя подсчитать, все элементы на диагонали должны быть положительны\n");
                return 0;
            }
            double constant = w / NEW[i][i];
            double prom = f1[i];
            for (int j = 0; j < N; j++) {
                prom -= NEW[i][j] * xs[j];
            }
            prom *= constant;
            xs[i] += prom;
        }
        for (int i = 0; i < N; i++) {
            double pp = 0.0;
            for (int j = 0; j < N; j++) {
                pp += NEW[i][j] * xs[j];
            }
            nevyaz[i] = pp - f1[i];
            discr += pow(nevyaz[i],2);
        }
        discr = sqrt(discr);
        t++;
    }
    printf("Последовательность сошлась за %d итераций\n", t);
    for (int i = 0; i < N; i++) {
        printf("x%d = %lf\n", i + 1, xs[i]);
    }
    return 0;
}
