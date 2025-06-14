﻿#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <variant>

using namespace std;

double* z, * r, * p, * t, * r1, * l, * l1, * f, vr, tp, t0, t1, * dt0v, * dtv;//vr глобально время сделала чтобы не передавть
double** tch = NULL; // Array for points
double* tr = NULL;  // Array for triangles
double* q1 = NULL;  // Array for triangles
double* q0 = NULL;  // Array for triangles
int* k1 = NULL;  // Array for k1
int* k2 = NULL;  // Array for k2
int* k3 = NULL;
double* T = NULL;
double** M_s3 = NULL;
double* F = NULL;// Array for k3
int* ig = NULL;
int* jg = NULL;
double* x = NULL;
double* x22 = NULL;
 double* x_direct = NULL;

double* di3 = NULL;
double* ggl3 = NULL;
double* ggu3 = NULL;
double* di = NULL;
double* ggl = NULL;
double* ggu = NULL;
double* di2 = NULL;
double* ggl2 = NULL;
double* ggu2 = NULL;
int n, n2, n3, n4, n5, n6;

double u0(double x, double y, int k) {//только для q0  и q1
    switch (k)
    {
        //case 0:   return sin(T[0]);// 
    case 0:   return x+y+T[0]*T[0];
        // case 1:   return sin(T[1]); //
    case 1:   return x+y+T[1]*T[1];
    }
    return 0;
}
double betta(int k) {

    return 2;
}
double ooo(double x, double y, int k) {
 
    return 4*y ;
}
void tochnoe()//u
{
    for (int i = 0; i < n; i++)
        x[i] = tch[i][0]  + tch[i][1]+vr*vr;//lkz 1-без время ,2 
}
double gamma(double x, double y) {//сигма

    return 5;
}
double lambda(double x, double y, int i) {
    //switch (i)
    //{
    //case 1:    return 10;
    //case 2:    return 1;
    //}
    return 2;
}
double func(double x, double y, int i)//f
{
 /*   switch (i)
    {
    case 1:    return -40;
    case 2:    return -4;
    }*/
    return 15*vr*vr;
    //return-5;
    //return 6 * vr - 8;
   // return cos(vr);
}
double func_kraev1(double x, double y, int k)
{
 
    return /*sin(vr);*/x+y+vr*vr*vr;
}
double func_kraev3(double* x, int k)
{
    switch (k)
    {
    case 1: return ( 2*x[0]* x[0] +32 );
    }
    return 0;
}
void input() {


    // Read points
    ifstream file("toch.txt");
    file >> n;
    tch = new double* [n];
    for (int i = 0; i < n; i++) {
        tch[i] = new double[2];
        for (int j = 0; j < 2; j++) {
            file >> tch[i][j];
        }
    }
    file.close();

    // Read triangles
    file.open("tr.txt");
    file >> n2;
    tr = new double[n2 * 4];
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < 4; j++) {
            file >> tr[i * 4 + j];
        }
    }
    file.close();

    // Read k1
    file.open("kr_1.txt");
    file >> n3;
    k1 = new int[n3 * 2];
    for (int i = 0; i < n3 * 2; i++) {
        file >> k1[i];
    }
    file.close();

    // Read k2 if k2 is initialized
    file.open("kr_2.txt");
    file >> n4;
    k2 = new int[n4 * 4];
    for (int i = 0; i < n4 * 4; i++) {
        file >> k2[i];
    }
    file.close();


    file.open("kr_3.txt");
    file >> n5;
    k3 = new int[n5 * 4];
    for (int i = 0; i < n5 * 4; i++) {
        file >> k3[i];
    }
    file.close();
    file.open("t.txt");
    double m, ma;
    file >>m>>ma>> n6;
    double hag = (ma - m) / (n6-1);

    T = new double[n6];
    for (int i = 0; i < n6; i++) {
         T[i]=i*hag;
         }
    file.close();
}


// Функции для LU-разложения (работают с вашими глобальными переменными)
void calcLU() {
    // Используем существующие глобальные массивы
    double* diL = new double[n]; // Диагональ L
    double* L = new double[ig[n]]; // Нижний треугольник L
    double* U = new double[ig[n]]; // Верхний треугольник U

    for (int i = 0; i < n; i++) {
        double sumDi = 0;
        int i0 = ig[i];
        int i1 = ig[i + 1];

        for (int k = i0; k < i1; k++) {
            double suml = 0, sumu = 0;
            int j = jg[k] - 1; // Ваш jg использует 1-based индексацию
            int j0 = ig[j];
            int j1 = ig[j + 1];

            // Поиск совпадающих элементов в строках i и j
            for (int ik = i0, kj = j0; ik < i1 && kj < j1; ) {
                if (jg[ik] > jg[kj])
                    kj++;
                else if (jg[ik] < jg[kj])
                    ik++;
                else {
                    suml += L[ik] * U[kj];
                    sumu += L[kj] * U[ik];
                    ik++;
                    kj++;
                }
            }

            L[k] = (ggl[k] - suml);
            U[k] = (ggu[k] - sumu) / diL[j];
            sumDi += L[k] * U[k];
        }
        diL[i] = di[i] - sumDi;
    }

    // Сохраняем L и U в существующие массивы
    for (int i = 0; i < ig[n]; i++) {
        ggl[i] = L[i];
        ggu[i] = U[i];
    }
    for (int i = 0; i < n; i++) {
        di[i] = diL[i];
    }

    delete[] diL;
    delete[] L;
    delete[] U;
}

void solveY(double* y) {
    // Прямая подстановка Ly = F
    for (int i = 0; i < n; i++) {
        double sum = 0;
        int i0 = ig[i];
        int i1 = ig[i + 1];

        for (int k = i0; k < i1; k++) {
            int j = jg[k] - 1; // 1-based индексация
            sum += ggl[k] * y[j];
        }

        y[i] = (F[i] - sum) / di[i];
    }
}

void solveX(double* y, double* x) {
    // Обратная подстановка Ux = y
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        int i0 = ig[i];
        int i1 = ig[i + 1];

        for (int k = i0; k < i1; k++) {
            int j = jg[k] - 1; // 1-based индексация
            y[j] -= x[i] * ggu[k];
        }
    }
}

void solveLU() {
    // Выделяем память для временных векторов
    double* y = new double[n];

    // Инициализация
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        x_direct[i] = 0.0;
    }

    // LU-разложение
    calcLU();

    // Решение системы
    solveY(y);
    solveX(y, x_direct);

    // Вывод результатов (можно сравнить с решением ЛОС)

    // Освобождаем память
    delete[] y;
}
int* copyToOneDimensionalArray(int** source, int n, int& newSize) {
    int* destination = new int[newSize];
    newSize = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < n; j++) { // Начинаем с 1, чтобы пропустить первый элемент
            if (source[i][j] != 0) {
                destination[newSize] = source[i][j];
                newSize++;
            }
        }
    }

    return destination; // Возвращаем указатель на новый массив
}

void sortArray(int** arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < n - 1; j++) { // Начинаем с 1, чтобы оставить первый элемент
            for (int k = 1; k < n - j; k++) {
                if (arr[i][k] != 0 && arr[i][k + 1] != 0 && arr[i][k] > arr[i][k + 1]) {
                    // Меняем местами
                    int temp = arr[i][k];
                    arr[i][k] = arr[i][k + 1];
                    arr[i][k + 1] = temp;
                }
            }
        }
    }
}
void portret()
{
    int i = 0;
    int j = 0;
    int** vrem;
    ig = new int(n + 1);
    vrem = new int* [n];
    for (int i = 0; i < n; i++) {
        vrem[i] = new int[n];
    }
    for (int i = 0; i < n; i++) {
        vrem[i][0] = (i + 1);
        for (int j = 1; j < n; j++) {
            vrem[i][j] = 0;
        }
    }
    for (int i = 0; i < n + 1; i++)
        ig[i] = 0;
    for (int i = 0; i < n2; i++)
    {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (tr[i * 4 + k] > tr[i * 4 + j]) {
                    int to = 0;
                    int kk = tr[i * 4 + k];
                    for (int l = 1; l < n; l++) {
                        if (vrem[kk - 1][l] == tr[i * 4 + j]) {
                            to = 1;
                            break;
                        }
                        if (vrem[kk - 1][l] == 0) {
                            break;
                        }
                    }
                    if (to == 0)
                        for (int l = 1; l < n; l++) {
                            if (vrem[kk - 1][l] == 0) {
                                vrem[kk - 1][l] = tr[i * 4 + j];
                                break;
                            }

                        }
                }

            }
        }
    }
    sortArray(vrem, n);
    for (int i = 1; i < n + 1; i++) {
        ig[i] = ig[i - 1];
        for (int j = 1; j < n; j++) {
            if (vrem[i - 1][j] > 0)
                ig[i]++;
        }
    }
    int newSize = ig[n];
    jg = new int[newSize];
    jg = copyToOneDimensionalArray(vrem, n, newSize);
}

void M_matrix(double* p1, double* p2, double* p3, double** M_matr, double* local_F, int num_of_area) {
    int i = 0;
    int j = 0;
    double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double mnoz = fabs(det) / 24;
    double* f = new double[3];
    double mnoz2 = gamma(p1[0], p1[1]) * mnoz;
    f[0] = mnoz * func(p1[0], p1[1], num_of_area);
    f[1] = mnoz * func(p2[0], p2[1], num_of_area);
    f[2] = mnoz * func(p3[0], p3[1], num_of_area);
    local_F[0] = 2 * f[0] + f[1] + f[2];
    local_F[1] = f[0] + 2 * f[1] + f[2];
    local_F[2] = f[0] + f[1] + 2 * f[2];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (i == j) {
                M_matr[i][j] = 2 * mnoz2;
            }
            else
                M_matr[i][j] = mnoz2;

}

void G_matrix(double* p1, double* p2, double* p3, double** G_matr, int k) {
    int i = 0;
    int j = 0;
    double ck = 0;
    double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double a11, a12, a21, a22, a31, a32;
    a11 = (p2[1] - p3[1]) / det;//y2-y3/det
    a12 = (p3[0] - p2[0]) / det;//x3-x2/det
    a21 = (p3[1] - p1[1]) / det;//y3-y1/det
    a22 = (p1[0] - p3[0]) / det;//x1-x3/det
    a31 = (p1[1] - p2[1]) / det;//y1-y2/det
    a32 = (p2[0] - p1[0]) / det;//x2-x1/det
    ck = lambda(p1[0], p1[1], k) * fabs(det) / 2;

    G_matr[0][0] = ck * (a11 * a11 + a12 * a12);
    G_matr[0][1] = ck * (a11 * a21 + a12 * a22);
    G_matr[0][2] = ck * (a11 * a31 + a12 * a32);

    G_matr[1][0] = G_matr[0][1];
    G_matr[1][1] = ck * (a21 * a21 + a22 * a22);
    G_matr[1][2] = ck * (a21 * a31 + a22 * a32);

    G_matr[2][0] = G_matr[0][2];
    G_matr[2][1] = G_matr[1][2];
    G_matr[2][2] = ck * (a31 * a31 + a32 * a32);

}


void addToElement(int row, int col, double value) {
    // Проверяем, что row и col находятся в допустимых пределах
    if (row < 0 || col < 0 || row >= n || col >= n) {
        cerr << "Индексы вне диапазона!" << std::endl;
        return;
    }

    // Находим позицию для добавления значения
    bool found = false;

    // Обновление верхней части матрицы (или главной диагонали)
    for (int i = ig[col]; i < ig[col + 1]; ++i) {
        if (jg[i] - 1 == row) {
            // Если элемент найден, обновляем его значение
            ggu3[i] += value;
            found = true;
            break;
        }
    }

    // Если не нашли в верхней части, проверяем нижнюю часть
    if (!found) {
        for (int i = ig[row]; i < ig[row + 1]; ++i) {
            if (jg[i] - 1 == col) {
                // Если элемент найден, обновляем его значение
                ggl3[i] += value;
                found = true;
                break;
            }
        }
    }

    // Если row и col совпадают, обновляем значение на диагонали
    if (row == col) {
        di3[row] += value;
    }
}

void zeroOutRow(int row) {
    for (int i = 0; i < ig[n]; ++i) {
        if (jg[i] - 1 == row) {
            ggu[i] = 0; // Обнуляем элементы ggu
        }
    }
    int  k = 0;
    for (int i = 0; i < n; ++i) {
        if (i == row)
            while (ig[i + 1] - (ig[i] + k) != 0) {
                ggl[ig[i] + k] = 0; // Обнуляем элементы ggu
                k++;
            }
    }

}

void pervoe_kraevoe(int vertex1, int vertex2, int form1, int form2) {
    int kol = 0;
    int lbeg;
    int lend;
    // Устанавливаем значения для трех вершин
    di[vertex1] = 1;
    di[vertex2] = 1;

    // Вычисляем значения для F
    F[vertex1] = func_kraev1(tch[vertex1][0], tch[vertex1][1], form1);
    F[vertex2] = func_kraev1(tch[vertex2][0], tch[vertex2][1], form2);
    zeroOutRow(vertex1);
    zeroOutRow(vertex2);
}

void tretie(int vertex1, int vertex2, int form1) {
    double** m = new double* [2];
    for (int i = 0; i < 2; i++) {
        m[i] = new double[2];
        for (int j = 0; j < 2; j++) {
            if (i == j)
                m[i][j] = 2;
            else
                m[i][j] = 1;
        }
    }
    double  a = func_kraev3(tch[vertex1], form1);
    double b = func_kraev3(tch[vertex2], form1);
    double h = betta(form1) * sqrt(pow((tch[vertex2][0] - tch[vertex1][0]), 2) + pow((tch[vertex2][1] - tch[vertex1][1]), 2)) / 6;
    F[vertex1] += h * (m[0][0] * a + m[0][1] * b);
    F[vertex2] += h * (m[1][0] * a + m[1][1] * b);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            m[i][j] *= h;
        }
    }

    addToElement(vertex1, vertex1, m[0][0]);
    addToElement(vertex1, vertex2, m[0][1]);
    addToElement(vertex2, vertex1, m[1][0]);
    addToElement(vertex2, vertex2, m[1][1]);
}

void vtoroe(int vertex1, int vertex2, int form1) {
    double** m = new double* [2];
    for (int i = 0; i < 2; i++) {
        m[i] = new double[2];
        for (int j = 0; j < 2; j++) {
            if (i == j)
                m[i][j] = 2;
            else
                m[i][j] = 1;
        }
    }
    double h = sqrt(pow((tch[vertex2][0] - tch[vertex1][0]), 2) + pow((tch[vertex2][1] - tch[vertex1][1]), 2)) / 6;
    F[vertex1] += h * (m[0][0] * ooo(tch[vertex1][0], tch[vertex1][1], form1) + m[0][1] * ooo(tch[vertex2][0], tch[vertex2][1], form1));
    F[vertex2] += h * (m[1][0] * ooo(tch[vertex1][0], tch[vertex1][1], form1) + m[1][1] * ooo(tch[vertex2][0], tch[vertex2][1], form1));

}

void local_matrix(int num_of_finit_element, double** local_matr, double** local_matrG, double* local_F, int io) {
    int ko = tr[io + 0] - 1;
    int l = tr[io + 1] - 1;
    int m = tr[io + 2] - 1;

    double** M_matr = new double* [3];
    double** G_matr = new double* [3];

    for (int i = 0; i < 3; i++) {
        M_matr[i] = new double[3]();
        G_matr[i] = new double[3]();
    }
    M_matrix(tch[ko], tch[l], tch[m], M_matr, local_F, tr[io + 3]);
    G_matrix(tch[ko], tch[l], tch[m], G_matr, tr[io + 3]);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            local_matr[i][j] = M_matr[i][j];
            local_matrG[i][j] = G_matr[i][j];
        }
    }
    for (int i = 0; i < 3; i++) {
        delete[] M_matr[i];
        delete[] G_matr[i];
    }
}

void Mult_A_Vectr(double* xn, double* vec)//перемножение исходной матрицы А на вектор
{
    long i, j, st;
    for (i = 0; i < n; i++) {
        vec[i] = di[i] * xn[i];
        for (j = ig[i]; j < ig[i + 1]; j++) {
            st = jg[j] - 1;
            vec[i] += ggl[j] * xn[st];
            vec[st] += ggu[j] * xn[i];
        }
    }
}

void global_matrix() {
    int i, j, k;
    int* L = new int[3];
    double* local_F = new double[3];
    double** local_matr = new double* [3];
    double** local_matrG = new double* [3];

    for (i = 0; i < 3; i++) {
        local_matr[i] = new double[3]();
        local_matrG[i] = new double[3]();
    }
    int t = 0;

    for (k = 0; k < n2; k++) {
        local_matrix(k, local_matr, local_matrG, local_F, k * 4);
        for (int i = 0; i < 3; i++) {
            L[i] = tr[k * 4 + i] - 1;
        }
        for (int i = 0; i < 3; ++i) {
            int ibeg = L[i]; // Начальный индекс строки
            di[ibeg] += local_matr[i][i]; // Обновление диагонали
            di2[ibeg] += local_matrG[i][i];
            for (int j = i + 1; j < 3; ++j) {
                int iend = L[j]; // Конечный индекс

                // Выбор подходящего диапазона и обновление
                int h;
                if (ibeg < iend) {
                    h = ig[iend];
                    while (jg[h] - 1 < ibeg) h++;
                    ggl[h] += local_matr[i][j]; // Нижний треугольник
                    ggu[h] += local_matr[j][i]; // Верхний треугольник
                    ggl2[h] += local_matrG[i][j]; // Нижний треугольник
                    ggu2[h] += local_matrG[j][i]; // Верхний треугольник
                }
                else {
                    h = ig[ibeg];
                    while (jg[h] - 1 < iend) h++;
                    ggl[h] += local_matr[i][j]; // Нижний треугольник
                    ggu[h] += local_matr[j][i]; // Верхний треугольник
                    ggl2[h] += local_matrG[i][j]; // Нижний треугольник
                    ggu2[h] += local_matrG[j][i]; // Верхний треугольник
                }
            }
        }

        // Обновление правой части
        for (int i = 0; i < 3; ++i) {
            F[L[i]] += local_F[i];
        }

    }

    for (i = 0; i < n4 * 2; i += 4)
        vtoroe(k2[i] - 1, k2[i + 2] - 1, k2[i + 1]);


    for (i = 0; i < n5 * 2; i += 4)
        tretie(k3[i] - 1, k3[i + 2] - 1, k3[i + 1]);
    //d
    Mult_A_Vectr(q0, dt0v);
    Mult_A_Vectr(q1, dtv);
    for (i = 0; i < n; i++)
        F[i] = F[i] - t0 * dt0v[i] + t1 * dtv[i];

    //a=tp_M+G//вроде верно
    for (i = 0; i < ig[n]; i++) {
        ggl[i] = (ggl[i] * tp) + ggl2[i];
        ggu[i] = (ggu[i] * tp) + ggu2[i];
    }
    for (i = 0; i < n; i++) {
        di[i] = (di[i] * tp) + di2[i];
    }

    //A+=M_s3
    for (i = 0; i < ig[n]; i ++) {
        ggl[i] += ggl3[i];
        ggu[i] += ggu3[i];
    }
    for (i = 0; i < n; i++)
        di[i] += di3[i];

    for (i = 0; i < n3 * 2; i = i + 4) {
        pervoe_kraevoe(k1[i] - 1, k1[i + 2] - 1, k1[i + 1], k1[i + 3]);
    }
    
}

double norma(double* w)//НОРМА ВЕКТОРА
{
    double s = 0;

    for (int i = 0; i < n; i++)
        s += w[i] * w[i];

    return sqrt(s);
}

double calc(int i, int j, double* gl, double* gu, int kl)//по формуле это сумма которую вычитаем
{
    double s = 0;
    int k, J = jg[kl] - 1, p;
    for (k = j; k > 0; k--)
        for (p = ig[J]; p < ig[J + 1]; p++)
            if (jg[p] - 1 == jg[kl - k] - 1)
                s += gl[kl - k] * gu[p];
    return s;
}

double calcD(int j, double* gl, double* gu, int kl)//аналогично только для диагонали
{
    double s = 0;
    for (int k = kl - j; k < kl; k++)
        s += gl[k] * gu[k];
    return s;
}

void lulu(double* gl, double* gu, double* gd)
{
    int i, j, kol, kl = 0, ku = 0;
    for (i = 0; i < n; i++) {
        kol = ig[i + 1] - ig[i];
        for (j = 0; j < kol; j++, kl++)
            gl[kl] = (ggl[kl] - calc(i, j, gl, gu, kl)) / gd[jg[kl] - 1];

        for (j = 0; j < kol; j++, ku++)
            gu[ku] = (ggu[ku] - calc(i, j, gu, gl, ku)) / gd[jg[ku] - 1];

        gd[i] = sqrt(di[i] - calcD(j, gu, gl, kl));
    }
}

void Mult_A_Vect(double* xn)//перемножение исходной матрицы А на вектор
{
    long i, j, st;
    for (i = 0; i < n; i++) {
        f[i] = di[i] * xn[i];
        for (j = ig[i]; j < ig[i + 1]; j++) {
            st = jg[j] - 1;
            f[i] += ggl[j] * xn[st];
            f[st] += ggu[j] * xn[i];
        }
    }
}//на выходе вектор f который является глобальным

double sk_pr(double* a, double* b)//скалярное произведение векторов.
{
    double s = 0;
    for (int i = 0; i < n; i++)
        s += a[i] * b[i];
    return s;
}

void LOC()
{
    double nvzk, alfa, beta, skp, eps = 9.999999682655226e-014;
    int i;

    double lastnvzk;
    for (int i = 0; i < n; i++) {
        f[i]= 0;
    }
    Mult_A_Vect(x);

    for (i = 0; i < n; i++)
        z[i] = r[i] = F[i] - f[i];

    Mult_A_Vect(z);

    for (i = 0; i < n; i++)
        p[i] = f[i];

    nvzk = sqrt(sk_pr(r, r)) / sqrt(sk_pr(F, F));

    for (int k = 1; k < 10000 && nvzk > eps; k++)
    {
        lastnvzk = nvzk;
        if (k % 10 == 0) {
            for (int i = 0; i < n; i++) {
                f[i] = 0;
            }
            Mult_A_Vect(x);

            for (i = 0; i < n; i++)
                z[i] = r[i] = F[i] - f[i];

            Mult_A_Vect(z);

            for (i = 0; i < n; i++)
                p[i] = f[i];

            nvzk = sqrt(sk_pr(r, r)) / sqrt(sk_pr(F, F));
        }
        skp = sk_pr(p, p);
        alfa = sk_pr(p, r) / skp;
        for (i = 0; i < n; i++) {
            x[i] += alfa * z[i];
            r[i] -= alfa * p[i];
        }
        for (int i = 0; i < n; i++) {
            f[i] = 0;
        }
        Mult_A_Vect(r);
        beta = -sk_pr(p, f) / skp;
        for (i = 0; i < n; i++) {
            z[i] = r[i] + beta * z[i];
            p[i] = f[i] + beta * p[i];
        }
        nvzk = sqrt(sk_pr(r, r)) / sqrt(sk_pr(F, F));
    }
}

int main() {
    double* tr;
    ofstream file11("1.1.txt");
    ofstream error_log("error_log.txt"); // Файл для записи ошибок
    input();
    x = new double[n];
    q0 = new double[n];
    q1 = new double[n];
    di = new double[n];
    di2 = new double[n];
    di3 = new double[n];
    F = new double[n];
    z = new double[n];
    r = new double[n];
    p = new double[n];
    l = new double[n];
    l1 = new double[n];
    f = new double[n];
    dtv = new double[n];
    dt0v = new double[n];
    tr = new double[n];
    x22 = new double[n];
    x_direct = new double[n];
    // Вычисление начальных условий
    for (int i = 0; i < n; i++) {
        q0[i] = u0(tch[i][0], tch[i][1], 0);
        q1[i] = u0(tch[i][0], tch[i][1], 1);
    }

    // Основной цикл по временным слоям
    for (int hi = 2; hi < n6; hi++) {
        vr = T[hi];

        double dt = T[hi] - T[hi - 2];
        double dt0 = T[hi] - T[hi - 1];
        double dt1 = T[hi - 1] - T[hi - 2];
        tp = (dt + dt0) / (dt * dt0);
        t0 = (dt0) / (dt * dt1);
        t1 = (dt) / (dt1 * dt0);

        portret();
        cout << endl;

        // Инициализация массивов
        F = new double[n];
        for (int i = 0; i < n; i++) {
           f[i]= di[i] = di3[i] = di2[i] = F[i] = x[i] =x22[i] = 0;
        }

        // Инициализация матриц
        ggu = new double[ig[n] - 1];
        ggl = new double[ig[n] - 1];
        ggu2 = new double[ig[n] - 1];
        ggl2 = new double[ig[n] - 1];
        ggu3 = new double[ig[n] - 1];
        ggl3 = new double[ig[n] - 1];
        for (int i = 0; i < ig[n]; i++) {
            ggu3[i] = ggl3[i] = ggu2[i] = ggl2[i] = ggu[i] = ggl[i] = 0;
        }

        global_matrix();
        LOC();

        // Сохранение решения
        for (int i = 0; i < n; i++)
            x22[i] = x[i];

        tochnoe();

        // Запись результатов
        file11 << vr << endl;
        for (int i = 0; i < n; i++)
            file11 << setprecision(20) << x[i] << "\t" << x22[i] << "\t" << x[i] - x22[i] << endl;
        // Вычисление погрешностей для последнего слоя

                // Обновление значений для следующего шага
        for (int i = 0; i < n; i++) {
            q0[i] = q1[i];
            q1[i] = x22[i];
        }

    }

    file11.close();
    error_log.close();
    return 0;
}