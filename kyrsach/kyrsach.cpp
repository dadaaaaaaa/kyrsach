#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

double* z, * r, * p, * t, * r1, * l, * l1, * f;
double** tch = NULL; // Array for points
double* tr = NULL;  // Array for triangles
int* k1 = NULL;  // Array for k1
int* k2 = NULL;  // Array for k2
int* k3 = NULL; 
double* F = NULL;// Array for k3
int* ig = NULL;
int* jg = NULL;
double* x = NULL;
double* di = NULL;
double* ggl = NULL;
double* ggu = NULL;
int n, n2, n3, n4, n5;


void input() {
    bool kf2 = false, kf3 = false;

    // Read points
    ifstream file("toch.txt");
    file >> n;
    tch = new double*[n]; // Allocate 1D array for points (2 coordinates per tch)
    for (int i = 0; i < n; i++) {
        tch[i] = new double[2];
        for (int j = 0; j < 2; j++) {
            file >> tch[i ][ j]; // Read into the allocated array
            cout << tch[i][ j] << " ";
        }
    }
    file.close();

    // Read triangles
    file.open("tr.txt");
    file >> n2;
    tr = new double[n2 * 7]; // Allocate 1D array for triangles (7 values per triangle)
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < 7; j++) {
            file >> tr[i * 7 + j]; // Read into the allocated array
            if (j == 6 && tr[i * 7 + j] == 3)
                kf3 = true;
            if (j == 6 && tr[i * 7 + j] == 2)
                kf2 = true;
            cout << tr[i * 7 + j] << " ";
        }
    }
    file.close();

    // Read k1
    file.open("kr_1.txt");
    file >> n3;
    k1 = new int[n3 * 2]; // Allocate 1D array for k1 (2 values per entry)
    for (int i = 0; i < n3 * 2; i++) {
        file >> k1[i]; // Read into the allocated array
        cout << k1[i] << " ";
    }
    file.close();

    // Read k2 if k2 is initialized
    if (kf2) {
        file.open("kr_2.txt");
        file >> n4;
        k2 = new int[n4 * 4]; // Allocate 1D array for k2 (4 values per entry)
        for (int i = 0; i < n4 * 4; i++) {
            file >> k2[i]; // Read into the allocated array
            cout << k2[i] << " ";
        }
        file.close();
    }

    // Read k3 if k3 is initialized
    if (kf3) {
        file.open("kr_3.txt");
        file >> n5;
        k3 = new int[n5 * 4]; // Allocate 1D array for k3 (4 values per entry)
        for (int i = 0; i < n5 * 4; i++) {
            file >> k3[i]; // Read into the allocated array
            cout << k3[i] << " ";
        }
        file.close();
    }
}
double resh(double x, double y, int k)
{
    switch (k)
    {
    case 0: return (y * y);
    case 1: return 20 * y - 19;
    }
    return 0;
}
void tochnoe()
{
    for (int i = 0; i < n; i++)
        x[i] = resh(tch[i][0], tch[i][1],0);
}

double lambda(int k, double y, double x) {
    // Пример задания коэффициента диффузии через линейные функции
    return (k == 0) ? 10.0 : 1.0; // Можно настроить на основании координат
}double func(double x, double y, int i)
{
    if (!i) return (-20);
    else return (0);
}
double func_kraev1(double* x, int k)
{
    switch (k)
    {
    case 0: return x[1] * x[1];
    case 1: return x[1];
    }
    return 0;
}
double func_kraev2(double* x, int k)
{
    switch (k)
    {
    case 0: return 20;
    case 1: return 0;
    case 2: return (2.);
    }
    return 0;
}
double func_kraev3(double* x, int k)
{
    switch (k)
    {
    case 0: return (20 * x[1] - 27);
    case 1: return 0;
    }
    return 0;
}

void sort(int* mas, int k)
{
    int l = 0;
    for (int i = k - 1; i >= 0; i--)
        for (int j = 0; j <= i; j++)
            if (mas[j] > mas[j + 1]) {
                l = mas[j];
                mas[j] = mas[j + 1];
                mas[j + 1] = l;
            }
}

void portret()
{
    int i = 0;
    int j = 0;
    int k = 0;
    int kk = 0;
    int key = 0;
    int o=0,pi = 0;
    int position = 0;//позиция в массиве jg, в которую надо добавлять
    int* mas = new int[n];
    struct List {
        int num;
        List* next;
    };

    List* list = new List[n+1];
    List* p = NULL;

    ig = new int[n + 1];

    for (i = 0; i < n+1; i++)
        list[i].next = NULL;

    //составление "массива", содержащего номер точки и смежные с ней 
    for (i = 0; i < n2; i++) {
        for (j = 0; j < 6; j++) { // Теперь 6 узлов на треугольник
            key = 0;
            k = tr[i+o+j];
            
            for (int m = 0; m < 6; m++) {
                if (m != j) {
                    kk = tr[i + pi + m];
                    if (k < kk) {
                        k += kk;
                        kk = k - kk;
                        k -= kk;
                    }
                    p = &list[k];
                    while (p->next) {
                        if (p->next->num == kk) {
                            key = 1;
                            break;
                        }
                        p = p->next;
                    }
                    if (!key) {
                        p->next = new List;
                        p->next->num = kk;
                        p->next->next = NULL;
                    }
                }
            }
        }
o += 6;
pi += 6;
    }

    //составление массива ig
    ig[0] = 0;
    for (i = 0; i < n; i++) {
        k = 0;
        p = &list[i];
        while (p = p->next)
            k++;
        ig[i + 1] = ig[i] + k;
    }

    jg = new int[ig[i] - 1];
    //составление массива jg
    for (i = 0; i < n; i++) {
        k = 0;
        key = 0;
        p = &list[i];
        while (p = p->next) {
            mas[k] = p->num;
            k++;
            key = 1;
        }
        if (key) {
            sort(mas, --k);//сортировка
            int ii = 0;//добавляет в jg
            int jj = 0;
            for (ii = position, jj = 0; ii <= k + position; ii++, jj++)
                jg[ii] = mas[jj];

            position += k + 1;
        }
    }
}
void M_matrix(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6, double** M_matr, double* local_F, int num_of_area) {
    // Вычисляем детерминант для площади треугольника
    double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double mnoz2 = (det / 24);//gammy
    double a = 1.0 / 60, b = 4.0 / 25, c = -1.0 / 360, d = 2.0 / 25;
    double mnoz = fabs(det) / 24;
    double* f = new double[5];
    double tt =  det / 360.0;//na gama domnohit
    M_matr[0][0] = a * mnoz2;
    M_matr[0][1] = c * mnoz2;
    M_matr[0][2] = c * mnoz2;
    M_matr[0][3] = 0;
    M_matr[0][4] = 0;
    M_matr[0][5] = 0;

    M_matr[1][0] = c * mnoz2;
    M_matr[1][1] = a * mnoz2;
    M_matr[1][2] = c * mnoz2;
    M_matr[1][3] = 0;
    M_matr[1][4] = 0;
    M_matr[1][5] = 0;

    M_matr[2][0] = c * mnoz2;
    M_matr[2][1] = c * mnoz2;
    M_matr[2][2] = a * mnoz2;
    M_matr[2][3] = 0;
    M_matr[2][4] = 0;
    M_matr[2][5] = 0;

    M_matr[3][0] = 0;
    M_matr[3][1] = 0;
    M_matr[3][2] = 0;
    M_matr[3][3] = b * mnoz2;
    M_matr[3][4] = d * mnoz2;
    M_matr[3][5] = d * mnoz2;

    M_matr[4][0] = 0;
    M_matr[4][1] = 0;
    M_matr[4][2] = 0;
    M_matr[4][3] = d * mnoz2;
    M_matr[4][4] = b * mnoz2;
    M_matr[4][5] = d * mnoz2;

    M_matr[5][0] = 0;
    M_matr[5][1] = 0;
    M_matr[5][2] = 0;
    M_matr[5][3] = d * mnoz2;
    M_matr[5][4] = d * mnoz2;
    M_matr[5][5] = b * mnoz2;

 


   

    f[0] = mnoz * func(p1[0], p1[1], num_of_area);
    f[1] = mnoz * func(p2[0], p2[1], num_of_area);
    f[2] = mnoz * func(p3[0], p3[1], num_of_area);
    f[3] = mnoz * func(p4[0], p4[1], num_of_area);
    f[4] = mnoz * func(p5[0], p5[1], num_of_area);
    f[5] = mnoz * func(p6[0], p6[1], num_of_area);
    tt = det / 360.0;

    local_F[0] = (6) * f[0] + (-1) * f[1] + (-1) * f[2] + (0) * f[3] + (-4) * f[4] + (0) * f[5];
    local_F[1] = (-1) * f[0] + (6) * f[1] + (-1) * f[2] + (0) * f[3] + (0) * f[4] + (-4) * f[5];
    local_F[2] = (-1) * f[0] + (-1) * f[1] + (6) * f[2] + (-4) * f[3] + (0) * f[4] + (0) * f[5];
    local_F[3] = (0) * f[0] + (0) * f[1] + (-4) * f[2] + (32) * f[3] + (16) * f[4] + (16) * f[5];
    local_F[4] = (-4) * f[0] + (0) * f[1] + (0) * f[2] + (16) * f[3] + (32) * f[4] + (16) * f[5];
    local_F[5] = (0) * f[0] + (-4) * f[1] + (0) * f[2] + (16) * f[3] + (16) * f[4] + (32) * f[5];
    for (int i = 0; i <6; i++)
        local_F[i] *= tt;

    cout << endl;
    cout << endl;

    cout << endl;

}

void G_matrix(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6, double** G_matr, int k) {
    double lamda[6] = { 0,1,2,3,4,5 };
    double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double mn = lamda[1] * det;
    double a11, a21, a12, a22, a31, a32;
    a11 = (p2[1] - p3[1]) ;
    a12 = (p3[0] - p2[0]) ;
    a21 = (p3[1] - p1[1]);
    a22 = (p1[0] - p3[0]) ;
    a31 = (p1[1] - p2[1]) ;
    a32 = (p2[0] - p1[0]) ;
    double h1 = lambda(k, p1[1], p1[0]);
    double h2 = lambda(k, p2[1], p2[0]);
    double h3 = lambda(k, p3[1], p3[0]);
    double tt = 1.0 / 30.0 / det;
    G_matr[0][0] = tt * (a11 * a11 + a12 * a12) * (h1 * (9) + h2 * (3) + h3 * (3));
    G_matr[1][0] = tt * (a21 * a11 + a12 * a22) * (h1 * (-2) + h2 * (-2) + h3 * (-1));
    G_matr[1][1] = tt * (a21 * a21 + a22 * a22) * (h1 * (3) + h2 * (9) + h3 * (3));
    G_matr[2][0] = tt * (a31 * a11 + a12 * a32) * (h1 * (-2) + h2 * (-1) + h3 * (-2));
    G_matr[2][1] = tt * (a31 * a21 + a22 * a32) * (h1 * (-1) + h2 * (-2) + h3 * (-2));
    G_matr[2][2] = tt * (a31 * a31 + a32 * a32) * (h1 * (3) + h2 * (3) + h3 * (9));
        G_matr[3][0] = tt * ((a11 * a11 + a12 * a12) * (h1 * (3) + h2 * (-2) + h3 * (-1)) + (a11 * a21 +
            a12 * a22) * (h1 * (14) + h2 * (3) + h3 * (3)));
    G_matr[3][1] = tt * ((a11 * a21 + a12 * a22) * (h1 * (3) + h2 * (14) + h3 * (3)) + (a21 * a21 + a22 * a22) * (h1 * (-
        2) + h2 * (3) + h3 * (-1)));
    G_matr[3][2] = tt * ((a11 * a31 + a12 * a32) * (h1 * (-1) + h2 * (-2) + h3 * (3)) + (a21 * a31 +
        a22 * a32) * (h1 * (-2) + h2 * (-1) + h3 * (3)));
    G_matr[4][0] = tt * ((a21 * a11 + a22 * a12) * (h1 * (3) + h2 * (-1) + h3 * (-2)) + (a31 * a11 +
        a32 * a12) * (h1 * (3) + h2 * (-2) + h3 * (-1)));
    G_matr[4][1] = tt * ((a21 * a21 + a22 * a22) * (h1 * (-1) + h2 * (3) + h3 * (-2)) + (a31 * a21 +
        a32 * a22) * (h1 * (3) + h2 * (14) + h3 * (3)));
    G_matr[4][2] = tt * ((a21 * a31 + a22 * a32) * (h1 * (3) + h2 * (3) + h3 * (14)) + (a31 * a31 + a32 * a32) * (h1 * (-
        1) + h2 * (-2) + h3 * (3)));
    G_matr[5][0] = tt * ((a11 * a11 + a12 * a12) * (h1 * (3) + h2 * (-1) + h3 * (-2)) + (a31 * a11 +
        a32 * a12) * (h1 * (14) + h2 * (3) + h3 * (3)));
    G_matr[5][1] = tt * ((a11 * a21 + a12 * a22) * (h1 * (-1) + h2 * (3) + h3 * (-2)) + (a31 * a21 +
        a32 * a22) * (h1 * (-2) + h2 * (3) + h3 * (-1)));
    G_matr[5][2] = tt * ((a11 * a31 + a12 * a32) * (h1 * (3) + h2 * (3) + h3 * (14)) + (a31 * a31 + a32 * a32) * (h1 * (-
        2) + h2 * (-1) + h3 * (3)));
    G_matr[3][3] = tt * ((a11 * a11 + a12 * a12) * (h1 * (8) + h2 * (24) + h3 * (8)) + (a11 * a21 +
        a12 * a22) * (h1 * (16) + h2 * (16) + h3 * (8)) + (a21 * a21 + a22 * a22) * (h1 * (24) + h2 * (8) + h3 * (8)));
    G_matr[4][3] = tt * ((a11 * a21 + a12 * a22) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a11 * a31 + a12 * a32) * (h1 * (8)
        + h2 * (24) + h3 * (8)) + (a21 * a21 + a22 * a22) * (h1 * (8) + h2 * (4) + h3 * (8)) + (a21 * a31 + a22 * a32) * (h1 * (8) + h2 * (8) +
            h3 * (4)));
    G_matr[4][4] = tt * ((a21 * a21 + a22 * a22) * (h1 * (8) + h2 * (8) + h3 * (24)) + (a21 * a31 +
        a22 * a32) * (h1 * (8) + h2 * (16) + h3 * (16)) + (a31 * a31 + a32 * a32) * (h1 * (8) + h2 * (24) + h3 * (8)));
    G_matr[5][3] = tt * ((a11 * a11 + a12 * a12) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a11 * a21 + a12 * a22) * (h1 * (8)
        + h2 * (4) + h3 * (8)) + (a11 * a31 + a12 * a32) * (h1 * (8) + h2 * (8) + h3 * (4)) + (a21 * a31 + a22 * a32) * (h1 * (24) + h2 * (8) +
            h3 * (8)));
    G_matr[5][4] = tt * ((a11 * a21 + a12 * a22) * (h1 * (24) + h2 * (8) + h3 * (8)) + (a11 * a31 +
        a12 * a32) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a21 * a31 + a22 * a32) * (h1 * (8) + h2 * (4) + h3 * (8)) + (a31 * a31 +
            a32 * a32) * (h1 * (8) + h2 * (8) + h3 * (4)));
    G_matr[5][5] = tt * ((a11 * a11 + a12 * a12) * (h1 * (8) + h2 * (8) + h3 * (24)) + (a11 * a31 +
        a12 * a32) * (h1 * (16) + h2 * (8) + h3 * (16)) + (a31 * a31 + a32 * a32) * (h1 * (24) + h2 * (8) + h3 * (8)));


}
void pervoe_kraevoe(int current_kraev, int form) {
    int kol = 0, m = 0;
    int lbeg;
    int lend;
    di[current_kraev] = 1;

    F[current_kraev] = func_kraev1(tch[current_kraev], form);

    kol = ig[current_kraev + 1] - ig[current_kraev];

    for (int i = 0; i < kol; i++)
        ggl[ig[current_kraev] + i] = 0;

    for (int i = current_kraev + 1; i < n; i++) {
        lbeg = ig[i];
        lend = ig[i + 1];
        for (int p = lbeg; p < lend; p++)
            if (jg[p] == current_kraev) {
                ggu[p] = 0;
                continue;
            }
    }
}
void local_matrix(int num_of_finit_element, double** local_matr, double* local_F,int io) {
    int ko = tr[io+0];
    int l = tr[io+1];
    int m = tr[io+2];
    int o = tr[io+3];
    int p = tr[io+4];
    int r = tr[io+5];
    // Выделяем матрицы размером 6x6
    double** M_matr = new double* [6];
    double** G_matr = new double* [6];

    for (int i = 0; i < 6; i++) {
        M_matr[i] = new double[6]();
        G_matr[i] = new double[6]();
    }

    // Заполнение матриц
    M_matrix(tch[ko], tch[l], tch[m], tch[o], tch[p], tch[r], M_matr, local_F, io+7);
    G_matrix(tch[ko], tch[l], tch[m], tch[o], tch[p], tch[r], G_matr, io + 7);

    // Объединяем результаты
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            local_matr[i][j] = M_matr[i][j] + G_matr[i][j];
        }
    }

    // Освобождение памяти
    for (int i = 0; i < 6; i++) {
        delete[] M_matr[i];
        delete[] G_matr[i];
    }
    delete[] M_matr;
    delete[] G_matr;
}

void global_matrix() {
    int i, j, k, p = 0, key = 0,o=0;
    int* L = new int[6];
    int* L2 = new int[2];
    int* K = new int[n / 2];

    double* local_F = new double[6];
    double** local_matr = new double* [6];
    for (i = 0; i < 6; i++) {
        local_matr[i] = new double[6]();
    }

    double* b = new double[2]();
    double** a = new double* [2];
    for (i = 0; i < 2; i++) {
        a[i] = new double[2]();
    }

    // Инициализация локальных матриц и правой части
    for (k = 0; k < n2; k++) {
        local_matrix(k, local_matr, local_F,k+o);
        for (int i = 0; i < 6; i++) {
            L[i] = tr[k+o+i];
        }
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++)
                cout << local_matr[i][j] << " ";
            cout << endl;
        }
        // Объединяем локальную матрицу в глобальную
        for (i = 0; i < 6; i++) {
            int ibeg = L[i];
            for (j = i + 1; j < 6; j++) {
                int iend = L[j];
                int h;

                if (ibeg < iend) {
                    h = ig[iend];
                    while (jg[h] < ibeg) h++;
                    ggl[h] += local_matr[i][j];
                    ggu[h] += local_matr[j][i];
                }
                else {
                    h = ig[ibeg];
                    while (jg[h] < iend) h++;
                    ggl[h] += local_matr[i][j];
                    ggu[h] += local_matr[j][i];
                }
            }
            di[ibeg] += local_matr[i][i];
        }

        // Добавляем правую часть
        for (i = 0; i < 6; i++) {
            F[L[i]] += local_F[i];
        }

        for (int i = 0; i < 6; i++) {
            cout << F[L[i]] << " ";
        }
        cout << endl;
        cout << endl;

        cout << endl;
        o += 6;
    }
    for (int i = 0; i < n; i++)
        cout << ggu[i]<<" ";
    cout << endl;
    cout << endl;
    cout << endl;
    for (int i = 0; i < n; i++) {
        cout << ggl[i] << " ";
    }
    cout << endl;
    for (int i = 0; i < n+1; i++)
        cout << ig[i] << " ";
    cout << endl;
    cout << endl;
    cout << endl;
    for (int i = 0; i < n; i++)
        cout << jg[i] << " ";
    cout << endl;
    cout << endl;
    cout << endl;
    for (int i = 0; i < n; i++)
        cout << di[i] << " ";
    cout << endl;
    cout << endl;
    cout << endl;

    // Обработка первого краевого условия
  /*  for (i = 0; i < n3; i++) {
            pervoe_kraevoe(k1[i],k1[i+1]);
    }*/

    // Освобождение памятиa
    delete[] L;
    delete[] L2;
    delete[] local_F;
    for (i = 0; i < 6; i++) {
        delete[] local_matr[i];
    }
    delete[] local_matr;
    for (i = 0; i < 2; i++) {
        delete[] a[i];
    }
    delete[] a;
    delete[] K;
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
    int k, J = jg[kl], p;
    for (k = j; k > 0; k--)
        for (p = ig[J]; p < ig[J + 1]; p++)
            if (jg[p] == jg[kl - k])
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
            gl[kl] = (ggl[kl] - calc(i, j, gl, gu, kl)) / gd[jg[kl]];

        for (j = 0; j < kol; j++, ku++)
            gu[ku] = (ggu[ku] - calc(i, j, gu, gl, ku)) / gd[jg[ku]];

        gd[i] = sqrt(di[i] - calcD(j, gu, gl, kl));
    }
}

void Mult_A_Vect(double* xn)//перемножение исходной матрицы А на вектор
{
    long i, j, st;
    for (i = 0; i < n; i++) {
        f[i] = di[i] * xn[i];
        for (j = ig[i]; j < ig[i + 1]; j++) {
            st = jg[j];
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
    double nvzk, alfa, beta, skp, eps = 9.999999682655226e-030;
    int i;

    double lastnvzk;

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
        skp = sk_pr(p, p);
        alfa = sk_pr(p, r) / skp;
        for (i = 0; i < n; i++) {
            x[i] += alfa * z[i];
            r[i] -= alfa * p[i];
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


    input();    x = new double[n];
    di = new double[n];
    F = new double[n];
    z = new double[n];
    r = new double[n];
    p = new double[n];
    l = new double[n];
    l1 = new double[n];
    f = new double[n];
    tr = new double[n];

    portret();
    cout << endl;
    F = new double[n];
    for (int i = 0; i < n; i++)
        di[i] = F[i] = x[i] = 0;

    ggu = new double[ig[n] - 1];
    ggl = new double[ig[n] - 1];

    for (int i = 0; i < ig[n]; i++)
        ggu[i] = ggl[i] = 0;

    global_matrix();
    LOC();
    for (int i = 0; i < n; i++)
        tr[i] = x[i];

    tochnoe();

    ofstream file11("1.txt");
    for (int i = 0; i < n; i++)
        file11 << setprecision(20) << x[i] << "	" << tr[i] << endl;
    file11.close();

}
