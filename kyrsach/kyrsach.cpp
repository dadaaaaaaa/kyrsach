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
    
}
double betta(int k) {
  
    return -1;
}
double ooo(double x, double y,int k) {
    switch (k)
    {
    case 1:    return 2*exp(x+y);
    case 2:    return -2*exp(x+y);
    }
    return 0;
}
double resh(double x, double y)
{
    return exp(x + y);
}
void tochnoe()
{
    for (int i = 0; i < n; i++)
        x[i] = resh(tch[i][0], tch[i][1]);
}
double gamma(int k) {
    switch (k)
    {
    case 1:    return 3;
    case 2:    return 3;
    }
    return 3;
}
double lambda(int k) {

    switch (k)
    {
    case 1:    return 2;
    case 2:    return 2;
    }
    return 0;
}
double func(double x, double y, int i)
{
    /*switch (i)
    {
    case 1:    return -20;
    case 2:    return 0;
    }*/
    return -exp(x + y);
}
double func_kraev1(double x, double y, int k)
{
    switch (k)
    {
    case 1: return exp(x);
    }
    return 0;
}
double func_kraev3(double* x, int k)
{
    switch (k)
    {
    case 1: return exp(x[0] + 1)-2*exp(x[0]+x[1]);
    }
    return 0;
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
    vrem = new int*[n];
    for (int i = 0; i < n; i++) {
        vrem[i] = new int [n];
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
                        if (vrem[kk-1][l] == tr[i * 4 + j]) {
                            to =1;
                            break;
                        }
                        if (vrem[kk - 1][l] == 0) {
                            break;
                        }
                    }
                    if(to==0)
                    for (int l = 1; l < n; l++) {
                        if (vrem[kk-1][l] == 0) {
                            vrem[kk-1][l] = tr[i * 4 + j];
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
            if(vrem[i - 1][j]>0)
            ig[i] ++;
        }
    }
    int newSize = ig[n];
    jg = new int[newSize];
    jg = copyToOneDimensionalArray(vrem, n, newSize);
    for (int i = 0; i < n + 1; i++)
        cout << ig[i] << " ";
    cout << endl;
        for (int i = 0; i < newSize; i++)
            cout << jg[i] << " ";
}

void M_matrix(double* p1, double* p2, double* p3, double** M_matr, double* local_F, int num_of_area) {
    int i = 0;
    int j = 0;
    double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double mnoz = fabs(det) / 24;
    double* f = new double[3];
    double mnoz2 = gamma(num_of_area)* mnoz;
    f[0] = mnoz * func(p1[0], p1[1], num_of_area);
    f[1] = mnoz * func(p2[0], p2[1], num_of_area) ;
    f[2] = mnoz * func(p3[0], p3[1], num_of_area);
    local_F[0] = 2 * f[0] + f[1] + f[2];
    local_F[1] = f[0] + 2 * f[1] + f[2];
    local_F[2] = f[0] + f[1] + 2 * f[2];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (i == j) {
                M_matr[i][j] = 2 * mnoz;
            }
            else
                M_matr[i][j] = mnoz2;


}

void G_matrix(double*p1, double* p2, double* p3, double** G_matr, int k) {
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
    double fi1 = (p2[0] * p3[1] - p3[0] * p2[1]) + (p2[1]-p3[1]) * p1[0] + (p3[0]-p2[0]) * p1[1];
    double fi2 = (p3[0] * p1[1] - p1[0] * p3[1]) + (p3[1]-p1[1]) * p2[0] + (p1[0]-p3[0]) * p2[1];
    double fi3 = (p1[0] * p2[1] - p1[0] * p2[1]) + (p1[1]-p2[1]) * p3[0] + (p2[0]-p1[0]) * p3[1];
    ck = lambda(k)*fabs(det) / 2;
    double ff = p1[0] * fi1 + p2[0] * fi2 + p3[0] * fi3;

    G_matr[0][0] = ck * (a11*a11+a12*a12);
    G_matr[0][1] = ck * (a11*a21+a12*a22);
    G_matr[0][2] =  ck * (a11*a31+a12*a32);

    G_matr[1][0] = G_matr[0][1];
    G_matr[1][1] = ck * (a21*a21+a22*a22);
    G_matr[1][2] = ck * (a21*a31+a22*a32);

    G_matr[2][0] = G_matr[0][2];
    G_matr[2][1] = G_matr[1][2];
    G_matr[2][2] = ck * (a31*a31+a32*a32);

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
            ggu[i] += value;
            found = true;
            break;
        }
    }

    // Если не нашли в верхней части, проверяем нижнюю часть
    if (!found) {
        for (int i = ig[row]; i < ig[row + 1]; ++i) {
            if (jg[i] - 1 == col) {
                // Если элемент найден, обновляем его значение
                ggl[i] += value;
                found = true;
                break;
            }
        }
    }

    // Если row и col совпадают, обновляем значение на диагонали
    if (row == col) {
                di[row] += value;
    }
}

void zeroOutRow(int row) {
    for (int i = 0; i < ig[n]; ++i) {
        if (jg[i]-1 == row) {
            ggu[i] = 0; // Обнуляем элементы ggu
        }
    }
    int  k=0;
    for (int i = 0; i < n; ++i) {
        if(i==row)
        while (ig[i+1]-( ig[i] +k) != 0) {
            ggl[ig[i]+k] = 0; // Обнуляем элементы ggu
            k++;
        }
    }

}

void pervoe_kraevoe(int vertex1, int vertex2,int form1, int form2) {
    int kol = 0;
    int lbeg;
    int lend;
    // Устанавливаем значения для трех вершин
    di[vertex1] = 1;
    di[vertex2] = 1;

    // Вычисляем значения для F
    F[vertex1] = func_kraev1(tch[vertex1][0], tch[vertex1][1], form1);
    F[vertex2] = func_kraev1(tch[vertex2][0], tch[vertex2][1], form2);
    zeroOutRow(vertex1 );
    zeroOutRow(vertex2 );
}

void tretie(int vertex1, int vertex2, int form1) {
    double** m = new double*[2];
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
    double h = betta(form1)* sqrt(pow((tch[vertex2][0] - tch[vertex1][0]), 2) + pow((tch[vertex2][1] - tch[vertex1][1]), 2))/6;
    F[vertex1] +=h*( m[0][0] * a + m[0][1]* b);
    F[vertex2] +=h*( m[1][0] * a + m[1][1] * b);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            m[i][j] *= h;
        }
    }

    addToElement(vertex1, vertex1,m[0][0]);
    addToElement(vertex1, vertex2,m[0][1]);
    addToElement(vertex2, vertex1,m[1][0]);
    addToElement(vertex2, vertex2,m[1][1]);
}

void vtoroe (int vertex1, int vertex2, int form1) {
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
    F[vertex1]+=h*( m[0][0] * ooo(tch[vertex1][0], tch[vertex1][1],form1) + m[0][1] * ooo(tch[vertex2][0], tch[vertex2][1], form1));
    F[vertex2] +=h*( m[1][0] * ooo(tch[vertex1][0], tch[vertex1][1],form1) + m[1][1] * ooo(tch[vertex2][0], tch[vertex2][1], form1));

}

void local_matrix(int num_of_finit_element, double** local_matr, double* local_F, int io) {
    int ko = tr[io + 0]-1;
    int l = tr[io + 1 ]-1;
    int m = tr[io + 2 ]-1;

    double** M_matr = new double* [3];
    double** G_matr = new double* [3];

    for (int i = 0; i <3; i++) {
        M_matr[i] = new double[3]();
        G_matr[i] = new double[3]();
    }
    M_matrix(tch[ko], tch[l], tch[m], M_matr, local_F, tr[io+3]);
    G_matrix(tch[ko], tch[l], tch[m], G_matr, tr[io + 3]);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            local_matr[i][j] = M_matr[i][j] + G_matr[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        delete[] M_matr[i];
        delete[] G_matr[i];
    }
}

void global_matrix() {
    int i, j, k;
    int* L = new int[3];
    double* local_F = new double[3];
    double** local_matr = new double* [3];
    for (i = 0; i < 3; i++) {
        local_matr[i] = new double[3]();
    }
    int t = 0;

    for (k = 0; k < n2; k++) {
        local_matrix(k, local_matr, local_F, k * 4);
        for (int i = 0; i < 3; i++) {
            L[i] = tr[k * 4 + i] - 1;
        }
        for (int i = 0; i < 3; ++i) {
            int ibeg = L[i]; // Начальный индекс строки
            di[ibeg] += local_matr[i][i]; // Обновление диагонали

            for (int j = i + 1; j < 3; ++j) {
                int iend = L[j]; // Конечный индекс

                // Выбор подходящего диапазона и обновление
                int h;
                if (ibeg < iend) {
                    h = ig[iend];
                    while (jg[h] - 1 < ibeg) h++;
                    ggl[h] += local_matr[i][j]; // Нижний треугольник
                    ggu[h] += local_matr[j][i]; // Верхний треугольник
                }
                else {
                    h = ig[ibeg];
                    while (jg[h] - 1 < iend) h++;
                    ggl[h] += local_matr[i][j]; // Нижний треугольник
                    ggu[h] += local_matr[j][i]; // Верхний треугольник
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

      for (i = 0; i < n5*2; i += 4) 
      tretie(k3[i] - 1, k3[i + 2] - 1, k3[i + 1]);

  
    for (i = 0; i < n3*2 ; i += 4) {
        pervoe_kraevoe(k1[i]-1, k1[i + 2]-1, k1[i + 1], k1[i + 3]);
    }      cout << endl;
      for (int i = 0; i < n; i++) {
          cout << F[i] << " ";
      }
      cout << endl;
      for (int i = 0; i < n; i++) {
          cout << di[i] << " ";
      }
      cout << endl;
      for (int i = 0; i < ig[n]; i++) {
          cout << ggu[i] << " ";
      }
      cout << endl;

      for (int i = 0; i < ig[n]; i++) {
          cout << ggl[i] << " ";

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
    int k, J = jg[kl]-1, p;
    for (k = j; k > 0; k--)
        for (p = ig[J]; p < ig[J + 1]; p++)
            if (jg[p]-1 == jg[kl - k]-1)
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
            gl[kl] = (ggl[kl] - calc(i, j, gl, gu, kl)) / gd[jg[kl]-1];

        for (j = 0; j < kol; j++, ku++)
            gu[ku] = (ggu[ku] - calc(i, j, gu, gl, ku)) / gd[jg[ku]-1];

        gd[i] = sqrt(di[i] - calcD(j, gu, gl, kl));
    }
}

void Mult_A_Vect(double* xn)//перемножение исходной матрицы А на вектор
{
    long i, j, st;
    for (i = 0; i < n; i++) {
        f[i] = di[i] * xn[i];
        for (j = ig[i]; j < ig[i + 1]; j++) {
            st = jg[j]-1;
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
        file11 << setprecision(20) << x[i] << "	" << tr[i]<< " "<< x[i]- tr[i] << endl;
    file11.close();

}
