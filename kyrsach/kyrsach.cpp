#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <cmath>
#include <iostream>

typedef double real;
typedef std::vector<real> vr;
typedef std::vector<vr> vvr;
typedef std::vector<int> vi;
typedef std::vector<vi> vvi;

struct Element {
    int region;
    vi _nodes;
};
typedef std::vector<Element> va;

class Problem {
public:
    // Пример методов, которые могут быть в классе Problem
    real lambda(int region, real x, real y) const {
        if (region == 1) return (10);
        else return (1);
    } // лямда в решениии(линейная)
    real gamma(int region) const {
        return 1;
    } // для решения(сохраняеться)
    real rightPart(int region, real x, real y) const {
        if (region) return (-20);
        else return (0);
    } // функция ф(от теста)
    real exactSolution(int region, real x, real y) const {
        switch (region)
        {
        case 0: return (y * y);
        case 1: return 20 * y - 19;
        }
        return 0;
    } //точное решении
    real ug(int region, real x, real y) const {
        switch (region)
        {
        case 0: return y * y;
        case 1: return 20 * y - 19;
        }
        return 0;
    } // первое к улс
    real theta(int region, real x, real y) const {
        switch (region)
        {
        case 0: return 20;
        case 1: return 0;
        }
        return 0;
    } // 2 ур услво├ λ ∂u/∂n┤|_(S_2 )=θ
    real beta(int region) const { return 2; } // для 3 к усл
    real ub(int region, real x, real y) const {
        switch (region)
        {
        case 0: return (20 * y - 27);
        case 1: return 0;
        }
        return 0;
    }//3 кр услов ├ λ ∂u/∂n┤|_(S_3 )+β(├ u┤|_(S_3 )-u_β )=0
};

class Solver {
public:
    Solver();
    void loadElements(const std::string& inputPath = "Elements.txt");
    void loadPoints(const std::string& inputPath = "Points.txt");
    void loadBoundary(const std::string& inputPath = "Boundary.txt");
    void solve();
    void writeAnswer(const std::string& outputPath = "output.txt") const;

private:
    void buildPortrait();
    void buildGlobalMatrix();
    void buildLocal(int i);
    void addLocal(int i);
    void addLocalBoundary(int b);
    void buildLocalStiffnessMatrix(int numberElement);
    void buildLocalMassMatrix(int numberElement);
    void buildLocalRightPart(int numberElement);
    void considerFirstBoundary();
    void considerSecondBoundary();
    void considerThirdBoundary();
    void losInitialize();
    void solveLos();
    void lusqDecomposition();
    void solveL(vr& x, const vr& b);
    void solveU(vr& x, const vr& b);
    void multy(const vr& x, vr& result);
    real dotProduct(const vr& x1, const vr& x2) const;
    real sum(int i, int j);

    Problem _problem;
    vr _xNode;
    vr _yNode;
    va _elements;
    va _mas;
    vvi _ltgm;
    vvi _boundary;
    vvi _sBoundary;
    vvi _ltgb;
    vvr _localMatrix;
    vvr _localMassMatrix;
    vvr _localStiffnessMatrix;
    vr _localRightPart;
    vr _localRightPartBoundary;
    vvr _localMatrixBoundary;
    int num_points = 0; // количество точек
    int countElements = 0; // количество конечных элементов
    int num_elements = 0;
    int num_kraev = 0;

    // LOS's data
    int _N;
    vr _ddi;
    vr _aal;
    vr _aau;
    vr _di;
    vr _al;
    vr _au;
    vr _ia;
    vr _ja;
    vr _x;
    vr _f;
    vr _r;
    vr _r0;
    vr _z;
    vr _p;
    vr _help;
    vr _help2;
    vr _help3;
    vr _M;
    int _maxIter;
    real _eps;
};

const int LOCAL_MATRIX_SIZE = 6;

Solver::Solver() :
    _localMatrix(LOCAL_MATRIX_SIZE, vr(LOCAL_MATRIX_SIZE)),
    _localMassMatrix(_localMatrix),
    _localRightPart(LOCAL_MATRIX_SIZE),
    _localRightPartBoundary(3),
    _localMatrixBoundary(3, vr(3)),
    _eps(1e-13),
    _maxIter(10000) {}

void Solver::loadElements(const std::string& inputPath) {
    std::ifstream input(inputPath.c_str());
    input >> countElements;
    _elements.resize(countElements);
    _mas.resize(countElements);
    _ltgm.resize(countElements);
    for (int i = 0; i < countElements; ++i) {
        _elements[i]._nodes.resize(6);
        _mas[i]._nodes.resize(6);
        _ltgm[i].resize(6);
        input >> _elements[i].region;
        for (int j = 0; j < LOCAL_MATRIX_SIZE; j++) {
            input >> _elements[i]._nodes[j];
            std::cout << _elements[i]._nodes[j];
        }
        _mas[i] = _elements[i];
        std::sort(_mas[i]._nodes.begin(), _mas[i]._nodes.begin() + 6);
        for (int j = 0; j < LOCAL_MATRIX_SIZE; j++) {
            int s;
            for (s = 0; s < LOCAL_MATRIX_SIZE; s++) {
                if (_elements[i]._nodes[s] == _mas[i]._nodes[j]) {
                    _ltgm[i][j] = s;
                    break;
                }
            }
        }
    }
}

void Solver::loadPoints(const std::string& inputPath) {
    std::ifstream input(inputPath.c_str());
    input >> num_points;
    _N = num_points;
    _xNode.resize(num_points);
    _yNode.resize(num_points);
    for (int i = 0; i < num_points; i++) {
        input >> _xNode[i];
        std::cout << _xNode[i];
        input >> _yNode[i];
    }
}

void Solver::loadBoundary(const std::string& inputPath) {
    std::ifstream input(inputPath.c_str());
    int n;
    input >> n;
    _boundary.resize(n, vi(5));
    _sBoundary.resize(n, vi(5));
    _ltgb.resize(n, vi(3));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++) {
            input >> _boundary[i][j];
            std::cout << _boundary[i][j];
        }
        _sBoundary[i] = _boundary[i];
        std::sort(_sBoundary[i].begin() + 2, _sBoundary[i].begin() + 5);
        for (int j = 2; j < 5; j++) {
            int s;
            for (s = 2; s < 5; s++) {
                if (_boundary[i][s] == _sBoundary[i][j]) {
                    _ltgb[i][j - 2] = s - 2;
                    break;
                }
            }
        }
    }
}

void Solver::solve() {
    loadElements();
    loadPoints();
    loadBoundary();
    buildPortrait();
    buildGlobalMatrix();
    solveLos();
    writeAnswer();
}

void Solver::writeAnswer(const std::string& outputPath) const {
    const int FIELD_SIZE = 25;
    std::ofstream output(outputPath.c_str());
    output << "x" << std::setw(FIELD_SIZE)
        << "y" << std::setw(FIELD_SIZE + 1)
        << "u*" << std::setw(FIELD_SIZE - 1)
        << "u" << std::setw(FIELD_SIZE + 7)
        << "|u* - u|" << std::endl;
    output.precision(15);
    output << std::scientific;
    int i = 0;
    real totalDelta = 0;
    real norm = 0;
    for (int j = 0; j < _N; j++) {
        int region = 0;
        real u = _problem.exactSolution(region, _xNode[j], _yNode[j]);
        real y = _x[i];
        real delta = fabs(u - y);
        totalDelta += delta * delta;
        norm += u * u;
        output << _xNode[j] << std::setw(FIELD_SIZE)
            << _yNode[j] << std::setw(FIELD_SIZE)
            << u << std::setw(FIELD_SIZE)
            << y << std::setw(FIELD_SIZE)
            << delta << std::endl;
        ++i;
    }
    output << "Total error: " << sqrt(totalDelta / norm);
}

void Solver::buildPortrait() {
    int position = 0; // позиция в массиве jg, в которую надо добавлять
    std::vector<std::set<int>> list(_N);
    // составление "массива", содержащего номер точки и смежные с ней 
    for (int i = 0; i < countElements; i++) {
        for (int j = 5; j >= 0; j--) {
            for (int m = j - 1; m >= 0; m--) {
                list[_mas[i]._nodes[j]].insert(_mas[i]._nodes[m]);
            }
        }
    }
    // составление массива ia
    _ia.resize(_N + 1, 0);
    for (int i = 1; i <= _N; ++i)
        _ia[i] = _ia[i - 1] + list[i - 1].size();
    _ia[0] = 0;
    _ja.resize(_ia[_N]);
    _al.resize(_ia[_N]);
    _au.resize(_ia[_N]);
    _aal.resize(_ia[_N]);
    _aau.resize(_ia[_N]);
    auto it = _ja.begin();
    for (int i = 0; i < _N; ++i) {
        std::copy(list[i].begin(), list[i].end(), it);
        it += list[i].size();
        list[i].clear();
    }
    _di.resize(_N);
    _ddi.resize(_N);
    _f.resize(_N);
}

void Solver::buildGlobalMatrix() {
    for (int i = 0; i < countElements; i++) {
        buildLocal(i);
        addLocal(i);
    }
    considerSecondBoundary();
    considerThirdBoundary();
    _au = _al;
    considerFirstBoundary();
}

void Solver::buildLocal(int numberElement) {

    buildLocalStiffnessMatrix(numberElement);
    buildLocalMassMatrix(numberElement);
    int i, j;
    for (i = 0; i < LOCAL_MATRIX_SIZE; ++i)
        for (j = 0; j < LOCAL_MATRIX_SIZE; ++j)
            _localMatrix[i][j] += _localMassMatrix[i][j];
    buildLocalRightPart(numberElement);
}

void Solver::addLocal(int numberElement) {
    int index[LOCAL_MATRIX_SIZE];
    vi ll(6);
    for (int k = 0; k < LOCAL_MATRIX_SIZE; k++) {
        index[k] = _mas[numberElement]._nodes[k];
        ll[k] = _ltgm[numberElement][k];
    }
    int i;
    for (i = 0; i < LOCAL_MATRIX_SIZE; ++i) {
        _di[index[i]] += _localMatrix[ll[i]][ll[i]]; // diagonal of global matrix
        _f[index[i]] += _localRightPart[ll[i]];// right part
    }
    // another part of global matrix
    int j;
    for (i = 0; i < LOCAL_MATRIX_SIZE; ++i) {
        int iBeg = _ia[index[i]];
        for (j = 0; j < i; ++j, ++iBeg) {
            int iEnd = _ia[index[i] + 1];
            while (_ja[iBeg] != index[j]) {
                int ind = (iBeg + iEnd) / 2;
                if (_ja[ind] <= index[j])
                    iBeg = ind;
                else
                    iEnd = ind;
            }
            _al[iBeg] += _localMatrix[ll[i]][ll[j]];
            _au[iBeg] += _localMatrix[ll[j]][ll[i]];
        }
    }
}

void Solver::addLocalBoundary(int b) {
    vi lb(3), index(3);
    lb = _ltgb[b];
    int i;
    for (i = 0; i < 3; ++i) {
        index[i] = _sBoundary[b][i + 2];
        _di[index[i]] += _localMatrixBoundary[lb[i]][lb[i]];
    }
    int j;
    for (i = 0; i < 3; ++i) {
        int iBeg = _ia[index[i]];
        for (j = 0; j < i; ++j, ++iBeg) {
            int iEnd = _ia[index[i] + 1];
            while (_ja[iBeg] != index[j]) {
                int ind = (iBeg + iEnd) / 2;
                if (_ja[ind] <= index[j])
                    iBeg = ind;
                else
                    iEnd = ind;
            }
            _al[iBeg] += _localMatrixBoundary[lb[i]][lb[j]];
            _au[iBeg] += _localMatrixBoundary[lb[j]][lb[i]];
        }
    }
}

void Solver::buildLocalStiffnessMatrix(int numberElement) {
    const int reg = _elements[numberElement].region;
    int n1 = _elements[numberElement]._nodes[0];
    int n2 = _elements[numberElement]._nodes[1];
    int n3 = _elements[numberElement]._nodes[2];
    const real det = fabs((_xNode[n2] - _xNode[n1]) * (_yNode[n3] - _yNode[n1]) - (_xNode[n3] - _xNode[n1]) * (_yNode[n2] - _yNode[n1]));
    const real h1 = _problem.lambda(reg, _xNode[n1], _yNode[n1]);
    const real h2 = _problem.lambda(reg, _xNode[n2], _yNode[n2]);
    const real h3 = _problem.lambda(reg, _xNode[n3], _yNode[n3]);
    const real a11 = _yNode[n2] - _yNode[n3];
    const real a12 = _xNode[n3] - _xNode[n2];
    const real a21 = _yNode[n3] - _yNode[n1];
    const real a22 = _xNode[n1] - _xNode[n3];
    const real a31 = _yNode[n1] - _yNode[n2];
    const real a32 = _xNode[n2] - _xNode[n1];
    const real tt = 1.0 / 30.0 / det;
    _localMatrix[0][0] = tt * (a11 * a11 + a12 * a12) * (h1 * (9) + h2 * (3) + h3 * (3));
    _localMatrix[1][0] = tt * (a21 * a11 + a12 * a22) * (h1 * (-2) + h2 * (-2) + h3 * (-1));
    _localMatrix[1][1] = tt * (a21 * a21 + a22 * a22) * (h1 * (3) + h2 * (9) + h3 * (3));
    _localMatrix[2][0] = tt * (a31 * a11 + a12 * a32) * (h1 * (-2) + h2 * (-1) + h3 * (-2));
    _localMatrix[2][1] = tt * (a31 * a21 + a22 * a32) * (h1 * (-1) + h2 * (-2) + h3 * (-2));
    _localMatrix[2][2] = tt * (a31 * a31 + a32 * a32) * (h1 * (3) + h2 * (3) + h3 * (9));
    _localMatrix[3][0] = tt * ((a11 * a11 + a12 * a12) * (h1 * (3) + h2 * (-2) + h3 * (-1)) + (a11 * a21 +
        a12 * a22) * (h1 * (14) + h2 * (3) + h3 * (3)));
    _localMatrix[3][1] = tt * ((a11 * a21 + a12 * a22) * (h1 * (3) + h2 * (14) + h3 * (3)) + (a21 * a21 + a22 * a22) * (h1 * (-
        2) + h2 * (3) + h3 * (-1)));
    _localMatrix[3][2] = tt * ((a11 * a31 + a12 * a32) * (h1 * (-1) + h2 * (-2) + h3 * (3)) + (a21 * a31 +
        a22 * a32) * (h1 * (-2) + h2 * (-1) + h3 * (3)));
    _localMatrix[4][0] = tt * ((a21 * a11 + a22 * a12) * (h1 * (3) + h2 * (-1) + h3 * (-2)) + (a31 * a11 +
        a32 * a12) * (h1 * (3) + h2 * (-2) + h3 * (-1)));
    _localMatrix[4][1] = tt * ((a21 * a21 + a22 * a22) * (h1 * (-1) + h2 * (3) + h3 * (-2)) + (a31 * a21 +
        a32 * a22) * (h1 * (3) + h2 * (14) + h3 * (3)));
    _localMatrix[4][2] = tt * ((a21 * a31 + a22 * a32) * (h1 * (3) + h2 * (3) + h3 * (14)) + (a31 * a31 + a32 * a32) * (h1 * (-
        1) + h2 * (-2) + h3 * (3)));
    _localMatrix[5][0] = tt * ((a11 * a11 + a12 * a12) * (h1 * (3) + h2 * (-1) + h3 * (-2)) + (a31 * a11 +
        a32 * a12) * (h1 * (14) + h2 * (3) + h3 * (3)));
    _localMatrix[5][1] = tt * ((a11 * a21 + a12 * a22) * (h1 * (-1) + h2 * (3) + h3 * (-2)) + (a31 * a21 +
        a32 * a22) * (h1 * (-2) + h2 * (3) + h3 * (-1)));
    _localMatrix[5][2] = tt * ((a11 * a31 + a12 * a32) * (h1 * (3) + h2 * (3) + h3 * (14)) + (a31 * a31 + a32 * a32) * (h1 * (-
        2) + h2 * (-1) + h3 * (3)));
    _localMatrix[3][3] = tt * ((a11 * a11 + a12 * a12) * (h1 * (8) + h2 * (24) + h3 * (8)) + (a11 * a21 +
        a12 * a22) * (h1 * (16) + h2 * (16) + h3 * (8)) + (a21 * a21 + a22 * a22) * (h1 * (24) + h2 * (8) + h3 * (8)));
    _localMatrix[4][3] = tt * ((a11 * a21 + a12 * a22) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a11 * a31 + a12 * a32) * (h1 * (8)
        + h2 * (24) + h3 * (8)) + (a21 * a21 + a22 * a22) * (h1 * (8) + h2 * (4) + h3 * (8)) + (a21 * a31 + a22 * a32) * (h1 * (8) + h2 * (8) +
            h3 * (4)));
    _localMatrix[4][4] = tt * ((a21 * a21 + a22 * a22) * (h1 * (8) + h2 * (8) + h3 * (24)) + (a21 * a31 +
        a22 * a32) * (h1 * (8) + h2 * (16) + h3 * (16)) + (a31 * a31 + a32 * a32) * (h1 * (8) + h2 * (24) + h3 * (8)));
    _localMatrix[5][3] = tt * ((a11 * a11 + a12 * a12) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a11 * a21 + a12 * a22) * (h1 * (8)
        + h2 * (4) + h3 * (8)) + (a11 * a31 + a12 * a32) * (h1 * (8) + h2 * (8) + h3 * (4)) + (a21 * a31 + a22 * a32) * (h1 * (24) + h2 * (8) +
            h3 * (8)));
    _localMatrix[5][4] = tt * ((a11 * a21 + a12 * a22) * (h1 * (24) + h2 * (8) + h3 * (8)) + (a11 * a31 +
        a12 * a32) * (h1 * (4) + h2 * (8) + h3 * (8)) + (a21 * a31 + a22 * a32) * (h1 * (8) + h2 * (4) + h3 * (8)) + (a31 * a31 +
            a32 * a32) * (h1 * (8) + h2 * (8) + h3 * (4)));
    _localMatrix[5][5] = tt * ((a11 * a11 + a12 * a12) * (h1 * (8) + h2 * (8) + h3 * (24)) + (a11 * a31 +
        a12 * a32) * (h1 * (16) + h2 * (8) + h3 * (16)) + (a31 * a31 + a32 * a32) * (h1 * (24) + h2 * (8) + h3 * (8)));
}

void Solver::buildLocalMassMatrix(int numberElement) {
    const int reg = _elements[numberElement].region;
    int n1 = _elements[numberElement]._nodes[0];
    int n2 = _elements[numberElement]._nodes[1];
    int n3 = _elements[numberElement]._nodes[2];
    const real det = fabs((_xNode[n2] - _xNode[n1]) * (_yNode[n3] - _yNode[n1]) - (_xNode[n3] -
        _xNode[n1]) * (_yNode[n2] - _yNode[n1]));
    const real g = _problem.gamma(reg);
    const real tt = g * det / 360.0;
    _localMassMatrix[0][0] = tt * (6);
    _localMassMatrix[1][0] = tt * (-1);
    _localMassMatrix[1][1] = tt * (6);
    _localMassMatrix[2][0] = tt * (-1);
    _localMassMatrix[2][1] = tt * (-1);
    _localMassMatrix[2][2] = tt * (6);
    _localMassMatrix[3][0] = 0;
    _localMassMatrix[3][1] = 0;
    _localMassMatrix[3][2] = tt * (-4);
    _localMassMatrix[4][0] = tt * (-4);
    _localMassMatrix[4][1] = 0;
    _localMassMatrix[4][2] = 0;
    _localMassMatrix[5][0] = 0;
    _localMassMatrix[5][1] = tt * (-4);
    _localMassMatrix[5][2] = 0;
    _localMassMatrix[3][3] = tt * (32);
    _localMassMatrix[4][3] = tt * (16);
    _localMassMatrix[4][4] = tt * (32);
    _localMassMatrix[5][3] = tt * (16);
    _localMassMatrix[5][4] = tt * (16);
    _localMassMatrix[5][5] = tt * (32);
}

void Solver::buildLocalRightPart(int numberElement) {
    int n1 = _elements[numberElement]._nodes[0];
    int n2 = _elements[numberElement]._nodes[1];
    int n3 = _elements[numberElement]._nodes[2];
    int n4 = _elements[numberElement]._nodes[3];
    int n5 = _elements[numberElement]._nodes[4];
    int n6 = _elements[numberElement]._nodes[5];
    const real det = fabs((_xNode[n2] - _xNode[n1]) * (_yNode[n3] - _yNode[n1]) - (_xNode[n3] -
        _xNode[n1]) * (_yNode[n2] - _yNode[n1]));
    const int reg = _elements[numberElement].region;
    const real tt = det / 360.0;
    const real f1 = _problem.rightPart(reg, _xNode[n1], _yNode[n1]);
    const real f2 = _problem.rightPart(reg, _xNode[n2], _yNode[n2]);
    const real f3 = _problem.rightPart(reg, _xNode[n3], _yNode[n3]);
    const real f4 = _problem.rightPart(reg, _xNode[n4], _yNode[n4]);
    const real f5 = _problem.rightPart(reg, _xNode[n5], _yNode[n5]);
    const real f6 = _problem.rightPart(reg, _xNode[n6], _yNode[n6]);
    _localRightPart[0] = (6) * f1 + (-1) * f2 + (-1) * f3 + (0) * f4 + (-4) * f5 + (0) * f6;
    _localRightPart[1] = (-1) * f1 + (6) * f2 + (-1) * f3 + (0) * f4 + (0) * f5 + (-4) * f6;
    _localRightPart[2] = (-1) * f1 + (-1) * f2 + (6) * f3 + (-4) * f4 + (0) * f5 + (0) * f6;
    _localRightPart[3] = (0) * f1 + (0) * f2 + (-4) * f3 + (32) * f4 + (16) * f5 + (16) * f6;
    _localRightPart[4] = (-4) * f1 + (0) * f2 + (0) * f3 + (16) * f4 + (32) * f5 + (16) * f6;
    _localRightPart[5] = (0) * f1 + (-4) * f2 + (0) * f3 + (16) * f4 + (16) * f5 + (32) * f6;
    for (int i = 0; i < LOCAL_MATRIX_SIZE; i++)
        _localRightPart[i] *= tt;
}

void Solver::considerFirstBoundary() {
    int i, n;
    real ug;
    for (i = 0, n = _boundary.size(); i < n; ++i) {
        if (_boundary[i][0] == 1) break;
    }
    const double BIG_NUM = 1e50;
    while (i < n && _boundary[i][0] == 1) {
        int region = _boundary[i][1];
        for (int j = 2; j < 5; j++) {
            ug = _problem.ug(region, _xNode[_boundary[i][j]], _yNode[_boundary[i][j]]);
            _di[_boundary[i][j]] = BIG_NUM;
            _f[_boundary[i][j]] = ug * BIG_NUM;
        }
        ++i;
    }
}

void Solver::considerSecondBoundary() {
    int i, n = _boundary.size(), region;
    int n1, n2, n3;
    real t1, t2, t3;
    real tt;
    for (i = 0; i < n; ++i) {
        if (_boundary[i][0] == 2)
            break;
    }
    while (i < n && _boundary[i][0] == 2)
    {
        region = _boundary[i][1];
        n1 = _boundary[i][2];
        n2 = _boundary[i][3];
        n3 = _boundary[i][4];
        t1 = _problem.theta(region, _xNode[n1], _yNode[n1]);
        t2 = _problem.theta(region, _xNode[n2], _yNode[n2]);
        t3 = _problem.theta(region, _xNode[n3], _yNode[n3]);
        tt = sqrt((_xNode[n1] - _xNode[n3]) * (_xNode[n1] - _xNode[n3]) +
            (_yNode[n1] - _yNode[n3]) * (_yNode[n1] - _yNode[n3])) / 30.0;
        _localRightPartBoundary[0] = tt * ((4) * t1 + (2) * t2 + (-1) * t3);
        _localRightPartBoundary[1] = tt * ((2) * t1 + (16) * t2 + (2) * t3);
        _localRightPartBoundary[2] = tt * ((-1) * t1 + (2) * t2 + (4) * t3);
        _f[n1] += _localRightPartBoundary[0];
        _f[n2] += _localRightPartBoundary[1];
        _f[n3] += _localRightPartBoundary[2];
        ++i;
    }

}

void Solver::considerThirdBoundary() {
    int i, n = _boundary.size(), region;
    int n1, n2, n3;
    real ub1, ub2, ub3;
    real tt;
    real beta;
    vi index(3);
    for (i = 0; i < n; ++i) {
        if (_boundary[i][0] == 3)
            break;
    }
    while (i < n && _boundary[i][0] == 3)
    {
        region = _boundary[i][1];
        beta = _problem.beta(region);
        n1 = _boundary[i][2];
        n2 = _boundary[i][3];
        n3 = _boundary[i][4];
        ub1 = _problem.ub(region, _xNode[n1], _yNode[n1]);
        ub2 = _problem.ub(region, _xNode[n2], _yNode[n2]);
        ub3 = _problem.ub(region, _xNode[n3], _yNode[n3]);
        tt = beta * sqrt((_xNode[n1] - _xNode[n3]) * (_xNode[n1] - _xNode[n3]) +
            (_yNode[n1] - _yNode[n3]) * (_yNode[n1] - _yNode[n3])) / 30.0;
        _localRightPartBoundary[0] = tt * ((4) * ub1 + (2) * ub2 + (-1) * ub3);
        _localRightPartBoundary[1] = tt * ((2) * ub1 + (16) * ub2 + (2) * ub3);
        _localRightPartBoundary[2] = tt * ((-1) * ub1 + (2) * ub2 + (4) * ub3);
        _f[n1] += _localRightPartBoundary[0];
        _f[n2] += _localRightPartBoundary[1];
        _f[n3] += _localRightPartBoundary[2];
        _localMatrixBoundary[0][0] = tt * 4;
        _localMatrixBoundary[0][1] = tt * 2;
        _localMatrixBoundary[0][2] = tt * (-1);
        _localMatrixBoundary[1][0] = tt * 2;
        _localMatrixBoundary[1][1] = tt * 16;
        _localMatrixBoundary[1][2] = tt * 2;
        _localMatrixBoundary[2][0] = tt * (-1);
        _localMatrixBoundary[2][1] = tt * 2;
        _localMatrixBoundary[2][2] = tt * 4;
        index[0] = n1;
        index[1] = n2;
        index[2] = n3;
        addLocalBoundary(i);
        ++i;
    }

}

void Solver::losInitialize() {
    _x.resize(_N);
    _r0.resize(_N);
    _r.resize(_N);
    _z.resize(_N);
    _p.resize(_N);
    _help.resize(_N);
    _help2.resize(_N);
    _help3.resize(_N);
}

void Solver::solveLos() {
    losInitialize();
    lusqDecomposition();
    multy(_x, _help);
    for (int i = 0; i < _N; ++i)
        _r[i] = _f[i] - _help[i];
    solveL(_r0, _r);
    _r = _r0;
    solveU(_z, _r);
    multy(_z, _help);
    solveL(_p, _help);
    real rr = 1;
    int k = 0;
    _eps *= _eps;
    while (rr > _eps && k < _maxIter)
    {
        real pp = dotProduct(_p, _p);
        double alpha = dotProduct(_p, _r) / pp;
        for (int i = 0; i < _N; ++i) {
            _x[i] += alpha * _z[i];
            _r[i] -= alpha * _p[i];
        }
        solveU(_help, _r);
        multy(_help, _help2);
        solveL(_help3, _help2);
        real beta = -dotProduct(_p, _help3) / pp;
        for (int i = 0; i < _N; ++i) {
            _z[i] = _help[i] + beta * _z[i];
            _p[i] = _help3[i] + beta * _p[i];
        }
        if (k % 100 != 0)
            rr -= alpha * alpha * pp;
        else
            rr = dotProduct(_r, _r);
        ++k;
    }
}

void Solver::lusqDecomposition() {
    _aal = _aau = _al;
    _ddi = _di;
    int j;
    for (int i = 0; i < _N; ++i) {
        for (j = _ia[i]; j < _ia[i + 1]; ++j)
        {
            _al[j] = (_aal[j] - sum(i, _ja[j]));
            _au[j] = (_aau[j] - sum(_ja[j], i)) / _di[_ja[j]];
        }
        _di[i] = _ddi[i] - sum(i, i);
    }
}

void Solver::solveL(vr& x, const vr& b) {
    int j;
    for (int i = 0; i < _N; ++i) {
        double s = b[i];
        int jEnd = _ia[i + 1];
        for (j = _ia[i]; j < jEnd; ++j)
            s -= _al[j] * x[_ja[j]];
        x[i] = s / _di[i];
    }
}

void Solver::solveU(vr& x, const vr& b) {
    x = b;
    int j;
    for (int i = _N - 1; i > -1; --i) {
        double xi = x[i];
        int jEnd = _ia[i + 1];
        for (j = _ia[i]; j < jEnd; ++j)
            x[_ja[j]] -= _au[j] * xi;
    }
}

void Solver::multy(const vr& x, vr& result) {
    int kk, i, j, kol_str, st;
    kk = 0;
    for (i = 0; i < _N; ++i)
        result[i] = _ddi[i] * x[i];
    for (i = 1; i < _N; ++i) {
        kol_str = _ia[i + 1] - _ia[i];
        for (j = 0; j < kol_str; ++j, ++kk)
        {
            st = _ja[kk];
            result[i] += _aal[kk] * x[st];
            result[st] += _aau[kk] * x[i];
        }
    }
}

real Solver::dotProduct(const vr& x1, const vr& x2) const {
    real s = 0;
    for (int i = 0, n = x1.size(); i < n; ++i)
        s += x1[i] * x2[i];
    return s;

}

real Solver::sum(int i, int j) {
    int k, l;
    bool find;
    real result = 0;
    if (i == j) {
        for (k = _ia[i]; k < _ia[i + 1]; ++k)
            result += _au[k] * _al[k];
    }
    else if (i > j) {
        // верхний треугольник
        for (k = _ia[j]; k < _ia[j + 1]; k++) {
            find = false;
            for (l = _ia[i]; l < _ia[i + 1]
                && !find; ++l) {
                if (_ja[l] == _ja[k]
                    ) {
                    result +=
                        _au[k] * _al[l];
                    find = true;
                }
            }
        }
    }
    else {
        // нижний треугольник
        for (l = _ia[i]; l < _ia[i + 1]; l++)
        {
            find = false;
            for (k = _ia[j]; k < _ia[j + 1]
                && !find; ++k) {
                if (_ja[l] == _ja[k]
                    ) {
                    result +=
                        _au[k] * _al[l];
                    find = true;
                }
            }
        }
    }
    return result;

}

int main() {
    Solver solver;
    solver.solve();
    return 0;
}
