#include <iostream>
#include <vector>
#include <cmath>

// Определение точки
struct Point {
    double x, y;
};

// Определение треугольника
struct Triangle {
    Point vertices[3];
};

// Граничные условия
enum BoundaryCondition { DIRICHLET, NEUMANN, ROBIN };

// Пример функции правой части уравнения
double f(double x, double y) {
    return -2.0; // Просто пример
}

// Пример граничного условия Дирихле
double g(double x, double y) {
    return x * y; // Просто пример
}

// Пример граничного условия Неймана
double h(double x, double y) {
    return 1.0; // Просто пример
}

// Функция для интегрирования по треугольнику
double integrate(const Triangle& tri, double (*func)(double, double)) {
    double area = 0.5 * std::abs(
        (tri.vertices[0].x * (tri.vertices[1].y - tri.vertices[2].y) +
            tri.vertices[1].x * (tri.vertices[2].y - tri.vertices[0].y) +
            tri.vertices[2].x * (tri.vertices[0].y - tri.vertices[1].y)));

    // Находим центр треугольника
    double center_x = (tri.vertices[0].x + tri.vertices[1].x + tri.vertices[2].x) / 3.0;
    double center_y = (tri.vertices[0].y + tri.vertices[1].y + tri.vertices[2].y) / 3.0;

    return func(center_x, center_y) * area;
}

// Основная функция
int main() {
    // Определяем треугольники
    std::vector<Triangle> triangles = {
        {{{0, 0}, {1, 0}, {0, 1}}},
        {{{1, 0}, {1, 1}, {0, 1}}}
        // Здесь можно добавить больше треугольников
    };

    // Элементы матрицы и вектора правой части
    int n_nodes = 3; // Пример: 2 треугольника с 3 узлами
    std::vector<std::vector<double>> A(n_nodes, std::vector<double>(n_nodes, 0.0));
    std::vector<double> b(n_nodes, 0.0);

    // Проход по всем треугольникам для сборки глобальной матрицы и вектора
    for (const auto& tri : triangles) {
        // Интегрируем матрицу и вектор правой части для каждого элемента
        double area = 0.5 * std::abs(
            (tri.vertices[0].x * (tri.vertices[1].y - tri.vertices[2].y) +
                tri.vertices[1].x * (tri.vertices[2].y - tri.vertices[0].y) +
                tri.vertices[2].x * (tri.vertices[0].y - tri.vertices[1].y)));

        // Пример: упрощенная сборка матрицы A
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                A[i][j] += area / 12.0; // Простая формула для примера
            }
            // Пример: добавление к вектору правой части
            b[i] += integrate(tri, f);
        }
    }

    // Вывод результатов
    std::cout << "Глобальная матрица A:\n";
    for (const auto& row : A) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

    std::cout << "Вектор правой части b:\n";
    for (const auto& val : b) {
        std::cout << val << "\n";
    }

    return 0;
}