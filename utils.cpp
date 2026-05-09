#include "utils.h"

#include <random>
#include <chrono>
#include <cmath>
#include <iostream>

// ====================
// Случайная матрица
// ====================
Matrix generate_matrix(int n) {
    Matrix A(n, Vector(n));

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(-10.0, 10.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = dist(gen);
        }
    }

    return A;
}

// ====================
// Случайный вектор
// ====================
Vector generate_vector(int n) {
    Vector v(n);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(-10.0, 10.0);

    for (int i = 0; i < n; i++) {
        v[i] = dist(gen);
    }

    return v;
}

// ====================
// Матрица Гильберта
// ====================
Matrix generate_hilbert(int n) {
    Matrix H(n, Vector(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }

    return H;
}

// ====================
// Замер времени
// ====================
double measure_time(std::function<void()> f) {
    auto start = std::chrono::high_resolution_clock::now();

    f();

    auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double>(end - start).count();
}

// ====================
// Норма вектора (евклидова)
// ====================
double norm(const Vector& v) {
    double sum = 0.0;

    for (double x : v) {
        sum += x * x;
    }

    return std::sqrt(sum);
}