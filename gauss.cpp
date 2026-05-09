#include "gauss.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// =====================
// Без выбора главного элемента
// =====================
Vector gauss_without_pivot(Matrix A, Vector b) {
    int n = A.size();

    // Прямой ход
    for (int i = 0; i < n; i++) {
        if (std::abs(A[i][i]) < 1e-12) {
            throw std::runtime_error("Zero pivot encountered");
        }

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];

            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }

            b[k] -= factor * b[i];
        }
    }

    // Обратный ход
    Vector x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }

        if (std::abs(A[i][i]) < 1e-12) {
            throw std::runtime_error("Zero on diagonal");
        }

        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

// =====================
// С выбором главного элемента (partial pivoting)
// =====================
Vector gauss_with_pivot(Matrix A, Vector b) {
    int n = A.size();

    // Прямой ход
    for (int i = 0; i < n; i++) {
        // поиск максимума в столбце
        int max_row = i;
        double max_val = std::abs(A[i][i]);

        for (int k = i + 1; k < n; k++) {
            if (std::abs(A[k][i]) > max_val) {
                max_val = std::abs(A[k][i]);
                max_row = k;
            }
        }

        // перестановка строк
        if (max_row != i) {
            std::swap(A[i], A[max_row]);
            std::swap(b[i], b[max_row]);
        }

        if (std::abs(A[i][i]) < 1e-12) {
            throw std::runtime_error("Matrix is singular");
        }

        // зануление
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];

            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }

            b[k] -= factor * b[i];
        }
    }

    // Обратный ход
    Vector x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }

        if (std::abs(A[i][i]) < 1e-12) {
            throw std::runtime_error("Zero on diagonal");
        }

        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}