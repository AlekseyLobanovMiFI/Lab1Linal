#include "LU.h"

#include <cmath>
#include <stdexcept>

// ====================
// LU-разложение (без перестановок)
// A = L * U
// ====================
void lu_decomposition(const Matrix& A, Matrix& L, Matrix& U) {
    int n = A.size();

    L = Matrix(n, Vector(n, 0.0));
    U = Matrix(n, Vector(n, 0.0));

    for (int i = 0; i < n; i++) {

        // ===== Верхняя матрица U =====
        for (int j = i; j < n; j++) {
            double sum = 0.0;

            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }

            U[i][j] = A[i][j] - sum;
        }

        // Проверка
        if (std::abs(U[i][i]) < 1e-12) {
            throw std::runtime_error("Zero pivot in LU decomposition");
        }

        // ===== Нижняя матрица L =====
        for (int j = i; j < n; j++) {
            if (i == j) {
                L[i][i] = 1.0;
            } else {
                double sum = 0.0;

                for (int k = 0; k < i; k++) {
                    sum += L[j][k] * U[k][i];
                }

                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }
}


// ====================
// Прямая подстановка: Ly = b
// ====================
Vector forward_substitution(const Matrix& L, const Vector& b) {
    int n = L.size();
    Vector y(n);

    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }

        // т.к. L[i][i] = 1, делить не нужно
        y[i] = b[i] - sum;
    }

    return y;
}


// ====================
// Обратная подстановка: Ux = y
// ====================
Vector backward_substitution(const Matrix& U, const Vector& y) {
    int n = U.size();
    Vector x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }

        if (std::abs(U[i][i]) < 1e-12) {
            throw std::runtime_error("Zero on diagonal in backward substitution");
        }

        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}


// ====================
// Удобная обёртка
// ====================
Vector solve_lu(const Matrix& L, const Matrix& U, const Vector& b) {
    Vector y = forward_substitution(L, b);
    return backward_substitution(U, y);
}