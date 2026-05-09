#include "tests.h"

#include <iostream>

#include "gauss.h"
#include "LU.h"
#include "utils.h"

// ====================
// Тест 1 — скорость
// ====================
void benchmark_test() {
    std::cout << "\n=== BENCHMARK TEST ===\n";

    int sizes[] = {100, 200, 400};

    for (int n : sizes) {
        Matrix A = generate_matrix(n);
        Vector b = generate_vector(n);

        double t1 = measure_time([&]() {
            gauss_without_pivot(A, b);
        });

        double t2 = measure_time([&]() {
            gauss_with_pivot(A, b);
        });

        Matrix L, U;

        double t3 = measure_time([&]() {
            lu_decomposition(A, L, U);
        });

        double t4 = measure_time([&]() {
            solve_lu(L, U, b);
        });

        std::cout << "n = " << n << "\n";
        std::cout << "Gauss: " << t1 << "\n";
        std::cout << "Gauss pivot: " << t2 << "\n";
        std::cout << "LU decomp: " << t3 << "\n";
        std::cout << "LU solve: " << t4 << "\n";
        std::cout << "----------------------\n";
    }
}

// ====================
// Тест 2 — много правых частей
// ====================
void multi_rhs_test() {
    std::cout << "\n=== MULTI RHS TEST ===\n";

    int n = 300;
    int ks[] = {1, 10, 50};

    Matrix A = generate_matrix(n);

    for (int k : ks) {
        std::vector<Vector> bs(k);

        for (int i = 0; i < k; i++) {
            bs[i] = generate_vector(n);
        }

        double t_gauss = measure_time([&]() {
            for (auto& b : bs) {
                gauss_with_pivot(A, b);
            }
        });

        Matrix L, U;
        lu_decomposition(A, L, U);

        double t_lu = measure_time([&]() {
            for (auto& b : bs) {
                solve_lu(L, U, b);
            }
        });

        std::cout << "k = " << k << "\n";
        std::cout << "Gauss total: " << t_gauss << "\n";
        std::cout << "LU total: " << t_lu << "\n";
        std::cout << "----------------------\n";
    }
}

// ====================
// Тест 3 — устойчивость (Гильберт)
// ====================
void hilbert_test() {
    std::cout << "\n=== HILBERT TEST ===\n";

    int sizes[] = {5, 8, 12};

    for (int n : sizes) {
        Matrix H = generate_hilbert(n);

        Vector x_true(n, 1.0);
        Vector b(n, 0.0);

        // b = H * x_true
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b[i] += H[i][j] * x_true[j];
            }
        }

        Vector x1 = gauss_without_pivot(H, b);
        Vector x2 = gauss_with_pivot(H, b);

        Matrix L, U;
        lu_decomposition(H, L, U);
        Vector x3 = solve_lu(L, U, b);

        std::cout << "n = " << n << "\n";
        std::cout << "||Gauss|| = " << norm(x1) << "\n";
        std::cout << "||Pivot|| = " << norm(x2) << "\n";
        std::cout << "||LU|| = " << norm(x3) << "\n";
        std::cout << "----------------------\n";
    }
}