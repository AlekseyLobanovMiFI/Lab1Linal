#include <iostream>

#include "gauss.h"
#include "LU.h"
#include "tests.h"
#include "utils.h"
#include "types.h"

// ====================
// вывод
// ====================
void print_vector(const Vector& v) {
    for (double x : v) std::cout << x << " ";
    std::cout << "\n";
}

void print_matrix(const Matrix& A) {
    for (const auto& row : A) {
        for (double x : row) std::cout << x << " ";
        std::cout << "\n";
    }
}

// ====================
// тесты
// ====================
void run_tests() {
    std::cout << "============= TESTS STARTED ==========\n";

    benchmark_test();
    multi_rhs_test();
    hilbert_test();

    std::cout << "\n============= TESTS ENDED =============\n";
}

// ====================
// меню
// ====================
void menu() {
    std::cout << "\n============= MENU =============\n";
    std::cout << "1. Compare methods\n";
    std::cout << "2. LU decomposition demo\n";
    std::cout << "3. Run tests\n";
    std::cout << "0. Exit\n";
    std::cout << "Choice: ";
}

// ====================
// MAIN
// ====================
int main() {
    bool work = true;

    while (work) {
        menu();

        int choice;
        std::cin >> choice;

        switch (choice) {

            // ====================
            // сравнение методов
            // ====================
            case 1: {
                int n;
                std::cout << "Matrix size: ";
                std::cin >> n;

                Matrix A = generate_matrix(n);
                Vector b = generate_vector(n);

                // КОПИИ (ВАЖНО)
                Matrix A1 = A, A2 = A, A3 = A;
                Vector b1 = b, b2 = b, b3 = b;

                try {
                    Vector x1 = gauss_without_pivot(A1, b1);
                    Vector x2 = gauss_with_pivot(A2, b2);

                    Matrix L, U;
                    lu_decomposition(A3, L, U);
                    Vector x3 = solve_lu(L, U, b3);

                    std::cout << "\nGauss:\n";
                    print_vector(x1);

                    std::cout << "\nGauss with pivot:\n";
                    print_vector(x2);

                    std::cout << "\nLU:\n";
                    print_vector(x3);

                    std::cout << "\nNorms:\n";
                    std::cout << "||x1|| = " << norm(x1) << "\n";
                    std::cout << "||x2|| = " << norm(x2) << "\n";
                    std::cout << "||x3|| = " << norm(x3) << "\n";
                }
                catch (...) {
                    std::cout << "Error during computation\n";
                }

                break;
            }

            // ====================
            // LU демонстрация
            // ====================
            case 2: {
                int n;
                std::cout << "Matrix size: ";
                std::cin >> n;

                Matrix A = generate_matrix(n);
                Matrix L, U;

                try {
                    lu_decomposition(A, L, U);

                    std::cout << "\nL matrix:\n";
                    print_matrix(L);

                    std::cout << "\nU matrix:\n";
                    print_matrix(U);
                }
                catch (...) {
                    std::cout << "LU failed\n";
                }

                break;
            }

            // ====================
            // тесты
            // ====================
            case 3:
                run_tests();
                break;

            // ====================
            // выход
            // ====================
            case 0:
                work = false;
                break;

            default:
                std::cout << "Invalid option\n";
                break;
        }
    }

    return 0;
}