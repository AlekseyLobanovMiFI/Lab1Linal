#pragma once

#include "types.h"

// LU-разложение матрицы A = L * U
void lu_decomposition(const Matrix& A, Matrix& L, Matrix& U);

// Прямая подстановка: Ly = b
Vector forward_substitution(const Matrix& L, const Vector& b);

// Обратная подстановка: Ux = y
Vector backward_substitution(const Matrix& U, const Vector& y);

// решение через готовые L и U
Vector solve_lu(const Matrix& L, const Matrix& U, const Vector& b);