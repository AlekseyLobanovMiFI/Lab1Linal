#pragma once
#include "types.h"
#include <functional>

Matrix generate_matrix(int n);
Vector generate_vector(int n);
Matrix generate_hilbert(int n);

double measure_time(std::function<void()> f);
double norm(const Vector& v);