#pragma once
#include "types.h"

Vector gauss_without_pivot(Matrix A, Vector b);
Vector gauss_with_pivot(Matrix A, Vector b);