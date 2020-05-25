#pragma once

#include <gmpxx.h>

#include <iostream>
#include <vector>

std::vector<mpz_class> factorization(const mpz_class& f, std::vector<mpz_class>& v);
