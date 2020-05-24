#pragma once

#include <iostream>
#include <gmpxx.h>
#include <vector>

using namespace std;

vector<mpz_class> factorization(const mpz_class& f, vector<mpz_class>& v);
