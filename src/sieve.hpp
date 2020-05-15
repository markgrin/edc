#pragma once

#include <gmpxx.h>

/**
 * @brief Алгоритм квадратичного решета
 * 
 * @param m - факторизумое число
 * @return нетривиальный делитель m, -1 в случае неудачи
 */
mpz_class quadratic_sieve_algorithm (const mpz_class& m);
