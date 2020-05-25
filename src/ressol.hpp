#ifndef HPP_D05352118F9245B992AAB52F7E7483B2

#define HPP_D05352118F9245B992AAB52F7E7483B2

#include <gmpxx.h>

/**
 * Алгоритм Тоннели-Шенкса, также известный как RESSOL
 *
 * @param p=(2^n) * q + 1
 * @param q нечетное
 * @param a - квадратичный вычет по модулю p
 *
 * @return x такое, что x^2 = a (mod p), -1 в случае ошибки
 */
mpz_class ressol (mpz_class p, mpz_class a);


#endif // HPP_D05352118F9245B992AAB52F7E7483B2
