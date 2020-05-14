#include "ressol.hpp"

#include "common.hpp"

mpz_class ressol (mpz_class p, mpz_class n, mpz_class q, mpz_class a) {
    if (q % 2 == 0) { // q должно быть нечетным
        return -1;
    }
    if (a % p == 0) { // a должно не делится на p
        return -1;
    }

    if (p % 4 == 3) {
        return powm(a, (p + 1) / 4, p);
    }
    if (p % 8 == 5) {
        mpz_class compare = powm(a, (p - 1) / 4, p);
        if (compare == 1) {
            return powm(a, (p + 3) / 8, p);
        } else if (compare == -1) {
            return 2 * a * powm(4 * a, (p - 5) / 8, p);
        }
    }
    mpz_class c = rand_in_interval(1, p); // Ищем квадратный невычет, случайным перебором
    while (legendre(c, p) != -1) {
        c = rand_in_interval(1, p);
    }

    mpz_class z = powm(c, q, p);
    mpz_class r = n;
    mpz_class t = powm(a, (q - 1) / 2, p);
    mpz_class x = a * t % p;
    mpz_class b = x * t % p;
    while (x % p != 1) {
        // TODO: Дописать алгоритм
        break ;
    }
    return 0;
}