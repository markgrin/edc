#include "ressol.hpp"

#include "common.hpp"

#include <iostream>

namespace {
    std::pair<mpz_class, mpz_class> get_q_n_from_prime_minus_1(mpz_class prime) {
        prime = prime - 1;
        mpz_class n = 0;
        while (prime % 2 == 0 && prime != 1) {
            prime /= 2;
            n = n + 1;
        }
        return {prime, n};
    }
    mpz_class pow_ui(mpz_class base, mpz_class exp) {
        if (exp < 1) {
            return 1;
        } else if (exp == 1) {
            return base;
        }
        auto exp_ui = exp.get_ui();
        mpz_class result;
        mpz_pow_ui(result.get_mpz_t(), base.get_mpz_t(), exp_ui);
        return result;
    }
}

mpz_class ressol (mpz_class p, mpz_class a) {
    a = a % p;
    auto [q, n] = get_q_n_from_prime_minus_1(p);
    if (n == 1) {
        return powm(a, (p + 1) / 4, p);
    }
    mpz_class z = 2;
    while (legendre(z, p) != -1) {
        z = z + 1;
    }

    mpz_class c = powm(z, q, p);
    mpz_class r = powm(a, (q + 1) / 2, p);
    mpz_class t = powm(a, q, p);
    mpz_class m = n;
    mpz_class t2 = 0;
    while ((t - 1) % p != 0) {
        t2 = (t * t) % p;
        mpz_class i = 0;
        for (i = 1; i < m; i++){
            if (t2 % p == 1) {
                break;
            }
            t2 = (t2 * t2) % p;
        }
        mpz_class b = powm(c, pow_ui(2, m - i - 1), p);
        r = (r * b) % p;
        c = (b * b) % p;
        t = (t * c) % p;
        m = i;
    }
    return r;
}
