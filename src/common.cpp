#include "common.hpp"

namespace {
    struct random_state__ {
        gmp_randstate_t state;
        random_state__() {
            gmp_randinit_mt(state);
        }
    };
    thread_local random_state__ state;
}

mpz_class powm(mpz_class base, mpz_class exp, mpz_class mod) {
    mpz_class result = 1;
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
    return result;
 }

mpz_class legendre(mpz_class a, mpz_class p) {
    return mpz_legendre(a.get_mpz_t(), p.get_mpz_t());
}

mpz_class rand_in_interval(mpz_class from, mpz_class to) {
    mpz_class result = 1;
    to = to + 1;
    mpz_urandomm(result.get_mpz_t(), state.state, to.get_mpz_t());
    return from + result;
}
