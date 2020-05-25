#include "sieve.hpp"
#include "common.hpp"
#include "ressol.hpp"

#include <vector>
#include <random>
#include <iostream>

namespace {

/**
 * @brief Функция f(x) = (x + h)^2 - m
 * 
 * @param x - аргумент функции
 * @param h - константа - рекомендуется sqrt(m) с округлением в меньшую сторону
 * @param m - факторизуемое число
 */
mpz_class f(const mpz_class& x, const mpz_class& h, const mpz_class& m) {
    //f(x) = (x + h)^2 - m
    mpz_class result = x + h;
    return result * result - m;
}

mpz_class up_sqrt(const mpz_class& x) {
    mpz_class result = sqrt(x);
    if (result * result == x) {
        return result;
    }
    return result + 1;
}

/**
 * @brief Применение к вектору X функции f(x_i)
 * 
 * @param T - вектор значений для f(x_i)
 * @param m - факторизуемое число
 */
void map_fx (std::vector<mpz_class>& T, const mpz_class& m, const mpz_class& h) {
    mpz_class i = 0;
    for (mpz_class& x : T) {    //Вычиляем для каждого значения f(x)
        x = f(i, h, m);
        i += 1;
    }
}

/**
 * @brief Построение факторной базы - массива простых числел x : legendre(x, p) == 1
 * 
 * @param B - Максимальное значения в фактор базе
 * @param m - факторизуемое число
 * @return std::vector<mpz_class>& 
 */
std::vector<mpz_class> init_factor_base (const mpz_class& B, const mpz_class& m) {
    std::vector<mpz_class> factor_base = {2};
    for (mpz_class i = 3; i <= B; i = nextprime(i)) {
        if (legendre(m, i) == 1) {
            factor_base.push_back(i);
        }
    }
    return factor_base;
}

void apply_solution(std::vector<mpz_class>& sieve, mpz_class x, mpz_class p) {
    for (mpz_class kp = x; kp < sieve.size(); kp = kp + p) {
        // Для всех k : T[x + k * p] = T[x + k * p] + ln(p), вместо ln пока sqrt
        sieve[kp.get_ui()] /= p;
    }
}

mpz_class do_sieve(const std::vector<mpz_class>& factor_base, std::vector<mpz_class>& sieve, mpz_class m, mpz_class h) {
    for (std::size_t i = 0; i < factor_base.size(); i++) {
        auto p = factor_base[i];
        // f(x) = (h + x)^2 - m
        // f(x) = 0 => x = sqrt(m) - h mod p
        mpz_class x1 = ressol(p, m);
        mpz_class x2 = p - x1;
        x1 = (x1 - (h % p) + p) % p;
        x2 = (x2 - (h % p) + p) % p;
        apply_solution(sieve, x1, p);
        if (x1 != x2) {
            apply_solution(sieve, x2, p);
        }
    }
}

std::vector<mpz_class> get_exponents(mpz_class number, const std::vector<mpz_class>& factors) {
    std::vector<mpz_class> result;
    for (const auto& factor : factors) {
        std::size_t i = 0;
        while (number % factor == 0) {
            number /= factor;
            i += 1;
        }
        result.push_back(i % 2);
    }
    return result;
}

}

mpz_class quadratic_sieve_algorithm (const mpz_class& m) {
    const mpz_class B = 29; // Максимальное число в факторной базе
    const std::size_t S = 100; // Размер решета

    //1) Инициализируем Факторную базу
    auto factor_base = init_factor_base(B, m);

    //2) Инициализируем массив решето
    std::vector<mpz_class> sieve(S);

    //3) Заполняем его значениями
    mpz_class h = up_sqrt(m);
    map_fx(sieve, m, h);

    //4) Просеиваем
    do_sieve(factor_base, sieve, m, h);

    //5) Из просеянного решета используем те значения, которые равны 1
    std::vector<mpz_class> left_factors;
    for (std::size_t i = 0; i < sieve.size(); i++) {
        if ( sieve[i] == 1) {
            left_factors.push_back(f(i, h, m));
        }
    }

    std::vector<mpz_class> left_exponents_matrix;
    for (const auto& x : left_factors) {
        auto line = get_exponents(x, factor_base);
        for (auto & y : line) {
        }
        left_exponents_matrix.insert(left_exponents_matrix.end(), line.begin(), line.end());
    }
    return 0;
}
