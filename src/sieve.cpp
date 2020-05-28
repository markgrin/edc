#include "sieve.hpp"
#include "common.hpp"
#include "ressol.hpp"
#include "gauss.hpp"

#include <vector>
#include <random>
#include <iostream>

namespace {

/**
 * @brief Функция f(x) = (x + h)^2 - m
 * 
 * @param x - аргумент функции
 * @param h - константа - рекомендуется sqrt(m) с округлением в большую сторону
 * @param m - факторизуемое число
 */
mpz_class f(const mpz_class& x, const mpz_class& h, const mpz_class& m) {
    //f(x) = (x + h)^2 - m
    mpz_class result = x + h;
    return result * result - m;
}

/**
 * @brief Корень с округлением вверх
 */
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
    std::vector<mpz_class> factor_base = {};
    if (m % 2 == 1) {
        factor_base.push_back(2);
    }
    std::cout << "FACTOR_BASE:" << "\n";
    for (mpz_class i = 3; i <= B; i = nextprime(i)) {
        if (legendre(m, i) == 1) {
            std::cout << i << ", ";
            factor_base.push_back(i);
        }
    }
    std::cout << "\n";
    return factor_base;
}

/**
 * @brief Применяет к решету найденное решение, просеивает
 * 
 * @param sieve Массив с решеткой
 * @param x Начальная позиция
 * @param p Простое число
 */
void apply_solution(std::vector<mpz_class>& sieve, mpz_class x, mpz_class p) {
    for (mpz_class kp = x; kp < sieve.size(); kp = kp + p) {
        sieve[kp.get_ui()] /= p;
    }
}

/**
 * @brief Просеивает решето через факторную базу
 * 
 * @param factor_base факторная база
 * @param sieve решето
 * @param m факторизуемое число
 * @param h параметр
 */
void do_sieve(const std::vector<mpz_class>& factor_base, std::vector<mpz_class>& sieve, mpz_class m, mpz_class h) {
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

/**
 * @brief Раскладывает число на указанные множители
 * 
 * @param number Число
 * @param factors Множетиле
 * @return std::vector<mpz_class> Показатели множетелей
 */
std::vector<mpz_class> get_exponents(mpz_class number, const std::vector<mpz_class>& factors) {
    std::vector<mpz_class> result;
    std::cout << number << " = ";
    for (const auto& factor : factors) {
        std::size_t i = 0;
        while (number % factor == 0) {
            number /= factor;
            i += 1;
        }
        if (i != 0) {
            std::cout << factor << "^" << i << " * ";
        }
        result.push_back(i % 2);
    }
    std::cout << "\n";
    return result;
}

}

mpz_class quadratic_sieve_algorithm (const mpz_class& m) {
    mpz_class B = 500; // Максимальное число в факторной базе
    mpz_class sm = sqrt(m) * 2;
    const std::size_t S = sm.get_ui(); // Размер решета


    if (is_prime(m)) { // Не раскладываем простые числа
        return m;
    }
    if (m % 2 == 0) { // Не раскалдываем четные
        return 2;
    }

    //1) Инициализируем Факторную базу
    auto factor_base = init_factor_base(B, m);

    //2) Инициализируем массив решето
    std::vector<mpz_class> sieve(S);

    //3) Заполняем его значениями
    mpz_class h = up_sqrt(m);
    map_fx(sieve, m, h);
    std::cout << "H = " << h << "\n";

    //4) Просеиваем
    do_sieve(factor_base, sieve, m, h);

    //5) Из просеянного решета используем те значения, которые равны 1
    std::vector<mpz_class> left_factors;
    std::vector<mpz_class> left_xs;
    std::cout << "Factors:\n";
    for (std::size_t i = 0; i < sieve.size(); i++) {
        if ( sieve[i] <= 1) {
            std::cout << "[" << i <<"] => " << f(i, h, m) << "\n";
            left_xs.push_back(i); // Запоминаем индексы - они понадабятся
            left_factors.push_back(f(i, h, m)); // Запоминаем y(x)
        }
    }

    //6) Получаем матрицу
    std::vector<mpz_class> left_exponents_matrix;
    for (const auto& x : left_factors) {
        auto line = get_exponents(x, factor_base); // Раскладываем на множители
        for (auto & y : line) {
        }
        left_exponents_matrix.insert(left_exponents_matrix.end(), line.begin(), line.end());
    }
    //7) Просеивание готово, есть множители, индексы и матрица с разложением этих множителей
    // Все это передаем на следующий этап - метод Гаусса-Жордана
    return gauss(left_exponents_matrix, left_factors.size(), factor_base.size(), left_factors, left_xs, h, m);
}
