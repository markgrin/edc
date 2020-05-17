#ifndef HPP_AB06F552D8C44A64962836BB3165C20E

#define HPP_AB06F552D8C44A64962836BB3165C20E

#include <gmpxx.h>

// Файл с полезными функциями

/**
 * Возведение в степень по модулю(с++ обертка)
 * @param base Основание
 * @param exp Показатель
 * @param mod Модуль
 * @return base^exp по модулю mod
 */
 mpz_class powm(mpz_class base, mpz_class exp, mpz_class mod);

 /**
  * Быстрое вычисление символа Лежандра (c++ обертка)
  * @param a не делится на p
  * @param p - простое
  * @return Символ Лежандра (a ; p)
  */
mpz_class legendre(mpz_class a, mpz_class p);

/**
 * Случайное число (с++ обертка)
 * @param from >= 0
 * @param to >= 0
 * @return Случайное число в интервал [from, to]
 */
mpz_class rand_in_interval(mpz_class from, mpz_class to);

/**
 * @brief Найти следующее простое число после a (c++ обертка)
 * 
 * @param a - от какого числа искать ближайшее простое число
 * @param res - куда поместить результат
 */
void nextprime (mpz_class& a, mpz_class& res);

/**
 * @brief Абсолютное значение a
 * 
 * @param a - аргумент, абсолютное значение которого будет вычислено и сохранено в нем же
 */
void abs (mpz_class& a);


#endif // HPP_AB06F552D8C44A64962836BB3165C20E
