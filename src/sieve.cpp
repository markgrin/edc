#include "sieve.hpp"
#include "common.hpp"
#include "ressol.hpp"

#include <vector>
#include <random>

namespace {

/**
 * @brief Заполнение вектора значениями [-1*(sqrt(m) - 1), sqrt(m) - 1]
 * 
 * @param X - пустой вектор
 * @param m - факторизуемое число
 */
void filling_x (std::vector<mpz_class>& X, const mpz_class& m) {

    mpz_class border = sqrt(m) - 1; //  Т.к. gmp не предоставляет функционал
                                    //  для возведения в дробную степень, возьмем 
                                    //  в качестве верхней/нижней границы максимум
                                    //  (его мы можем вычислить средствами gmp)

    X.resize(border.get_ui()*2 + 1);
    for (mpz_class& i : X) {        // Заполняем вектор значениями [-border, border]
        i = border++;
    }
}

/**
 * @brief Функция для генерации квадратичных вычетов
 * 
 * @param x - аргумент функции
 * @param h - константа sqrt(m) с округлением в меньшую сторону
 * @param m - факторизуемое число
 */
void f(mpz_class& x, const mpz_class& h, const mpz_class& m) {
    //f(x) = (x + h)^2 - m
    x += h;
    x *= x;
    x -= m;
}

/**
 * @brief Применение к вектору X функции f(x)
 * 
 * @param X - вектор значений для f(x)
 * @param m - факторизуемое число
 */
void map_fx (std::vector<mpz_class>& X, const mpz_class& m) {
    const mpz_class h = sqrt(m);  //Вычисляем заранее константу h
    for (mpz_class& i : X) {    //Вычиляем для каждого значения f(x)
        f(i, h, m);
    }
}

void init_factor_base (std::vector<mpz_class>& B, const mpz_class& m) {
    mpz_class border = sqrt(m);     //Возьмем верхнюю оценку факторной базы равной sqrt(m)
    B.insert(begin(B), {-1, 2});    //Заполним базу начальными значениями
    mpz_class nxtprime;
    nextprime(nxtprime, B.back());  
    while (nxtprime <= border) {    //Перебираем все простые числа до границы
        if (legendre(m, nxtprime) == 1) {   //Если Cимвол Лежандра == 1, то значение подходит
            B.push_back(nxtprime);
        }
        nextprime(nxtprime, nxtprime);  //Берем следующее простое число
    }
}

/**
 * @brief Увеличивает все элементы с подходящими индексами  на величину ln(p)
 * 
 * @param T - ветор с элементами для увеличения на ln(p)
 * @param x - полученный корень
 * @param p - элемент факторной базы
 * @param border - верхняя/нижняя граница множества X
 */
void distribution_xi(std::vector<mpz_class> T, const mpz_class& x, const mpz_class& p, const mpz_class& border) {
    mpz_class k = 0;
    mpz_class nxt_index = x;
    while (nxt_index <= border && nxt_index < T.size()) {
        //T[nxt_index.get_ui()] += ln(p);   //Нет ln() в gmp
        ++k;
        nxt_index = x + k * p;
    }
    k = -1;
    nxt_index = x - p;
    while (nxt_index >= -border && nxt_index >= 0) {
        //T[nxt_index.get_ui()] += ln(p);   //Нет ln() в gmp
        --k;
        nxt_index = x + k * p;
    }
}

/**
 * @brief Заполняет вектор T значениями по вычисленным индексам
 * 
 * @param T - пустой проициниализированный вектор
 * @param fX - вектор значений f(X)
 * @param B - факторная база
 * @param m - факторизируемое число
 */
void filling_t(
    std::vector<mpz_class> T, 
    const std::vector<mpz_class>& fX, 
    const std::vector<mpz_class>& B, 
    const mpz_class& m) {

    mpz_class border = sqrt(m) - 1; //Снова найдем граници X
    for (const mpz_class& x : fX) {     //Для каждого fxi
        for (const mpz_class& p : B) {  //Для каждого pi
            //auto [x1, x2] = ressol(x, 0, p);  //H^2 == 0 mod p  //Убрать с готовностью ressol
            //distribution_xi(x1, p, border);
            //distribution_xi(x2, p, border);
        }
    }
}

/**
 * @brief Поиск подходящих значений (abs(T[i] - ln(X[i])) <= criterion 
 *        и попытка разложить их на элементы факторной базы
 * 
 * @param T - вектор c результатами просеивания
 * @param B - факторная база
 * @return нетривиальный делитель m, -1 в случае неудачи
 */
mpz_class check_t(
    const std::vector<mpz_class> T, 
    const std::vector<mpz_class>& B) {

        const mpz_class criterion = 0; // Нужно выбрать допустимое отклонение от ln(f(x)) для опробования
        for (size_t i = 0; i < T.size(); ++i) {
            /*if ( abs(T[i] - ln(X[i])) <= criterion) { //Смотрим разницу и если она нас устраивает
                vector<y> res = Функция_разложения_на_множители(T[i], B); //Раскладываем на множители и сверяемся с факторной базой
                if (res.size()) {// Если функция вернула значения y, 
                                    // то значит разложение T[i] входит в
                                    // факторную базу (договоренность о возвращаемом значении)
                    return T[i];
                }
            }
            */
        }
        return -1;
    }
}

mpz_class quadratic_sieve_algorithm (const mpz_class& m) {
    // Алгоритм:
    // 1) Инициализируем множество X
    // 2) Строим вектор f(X)
    // 3) Инициализируем факторную базу
    // 4) Инициализируем пустой вектор T
    // 5) Заполняем вектор Т
    // 6) Ищем Т* и раскладываем каждый элемент Т* на множители, сверяя с факторной базой

    //1) Инициализируем множество X
    std::vector<mpz_class> X;
    filling_x(X, m);

    //2) Строим вектор f(X)
    map_fx(X, m);

    //3) Инициализируем факторную базу
    std::vector<mpz_class> B;
    init_factor_base(B, m);

    //4) Инициализируем пустой вектор T длины |X|
    std::vector<mpz_class> T;
    T.assign(X.size(), 0);  //Перестраховка, альтернатива T(X.size())

    //5) Заполняем вектор Т
    filling_t(T, X, B, m);

    //6) Ищем Т* и раскладываем каждый элемент Т* на множители, сверяя с факторной базой
    return check_t(T, B);

}
