#ifndef HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#define HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#include <gmpxx.h>

#include <vector>

/**
 * @brief Этап разложения - Метод Гаусса-Жордана.
 * 
 * @param matrix Матрица с разложением чисел из factors. Она рассматривается транспонированной относительного того,
 * как она вычисляется в этапе решета. 
 * @param width Её ширина
 * @param height Её высота
 * @param factors Сами множители
 * @param lefts Индексы множителей
 * @param h Параметр
 * @param m Раскладываемое число
 * @return mpz_class множитель
 */
mpz_class gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors,
          std::vector<mpz_class>& lefts, mpz_class h, mpz_class m);

#endif // HPP_976D9043F68047C7BC9DC0E1FC46AFEE