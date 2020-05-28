#include "gauss.hpp"

#include <iostream>
#include <utility>

namespace {

/**
 * Функция для получения числа по индексу. Написана так, чтобы выдавать числа из транспонированной матрицы.
 */
mpz_class& get_num(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t x, std::size_t y) {
    return matrix[x * height + y];
}

/**
 * Вычитает из i строки j
 */
void sub_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j) {
    for (std::size_t x = 0; x < width; x++) {
        mpz_class& num = get_num(matrix, width, height, x, i);
        num = (num + get_num(matrix, width, height, x, j)) % 2;
    }
}

/**
 * Меняет i, j строки
 */
void swap_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j) {
    for (std::size_t x = 0; x < width; x++) {
        std::swap(get_num(matrix, width, height, x, i), get_num(matrix, width, height, x, j));
    }
}

/**
 * В столбце i находит pivot элемент начиная с clear_y позиции
 */
std::size_t find_1_index(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t clear_y) {
    for (std::size_t y = clear_y; y < height; y++) {
        if (get_num(matrix, width, height, i, y)) {
            return y;
        }
    }
    return height + 1;
}

/**
 * Печать матрицы на экран. Так как матрица mod 2 - получается красиво.
 */
void print(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height) {
    for (std::size_t y = 0; y < height; y++) {
        if (y <= 9) {
            std::cout << " ";
        }
        std::cout << y << "]  ";
        for (std::size_t x = 0; x < width; x++) {
            if (x) {
                std::cout << ",";
            }
            std::cout << get_num(matrix, width, height, x, y) << " ";
        }
        std::cout << "\n";
    }
}

/**
 * Применяет метод Гаусса-Жордана к матрице
 */
void solve_gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height) {
    std::size_t clear_y = 0;
    for (std::size_t x = 0, y = 0; x < width; x++) {
        std::size_t keep_y = find_1_index(matrix, width, height, x, clear_y);

        if (keep_y > height) {
            continue;
        }
        swap_lines(matrix, width, height, clear_y, keep_y);
        clear_y++;
        for (std::size_t y = clear_y; y < height; y++) {
            if (!get_num(matrix, width, height, x, y)) {
                continue;
            }
            sub_lines(matrix, width, height, y, clear_y - 1);
        }
    }
    for (std::size_t y = height - 1; y >= 1; y--) {
        std::size_t x_one = width + 1;
        for (std::size_t x = 0; x < width; x++) {
            if (get_num(matrix, width, height, x, y)) {
                x_one = x;
                break;
            }
        }
        if (x_one > width) {
            continue;
        }
        for (std::size_t sub_y = y - 1; ; sub_y--) {
            if (get_num(matrix, width, height, x_one, sub_y)) {
                sub_lines(matrix,width, height, sub_y, y);
            }
            if (!sub_y) {
                break;
            }
        }
    }
    if (width < 40 && height < 40) { // Если уместится печатаем
        print(matrix, width, height);
    }
}

/**
 * Находит линейную оболочку и вычисляет суммы для каждого из неизвестных
 */
std::pair<std::vector<int>, std::vector<std::vector<int>>> get_span_and_summs(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height) {
    std::vector<int> span(width, 0);
    std::vector<std::vector<int>> summs(width);
    for (std::size_t y = 0; y < height; y++) {
        std::vector<std::size_t> found;
        for (std::size_t x = 0; x < width; x++) {
            if (get_num(matrix, width, height, x, y)) {
                found.push_back(x);
            }
        }
        if (found.size() == 1) {
            span[found.back()] = 3;
        }
        int first = -1;
        for (auto x : found) {
            if (first == -1) {
                if (span[x] < 3) {
                    span[x] = 0;
                    first = x;
                }
            } else {
                if (span[x] < 3) {
                    summs[first].push_back(x);
                    span[x] = 1;
                }
            }
        }
    }
    return {span, summs};
}

/**
 * Находит индексы неизвестных, которые отвечают за линейную оболочку
 */
std::vector<std::size_t> get_span_offsets(std::vector<int>& span) {
    std::vector<std::size_t> offsets;
    for (std::size_t i = 0; i < span.size(); i++) {
        if (span[i] == 1) {
            offsets.push_back(i);
        }
    }
    return offsets;
}

/**
 * По данной оболочке и суммам, вычисляет зависимые элементы 
 */
void fill_by_span(std::vector<int>& result, std::vector<std::vector<int>>& summs, std::vector<int>& span) {
    for (std::size_t i = 0; i < result.size(); i++) {
        if (span[i] != 0) {
            continue ;
        }
        const auto& sum = summs[i];
        int summing = 0;
        for (const auto& x : sum) {
            summing = (summing + result[x]) % 2;
        }
        result[i] = summing;
    }
}

} // namespace

mpz_class gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors,
    std::vector<mpz_class>& lefts, mpz_class h, mpz_class m) {
    std::cout << "MATRIX: " << width << "x" << height << "\n";
    solve_gauss(matrix, width, height); // Применяем оболочку
    auto [span, sums] = get_span_and_summs(matrix, width, height); // Получаем оболочку и выражения для зависимых переменных
    auto offsets = get_span_offsets(span);
    std::vector<int> adds(offsets.size(), 0); // Значения для элементов оболочки
    while(true) {
        bool stop = true; // Остоновимся, если все элементы оболочки 1
        for (auto a : adds) {
            if (!a) {
                stop = false;
                break ;
            }
        }
        std::vector<int> result(span.size(), 0); // Результат - те элементы, которые нужно перемножать
        for (std::size_t i = 0; i < adds.size(); i++) {
            result[offsets[i]] = adds[i];
        }
        fill_by_span(result, sums, span); // Получаем его на основе линейной оболочки
        mpz_class left = 1;
        mpz_class right = 1;
        for (std::size_t j = 0; j < result.size(); j++) { // Перемножаем нужные элементы и индексы
            if (result[j]) {
                right = right * factors[j];
                left = left * (lefts[j] + h);
            }
        }
        if (sqrt(right) * sqrt(right) == right && right != 1) { // Если реально получили квадрат
            std::cout << "FOUND RIGHT:" << right << " SQRT:" << sqrt(right) << "\n";
            std::cout << "LEFT:" << left << "\n";
            right = sqrt(right); // Выполним ряд проверок, так иногда, например, (x - y(x)) = 0
            mpz_class result = gcd(left - right, m);
            std::cout << "GCD:" << result << "\n";
            if (result != 1 && result != m) {
                return result;
            }
            result = gcd(left + right, m);
            std::cout << "GCD:" << result << "\n";
            if (result != 1 && result != m) {
                return result;
            }
            std::cout << "BOTH BAD, CONTINUE:\n";
        }
        if (stop) {
            break ;
        }
        bool carry = true;
        for (std::size_t i = 0; i < adds.size(); i++) { // Получаем следующее состояние лин. оболочки
            int& a = adds[i];
            if (!a) {
                a = 1;
                break;
            }
            a = 0;
            carry = true;
        }
    }
  
    return m;
}
