#include "gauss.hpp"

#include <iostream>
#include <utility>

namespace {

mpz_class& get_num(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t x, std::size_t y) {
    //return matrix[x * height + y];
    return matrix[y * width + x];
}

void sub_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j) {
    for (std::size_t x = 0; x < width; x++) {
        mpz_class& num = get_num(matrix, width, height, x, i);
        num = (num + get_num(matrix, width, height, x, j)) % 2;
    }
}

void swap_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j, std::vector<mpz_class>& factors) {
    for (std::size_t x = 0; x < width; x++) {
        std::swap(get_num(matrix, width, height, x, i), get_num(matrix, width, height, x, j));
    }
    std::swap(factors[i], factors[j]);
}

std::size_t find_1_index(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i) {
    for (std::size_t y = i; y < height; y++) {
        if (get_num(matrix, width, height, i, y)) {
            return y;
        }
    }
    return height + 1;
}

bool is_ones (std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t y) {
    for (std::size_t x = 0; x < width; x++) {
        if (get_num(matrix, width, height, x, y)) {
            return true;
        }
    }
    return false;
}

void print(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height) {
    for (std::size_t y = 0; y < height; y++) {
        for (std::size_t x = 0; x < width; x++) {
            std::cout << get_num(matrix, width, height, x, y) << " ";
        }
        std::cout << "\n";
    }
}

}

mpz_class gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors) {
    std::cout << "BEFORE:" << "\n";
    print(matrix, width, height);
    for (std::size_t x = 0; x < width; x++) {
        std::size_t keep_y = find_1_index(matrix, width, height, x);
        if (keep_y > height) {
            continue;
        }
        swap_lines(matrix, width, height, x, keep_y, factors);
        for (std::size_t y = x + 1; y < height; y++) {
            if (!get_num(matrix, width, height, x, y)) {
                continue;
            }
            sub_lines(matrix, width, height, y, x);
        }
    }
    mpz_class result = 1;
    for (std::size_t y = 0; y < height; y++) {
        if (is_ones(matrix, width, height, y)) {
            result = result * factors[y];
        }
    }
    std::cout << "AFTER:" << "\n";
    print(matrix, width, height);
    return result;
}
