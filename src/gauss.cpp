#include "gauss.hpp"

#include <iostream>
#include <utility>

namespace {


mpz_class slow_solver(std::vector<mpz_class>& factors, std::vector<mpz_class>& lefts, mpz_class h, mpz_class m) {
    std::vector<int> adds(factors.size(), 0);
    mpz_class result = 1;
    while(true) {
        bool stop = true;
        for (auto a : adds) {
            if (!a) {
                stop = false;
                break ;
            }
        }
        mpz_class right = 1;
        mpz_class left = 1;
        for (std::size_t j = 0; j < factors.size(); j++) {
            if (adds[j]) {
                right = right * factors[j];
                left = left * (lefts[j] + h);
            }
        }
        if (sqrt(right) * sqrt(right) == right && right != 1) {
            ////std::cout << "FOUND RIGHT:" << right << " SQRT:" << sqrt(right) << "\n";
            ////std::cout << "LEFT:" << left << "\n";
            right = sqrt(right);
            mpz_class result = gcd(left - right, m);
            ////std::cout << "GCD:" << result << "\n";
            if (result != 1 && result != m) {
                return result;
            }
            result = gcd(left + right, m);
            ////std::cout << "GCD:" << result << "\n";
            if (result != 1 && result != m) {
                return result;
            }
            ////std::cout << "BOTH BAD, CONTINUE:\n";
            //for (std::size_t x = 0; x < adds.size(); x++) {
            //////    std::cout << adds[x] << " ";
            //}
        }
        bool carry = true;
        if (stop) {
            break ;
        }
        for (std::size_t i = 0; i < adds.size(); i++) {
            int& a = adds[i];
            if (!a) {
                a = true;
                break;
            }
            a = false;
            carry = true;
        }
    }
    return 1;
}

mpz_class& get_num(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t x, std::size_t y) {
    return matrix[x * height + y];
    //return matrix[y * width + x];
}

void sub_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j) {
    for (std::size_t x = 0; x < width; x++) {
        mpz_class& num = get_num(matrix, width, height, x, i);
        num = (num + get_num(matrix, width, height, x, j)) % 2;
    }
}

void swap_lines(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t j){//, std::vector<mpz_class>& factors) {
    for (std::size_t x = 0; x < width; x++) {
        std::swap(get_num(matrix, width, height, x, i), get_num(matrix, width, height, x, j));
    }
   // std::swap(factors[i], factors[j]);
}

std::size_t find_1_index(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::size_t i, std::size_t clear_y) {
    for (std::size_t y = clear_y; y < height; y++) {
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
////    std::cout << "[\n";
    for (std::size_t y = 0; y < height; y++) {
        if (y <= 9) {
            std::cout << " ";
        }
        std::cout << y << "]  ";
////        std::cout << "[\n";
        for (std::size_t x = 0; x < width; x++) {
            if (x) {
                std::cout << ",";
            }
            std::cout << get_num(matrix, width, height, x, y) << " ";
        }
////        std::cout << "]";
        if (y != height - 1) {
////            std::cout << ",";
        }
        std::cout << "\n";
    }
////    std::cout << "]\n";
}

void solve_gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height) {
    //std::cout << "\nBEFORE:" << "\n";
    //print(matrix, width, height);
    std::size_t clear_y = 0;
    for (std::size_t x = 0, y = 0; x < width; x++) {
        std::size_t keep_y = find_1_index(matrix, width, height, x, clear_y);

        //std::cout << x << " " << clear_y << "\n";
        if (keep_y > height) {
            continue;
        }
        swap_lines(matrix, width, height, clear_y, keep_y);// {});//factors);
        //std::cout << "SWAP:" << clear_y << " " << keep_y << "\n";
        clear_y++;
        for (std::size_t y = clear_y; y < height; y++) {
            if (!get_num(matrix, width, height, x, y)) { // clear_y == x
                continue;
            }
            //std::cout << "SUB:" << clear_y - 1 << " " << y << "\n";
            sub_lines(matrix, width, height, y, clear_y - 1);
        }
    //std::cout << "\nSTEP:" << "\n";
    //print(matrix, width, height);
    }
    //std::cout << "\nBACK STEP:" << "\n";
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
    if (width < 40 && height < 40) {
        print(matrix, width, height);
    }
}

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

std::vector<std::size_t> get_span_offsets(std::vector<int>& span) {
    std::vector<std::size_t> offsets;
    for (std::size_t i = 0; i < span.size(); i++) {
        if (span[i] == 1) {
            offsets.push_back(i);
        }
    }
    return offsets;
}

void fill_by_span(std::vector<int>& result, std::vector<std::vector<int>>& summs, std::vector<int>& span) {
    for (std::size_t i = 0; i < result.size(); i++) {
        if (span[i] != 0) {
            continue ;
        }
        const auto& sum = summs[i];
        //std::cout << "SUMMSSIZE:" << sum.size() << "\n";
        int summing = 0;
        for (const auto& x : sum) {
        //std::cout << "SUMMADD:" << x << "=" <<  result[x] << "\n";
            summing = (summing + result[x]) % 2;
        }
        result[i] = summing;
    }
}

} // namespace

mpz_class gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors,
    std::vector<mpz_class>& lefts, mpz_class h, mpz_class m) {
    std::cout << "MATRIX: " << width << "x" << height << "\n";
    //std::cout << "GAUSS FACTORS:" << "\n";
    for (auto factor : factors) {
        //std:: cout << factor << ", ";
    }
    solve_gauss(matrix, width, height);
    auto [span, sums] = get_span_and_summs(matrix, width, height);
    auto offsets = get_span_offsets(span);
    std::vector<int> adds(offsets.size(), 0);
    while(true) {
        bool stop = true;
        for (auto a : adds) {
            if (!a) {
                stop = false;
                break ;
            }
        }
        std::vector<int> result(span.size(), 0);
        for (std::size_t i = 0; i < adds.size(); i++) {
            result[offsets[i]] = adds[i];
        }
        fill_by_span(result, sums, span);
        mpz_class left = 1;
        mpz_class right = 1;
        for (std::size_t j = 0; j < result.size(); j++) {
            if (result[j]) {
                right = right * factors[j];
                left = left * (lefts[j] + h);
            }
        }
        if (sqrt(right) * sqrt(right) == right && right != 1) {
            std::cout << "FOUND RIGHT:" << right << " SQRT:" << sqrt(right) << "\n";
            std::cout << "LEFT:" << left << "\n";
            right = sqrt(right);
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
            //for (std::size_t x = 0; x < adds.size(); x++) {
            //////    std::cout << adds[x] << " ";
            //}
        }
        bool carry = true;
        if (stop) {
            break ;
        }
        for (std::size_t i = 0; i < adds.size(); i++) {
            int& a = adds[i];
            if (!a) {
                a = 1;
                break;
            }
            a = 0;
            carry = true;
        }
    }
  
    //////std::cout << std::endl << "\n";
    //std::vector<int> result(span.size(), 0);
 
    //////std::cout << "RESULT:\n" << std::endl;
    //mpz_class right = 1;
    //for (std::size_t i = 0; i < result.size(); i++) {
    //////    std::cout << result[i] << " ,";
    //    if (result[i]) {
    //        right *= factors[i];
    //    }
    //}
    //////std::cout << "\nAFTER:" << "\n";
    //print(matrix, width, height);
    return m;
}
