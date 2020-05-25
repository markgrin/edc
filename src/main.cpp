#include "sieve.hpp"

#include <gmpxx.h>
#include <iostream>


int help () {
    std::cout << "factorization programm.\n"
                 "USAGE: factorization [number]\n"
                 "number - number to factorize\n"
                 "\n"
                 "If called with no arguments they will be promted at runtime\n";
    return 1;
}

int main (int argv, char ** argc) {
    try {
        std::string number;
        if (argv == 1) { // Аргументов нет, запрашиваем число
            std::cout << "Input number:";
            std::cin >> number;
        } else if (argv == 2) { // Берем число из аргументов
            number = argc[1];
        } else { // Аргументов больше чем нужно, завершаемся
            return help();
        }

        auto input_number = mpz_class(number.data());
        if (input_number < 1) {
            std::cout << "Bad input\n";
            return 1;
        } else if (input_number <= 3) {
            std::cout << input_number << "\n";
            return 0;
        }

        auto result = quadratic_sieve_algorithm(input_number);
        std::cout << "RESULT:" << result << "\n";
        if (input_number % result != 0) {
            std::cout << "INCORRECT\n";
            return 1;
        } else {
            std::cout << "CORRECT\n";
            return 0;
        }

    } catch (...) {
        std::cout << "Error!\n";
        return 1;
    }
    return 0;
    
}