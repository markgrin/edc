#include <gmpxx.h>
#include <iostream>

#include "ressol.hpp"

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
            std::cout << "Input numbers:";
            std::cin >> number;
        } else if (argv == 2) { // Берем число из аргументов
            number = argc[2];
        } else { // Аргументов больше чем нужно, завершаемся
            return help();
        }

        auto input_number = mpz_class(number.data());
        if (input_number < 1) {
            std::cout << "Bad input\n";
            return 1;
        }

        // В будущем здесь можно вызывать функцию факторизации, но пока алгоритм Тонелли-Шенкса

    } catch (...) {
        std::cout << "Error!\n";
        return 1;
    }
    return 0;
    
}