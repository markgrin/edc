#ifndef HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#define HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#include <gmpxx.h>

#include <vector>

mpz_class slow_solver(std::vector<mpz_class>& factors, std::vector<mpz_class>& lefts, mpz_class h, mpz_class m);

std::pair<mpz_class, std::vector<int>> gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors,
          std::vector<mpz_class>& lefts, mpz_class h, mpz_class m);

#endif // HPP_976D9043F68047C7BC9DC0E1FC46AFEE