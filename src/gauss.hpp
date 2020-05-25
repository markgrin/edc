#ifndef HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#define HPP_976D9043F68047C7BC9DC0E1FC46AFEE

#include <gmpxx.h>

#include <vector>

mpz_class gauss(std::vector<mpz_class>& matrix, std::size_t width, std::size_t height, std::vector<mpz_class>& factors);

#endif // HPP_976D9043F68047C7BC9DC0E1FC46AFEE