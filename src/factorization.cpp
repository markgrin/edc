#include "factorization.hpp"

std::vector<mpz_class> factorization(const mpz_class& f, std::vector<mpz_class>& v)
{
    std::vector<mpz_class> result = {};
    mpz_class current;
    mpz_class nf = f;

    for (mpz_class& n : v)
    {
        current = 0;
        while (nf % n == 0 && abs(nf) >= n)
        {
            current++;
            nf = nf / n;
        }
        result.push_back(current);
    }
    if(!(nf==1||nf==-1))
    {
        result.erase(result.begin(), result.end());
        result.push_back(0);
    }
    else
    {
        result[0] = result[0] * nf;
    }
    return result;
} 
