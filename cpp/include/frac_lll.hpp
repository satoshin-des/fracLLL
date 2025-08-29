#ifndef FRAC_LLL_HPP
#define FRAC_LLL_HPP

#include <iostream>

#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#define NUM 0
#define DEN 1

class FracLLL
{
private:
    NTL::mat_ZZ m_G;
    NTL::mat_ZZ m_nu[2];
    NTL::vec_ZZ m_C[2];
    NTL::ZZ m_d;

public:
    long n;
    NTL::mat_ZZ basis;
    NTL::mat_ZZ mu[2];
    NTL::vec_ZZ B[2];

    void svpChallenge(const long dim, const long seed);

    /**
     * @brief fracGSO
     *
     */
    void fracGSO();

    /**
     * @brief fracSize
     *
     * @param i index
     * @param j index
     */
    void fracSize(const long i, const long j);

    /**
     * @brief update GSO-information by swapping b_{k - 1} and b_k
     *
     * @param k index
     */
    void updateGSOSwap(const long k);

    /**
     * @brief 
     * 
     * @param i 
     * @param k 
     */
    void updateGSODeepInsertion(const long i, const long k);

    /**
     * @brief fracLLL
     *
     * @param a numerator of reduction parameter
     * @param b denominator of reduction parameter
     */
    void fracLLL(const long a, const long b);

    /**
     * @brief fracDeepLLL
     * 
     * @param a numerator of reduction parameter
     * @param b denominator of reduction parameter
     */
    void fracDeepLLL(const long a, const long b);
};

#endif // !FRAC_LLL_H
