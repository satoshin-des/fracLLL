#include "frac_lll.hpp"

#include <NTL/ZZ.h>

void FracLLL::updateGSOSwap(const long k)
{
    long i, j;
    this->m_C[NUM] = this->B[NUM];
    this->m_C[DEN] = this->B[DEN];
    this->m_nu[NUM] = this->mu[NUM];
    this->m_nu[DEN] = this->mu[DEN];

    this->m_C[NUM][k - 1] = this->B[NUM][k] * this->B[DEN][k - 1] * this->mu[DEN][k][k - 1] * this->mu[DEN][k][k - 1] + this->mu[NUM][k][k - 1] * this->mu[NUM][k][k - 1] * this->B[DEN][k] * this->B[NUM][k - 1];
    this->m_C[DEN][k - 1] = this->B[DEN][k] * this->B[DEN][k - 1] * this->mu[DEN][k][k - 1] * this->mu[DEN][k][k - 1];
    m_d = NTL::GCD(this->m_C[NUM][k - 1], this->m_C[DEN][k - 1]);
    this->m_C[NUM][k - 1] /= m_d;
    this->m_C[DEN][k - 1] /= m_d;
    this->m_C[NUM][k] = this->B[NUM][k - 1] * this->B[NUM][k] * this->m_C[DEN][k - 1];
    this->m_C[DEN][k] = this->B[DEN][k - 1] * this->B[DEN][k] * this->m_C[NUM][k - 1];
    m_d = NTL::GCD(this->m_C[NUM][k], this->m_C[DEN][k]);
    this->m_C[NUM][k] /= m_d;
    this->m_C[DEN][k] /= m_d;

    this->m_nu[NUM][k][k - 1] = this->mu[NUM][k][k - 1] * this->B[NUM][k - 1] * this->m_C[DEN][k - 1];
    this->m_nu[DEN][k][k - 1] = this->mu[DEN][k][k - 1] * this->B[DEN][k - 1] * this->m_C[NUM][k - 1];
    m_d = NTL::GCD(this->m_nu[NUM][k][k - 1], this->m_nu[DEN][k][k - 1]);
    this->m_nu[NUM][k][k - 1] /= m_d;
    this->m_nu[DEN][k][k - 1] /= m_d;
    for (j = 0; j < k - 1; ++j)
    {
        this->m_nu[NUM][k - 1][j] = this->mu[NUM][k][j];
        this->m_nu[DEN][k - 1][j] = this->mu[DEN][k][j];
        m_d = NTL::GCD(this->m_nu[NUM][k - 1][j], this->m_nu[DEN][k - 1][j]);
        this->m_nu[NUM][k - 1][j] /= m_d;
        this->m_nu[DEN][k - 1][j] /= m_d;
        this->m_nu[NUM][k][j] = this->mu[NUM][k - 1][j];
        this->m_nu[DEN][k][j] = this->mu[DEN][k - 1][j];
        m_d = NTL::GCD(this->m_nu[NUM][k][j], this->m_nu[DEN][k][j]);
        this->m_nu[NUM][k][j] /= m_d;
        this->m_nu[DEN][k][j] /= m_d;
    }
    for (i = k + 1; i < this->n; ++i)
    {
        this->m_nu[NUM][i][k] = this->mu[NUM][i][k - 1] * this->mu[DEN][k][k - 1] * this->mu[DEN][i][k] - this->mu[NUM][k][k - 1] * this->mu[NUM][i][k] * this->mu[DEN][i][k - 1];
        this->m_nu[DEN][i][k] = this->mu[DEN][i][k - 1] * this->mu[DEN][k][k - 1] * this->mu[DEN][i][k];
        m_d = NTL::GCD(this->m_nu[NUM][i][k], this->m_nu[DEN][i][k]);
        this->m_nu[NUM][i][k] /= m_d;
        this->m_nu[DEN][i][k] /= m_d;
        this->m_nu[NUM][i][k - 1] = this->mu[NUM][i][k] * this->m_nu[DEN][k][k - 1] * this->m_nu[DEN][i][k] + this->mu[DEN][i][k] * this->m_nu[NUM][k][k - 1] * this->m_nu[NUM][i][k];
        this->m_nu[DEN][i][k - 1] = this->mu[DEN][i][k] * this->m_nu[DEN][k][k - 1] * this->m_nu[DEN][i][k];
        m_d = NTL::GCD(this->m_nu[NUM][i][k - 1], this->m_nu[DEN][i][k - 1]);
        this->m_nu[NUM][i][k - 1] /= m_d;
        this->m_nu[DEN][i][k - 1] /= m_d;
    }

    this->B[NUM] = this->m_C[NUM];
    this->B[DEN] = this->m_C[DEN];
    this->mu[NUM] = this->m_nu[NUM];
    this->mu[DEN] = this->m_nu[DEN];
}

void FracLLL::fracLLL(const long a, const long b)
{
    NTL::vec_ZZ t;

    this->fracGSO();

    for (long k = 1, j; k < this->n;)
    {
        for (j = k - 1; j > -1; --j)
        {
            fracSize(k, j);
        }

        if (this->B[NUM][k] * this->B[DEN][k - 1] * b * this->mu[DEN][k][k - 1] * this->mu[DEN][k][k - 1] >= (a * this->mu[DEN][k][k - 1] * this->mu[DEN][k][k - 1] - b * this->mu[NUM][k][k - 1] * this->mu[NUM][k][k - 1]) * this->B[NUM][k - 1] * this->B[DEN][k])
        {
            ++k;
        }
        else
        {
            this->basis[k].swap(this->basis[k - 1]);
            this->updateGSOSwap(k);
            // this->fracGSO();

            if (k > 2)
            {
                --k;
            }
            else
            {
                k = 1;
            }
        }
    }
}
