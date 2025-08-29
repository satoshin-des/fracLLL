#include "frac_lll.hpp"

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

void FracLLL::updateGSODeepInsertion(const long i, const long k)
{
    this->m_C[NUM] = this->B[NUM];
    this->m_C[DEN] = this->B[DEN];
    this->m_nu[NUM] = this->mu[NUM];
    this->m_nu[DEN] = this->mu[DEN];

    for (long l = k, j; l > i; --l)
    {
        this->m_C[NUM][l - 1] = this->B[NUM][l] * this->B[DEN][l - 1] * this->mu[DEN][l][l - 1] * this->mu[DEN][l][l - 1] + this->mu[NUM][l][l - 1] * this->mu[NUM][l][l - 1] * this->B[DEN][l] * this->B[NUM][l - 1];
        this->m_C[DEN][l - 1] = this->B[DEN][l] * this->B[DEN][l - 1] * this->mu[DEN][l][l - 1] * this->mu[DEN][l][l - 1];
        m_d = NTL::GCD(this->m_C[NUM][l - 1], this->m_C[DEN][l - 1]);
        this->m_C[NUM][l - 1] /= m_d;
        this->m_C[DEN][l - 1] /= m_d;
        this->m_C[NUM][l] = this->B[NUM][l - 1] * this->B[NUM][l] * this->m_C[DEN][l - 1];
        this->m_C[DEN][l] = this->B[DEN][l - 1] * this->B[DEN][l] * this->m_C[NUM][l - 1];
        m_d = NTL::GCD(this->m_C[NUM][l], this->m_C[DEN][l]);
        this->m_C[NUM][l] /= m_d;
        this->m_C[DEN][l] /= m_d;

        this->m_nu[NUM][l][l - 1] = this->mu[NUM][l][l - 1] * this->B[NUM][l - 1] * this->m_C[DEN][l - 1];
        this->m_nu[DEN][l][l - 1] = this->mu[DEN][l][l - 1] * this->B[DEN][l - 1] * this->m_C[NUM][l - 1];
        m_d = NTL::GCD(this->m_nu[NUM][l][l - 1], this->m_nu[DEN][l][l - 1]);
        this->m_nu[NUM][l][l - 1] /= m_d;
        this->m_nu[DEN][l][l - 1] /= m_d;
        for (j = 0; j < l - 1; ++j)
        {
            this->m_nu[NUM][l - 1][j] = this->mu[NUM][l][j];
            this->m_nu[DEN][l - 1][j] = this->mu[DEN][l][j];
            m_d = NTL::GCD(this->m_nu[NUM][l - 1][j], this->m_nu[DEN][l - 1][j]);
            this->m_nu[NUM][l - 1][j] /= m_d;
            this->m_nu[DEN][l - 1][j] /= m_d;
            this->m_nu[NUM][l][j] = this->mu[NUM][l - 1][j];
            this->m_nu[DEN][l][j] = this->mu[DEN][l - 1][j];
            m_d = NTL::GCD(this->m_nu[NUM][l][j], this->m_nu[DEN][l][j]);
            this->m_nu[NUM][l][j] /= m_d;
            this->m_nu[DEN][l][j] /= m_d;
        }
        for (j = l + 1; j < this->n; ++j)
        {
            this->m_nu[NUM][j][l] = this->mu[NUM][j][l - 1] * this->mu[DEN][l][l - 1] * this->mu[DEN][j][l] - this->mu[NUM][l][l - 1] * this->mu[NUM][j][l] * this->mu[DEN][j][l - 1];
            this->m_nu[DEN][j][l] = this->mu[DEN][j][l - 1] * this->mu[DEN][l][l - 1] * this->mu[DEN][j][l];
            m_d = NTL::GCD(this->m_nu[NUM][j][l], this->m_nu[DEN][j][l]);
            this->m_nu[NUM][j][l] /= m_d;
            this->m_nu[DEN][j][l] /= m_d;
            this->m_nu[NUM][j][l - 1] = this->mu[NUM][j][l] * this->m_nu[DEN][l][l - 1] * this->m_nu[DEN][j][l] + this->mu[DEN][j][l] * this->m_nu[NUM][l][l - 1] * this->m_nu[NUM][j][l];
            this->m_nu[DEN][j][l - 1] = this->mu[DEN][j][l] * this->m_nu[DEN][l][l - 1] * this->m_nu[DEN][j][l];
            m_d = NTL::GCD(this->m_nu[NUM][j][l - 1], this->m_nu[DEN][j][l - 1]);
            this->m_nu[NUM][j][l - 1] /= m_d;
            this->m_nu[DEN][j][l - 1] /= m_d;
        }

        this->B[NUM] = this->m_C[NUM];
        this->B[DEN] = this->m_C[DEN];
        this->mu[NUM] = this->m_nu[NUM];
        this->mu[DEN] = this->m_nu[DEN];
    }
}

void FracLLL::fracDeepLLL(const long a, const long b)
{
    NTL::ZZ prod, sum, g;
    NTL::vec_ZZ t;

    this->fracGSO();

    for (long k = 1, i, j; k < this->n;)
    {
        printf("k = %ld\n", k);
        for (j = k - 1; j > -1; --j)
        {
            fracSize(k, j);
        }

        prod = 1;
        sum = 0;
        NTL::InnerProduct(g, this->basis[k], this->basis[k]);
        for (i = 0; i < k;)
        {
            if ((b * g * this->B[DEN][i] - a * this->B[NUM][i]) * prod >= b * this->B[DEN][i] * sum)
            {
                prod *= this->mu[DEN][k][i] * this->mu[DEN][k][i] * this->B[DEN][i];
                ++i;
                sum = 0;
                for (j = 0; j < i; ++j)
                {
                    sum += (this->mu[NUM][k][j] * this->B[NUM][j] * prod) / (this->mu[DEN][k][j] * this->mu[DEN][k][j] * this->B[DEN][j]);
                }
            }
            else
            {
                t = this->basis[k];
                for (j = k; j > i; --j)
                {
                    this->basis[j] = this->basis[j - 1];
                }
                this->basis[i] = t;
                this->updateGSODeepInsertion(i, k);

                if (i > 1)
                {
                    k = i - 1;
                }
                else
                {
                    k = 0;
                }
            }
        }
        ++k;
    }
}
