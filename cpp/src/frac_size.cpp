#include "frac_lll.hpp"

#include <NTL/ZZ.h>
#include <NTL/RR.h>

void FracLLL::fracSize(const long i, const long j)
{
    if ((NTL::abs(this->mu[NUM][i][j]) << 1) > NTL::abs(this->mu[DEN][i][j]))
    {
        const NTL::ZZ q = NTL::RoundToZZ(NTL::to_RR(this->mu[NUM][i][j]) / NTL::to_RR(this->mu[DEN][i][j]));
        this->basis[i] -= q * this->basis[j];
        
        for (long h = 0; h <= j; ++h)
        {
            this->mu[NUM][i][h] = this->mu[NUM][i][h] * this->mu[DEN][j][h] - q * this->mu[DEN][i][h] * this->mu[NUM][j][h];
            this->mu[DEN][i][h] *= this->mu[DEN][j][h];
            this->m_d = NTL::GCD(this->mu[NUM][i][h], this->mu[DEN][i][h]);
            this->mu[NUM][i][h] /= this->m_d;
            this->mu[DEN][i][h] /= this->m_d;
        }
    }
}
