#ifndef _ACCUMULATOR_HPP
#define _ACCUMULATOR_HPP
// Original source location: http://roth.cs.kuleuven.be/w-ess/index.php/File:Variance.tar.bz2
// supporting description at http://roth.cs.kuleuven.be/w-ess/index.php/Accurate_variance_and_mean_calculations_in_C%2B%2B11#Source_code

//////////////////////////////////////////////////////////////////
// ESS Seminar #24:
//   Accurate variance and mean calculations in C++11
//   Dirk Nuyens
//////////////////////////////////////////////////////////////////

template <typename T, typename T2=T>
struct accumulator
{
    T2 sum;
    T S;
    T M;
    size_t N;

    accumulator() : sum(0), S(0), M(0), N(0) { }

    T2 operator()(const T& x)
    {
        ++N;
        sum += x;
        T Mprev = M;
        M += (x - Mprev) / N;
        S += (x - Mprev) * (x - M);
        return sum;
    }

    T mean() const
    {
        return sum / N;
    }

    T variance() const
    {
        return S / (N - 1);
    }

    friend std::ostream& operator<<(std::ostream& out,
            const accumulator& a)
    {
        out << " N         = " << a.N << std::endl
            << " sum       = " << a.sum << std::endl
            << " mean      = " << a.mean() << std::endl
            << " variance  = " << a.variance() << std::endl;
        return out;
    }

};

#endif

