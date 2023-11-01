#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <set>

class DegreeFrac
{
private:
    int n, m;

    void Reduction();
public:
    DegreeFrac(int _n = 1, int _m = 1) : n(_n), m(_m) { Reduction(); }

    inline int GetN() const { return n; }
    inline int GetM() const { return m; }

    inline double ToDouble() const { return double(n) / m; }

    DegreeFrac& operator+=(const DegreeFrac& df);
    DegreeFrac operator+(const DegreeFrac& df) const;
    DegreeFrac operator-() const;
    DegreeFrac& operator*=(const DegreeFrac& df);
    DegreeFrac operator*(const DegreeFrac& df) const;

    bool operator<=(const DegreeFrac& df) const;
    inline bool operator>=(const DegreeFrac& df) const { return df <= *this; }
    inline bool operator==(const DegreeFrac& df) const { return *this <= df && df <= *this; }
    inline bool operator!=(const DegreeFrac& df) const { return !(df == *this); }
    inline bool operator<(const DegreeFrac& df) const { return *this <= df && df != *this; }
    inline bool operator>(const DegreeFrac& df) const { return *this >= df && df != *this; };

    friend std::ostream& operator<<(std::ostream& out, const DegreeFrac& df);
    friend std::istream& operator>>(std::istream& in, DegreeFrac& df);
};

std::ostream& operator<<(std::ostream& out, const DegreeFrac& df);
std::istream& operator>>(std::istream& in, DegreeFrac& df);

class PrimeDegreeFrac
{
private:
    int n;
    DegreeFrac df;
public:
    PrimeDegreeFrac(int _n, DegreeFrac _df) : n(_n), df(_df) {};

    inline int GetPrime() const { return n; }
    inline const DegreeFrac& GetDegreeFrac() const { return df; }

    inline double ToDouble() const { return std::pow(static_cast<double>(n), df.ToDouble()); }

    PrimeDegreeFrac& operator*=(const PrimeDegreeFrac& pdf);
    PrimeDegreeFrac operator*(const PrimeDegreeFrac& pdf) const;
    PrimeDegreeFrac& operator^=(const DegreeFrac& _df);
    PrimeDegreeFrac operator^(const DegreeFrac& _df) const;

    bool operator<=(const PrimeDegreeFrac& pdf) const;
    inline bool operator==(const PrimeDegreeFrac& pdf) const { return *this <= pdf && pdf <= *this; }
    inline bool operator!=(const PrimeDegreeFrac& pdf) const { return !(pdf == *this); }
    inline bool operator<(const PrimeDegreeFrac& pdf) const { return *this <= pdf && pdf != *this; };
};

class TropicoFrac
{
private:
    std::set<PrimeDegreeFrac> n;

    void Convert(int _n, const DegreeFrac& delta);
    bool ContainPrime(int _n) const;
    void AddDegree(int prime, DegreeFrac df);
    void PositiveDegree(TropicoFrac& tf1, TropicoFrac& tf2) const;

    DegreeFrac GetDegreePrime(int _n) const;
public:
    TropicoFrac(int _n = 0, int _m = 1);
    TropicoFrac(const PrimeDegreeFrac _prime) { n.insert(_prime); }
    TropicoFrac(const std::set<PrimeDegreeFrac>& _n) : n(_n) {}

    unsigned long long ValueIntPositive() const;

    double ToDouble() const;

    bool operator<=(const TropicoFrac& tf) const;
    inline bool operator>=(const TropicoFrac& tf) const { return tf <= *this; }
    inline bool operator==(const TropicoFrac& tf) const { return tf >= *this && tf <= *this; }
    inline bool operator!=(const TropicoFrac& tf) const { return !(tf == *this); }
    inline bool operator<(const TropicoFrac& tf) const { return *this <= tf && tf != *this; };
    inline bool operator>(const TropicoFrac& tf) const { return *this >= tf && tf != *this; };

    TropicoFrac& operator+=(const TropicoFrac& tf);
    TropicoFrac operator+(const TropicoFrac& tf) const;
    TropicoFrac& operator*=(const TropicoFrac& tf);
    TropicoFrac operator*(const TropicoFrac& tf) const;
    TropicoFrac& operator^=(const DegreeFrac& df);
    TropicoFrac operator^(const DegreeFrac& df) const;

    friend std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tf);
    friend std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf);
    friend std::istream& operator>>(std::istream& in, TropicoFrac& tf);
};

std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tm);
std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf);
std::istream& operator>>(std::istream& in, TropicoFrac& tf);

class TropicoMatrix
{
private:
    std::vector<std::vector<TropicoFrac>> matr;
public:
    TropicoMatrix(size_t height = 0, size_t width = 0, TropicoFrac elem = TropicoFrac());

    inline size_t GetHeight() const { return matr.size(); }
    size_t GetWidth() const;

    static TropicoMatrix GetI(size_t size);

    void Resize(size_t height, size_t width);

    TropicoMatrix SubMatrix(size_t x, size_t y, size_t height, size_t width) const;
    TropicoMatrix MultiConjTrans() const;
    void Standardization();

    std::vector<std::vector<double>> ToDouble() const;

    const std::vector<TropicoFrac>& operator[](size_t i) const { return matr[i]; }
    std::vector<TropicoFrac>& operator[](size_t i) { return matr[i]; }

    TropicoMatrix& operator+=(const TropicoMatrix& tm);
    TropicoMatrix operator+(const TropicoMatrix& tm) const;
    TropicoMatrix& operator*=(const TropicoMatrix& tm);
    TropicoMatrix operator*(const TropicoMatrix& tm) const;
    TropicoMatrix& operator*=(const TropicoFrac& tf);
    TropicoMatrix operator*(const TropicoFrac& tf) const;

    TropicoFrac Trace() const;
    TropicoFrac SpectralRadius() const;
    TropicoMatrix Kleene() const;

    friend bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2);

    friend std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm);
    friend std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm);
    friend std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm);
    friend std::istream& operator>>(std::istream& in, TropicoMatrix& tm);
};

bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2);

inline TropicoMatrix operator*(const TropicoFrac& tf, const TropicoMatrix& tm) { return tm * tf; }

std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm);
std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm);
std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm);
std::istream& operator>>(std::istream& in, TropicoMatrix& tm);

class TropicoSolve
{
private:
    TropicoMatrix tm;
    TropicoFrac spectral;
    TropicoMatrix kleene;
    TropicoMatrix worst;
    TropicoMatrix best;
    TropicoFrac delta;
    TropicoMatrix P;
    TropicoMatrix Plk;

    std::string name;
public:
    TropicoSolve(TropicoMatrix _tm = TropicoMatrix());

    TropicoMatrix WorstSolve();
    TropicoMatrix MakeP() const;
    TropicoMatrix MakePlk() const;
    TropicoMatrix BestSolve();

    inline void GiveName(std::string _name) { name = _name; }

    inline const TropicoMatrix& GetMatrix() const { return tm; }
    inline const TropicoFrac& GetSpectral() const { return spectral; }
    inline const TropicoMatrix& GetKleene() const { return kleene; }
    inline const TropicoMatrix& GetWorst() const { return worst; }
    inline const TropicoMatrix& GetBest() const { return best; }
    inline const TropicoFrac& GetDelta() const { return delta; }
    inline const TropicoMatrix& GetP() const { return P; }
    inline const TropicoMatrix& GetPlk() const { return Plk; }

    inline const std::string& GetName() const { return name; }

    void Solve();

    friend std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts);
};

std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts);

class TropicoMultiSolve
{
private:
    TropicoMatrix C;
    TropicoMatrix BWorst, BBest;
    TropicoSolve tsC, tsBW, tsBB;
    std::vector<TropicoMatrix> A;
public:
    TropicoMultiSolve() = default;
    TropicoMultiSolve(TropicoMatrix& _C, std::vector<TropicoMatrix>& _A) : C(_C), A(_A) {}
    
    void Solve();

    friend std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms);
    friend std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms);
};

std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms);
std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms);