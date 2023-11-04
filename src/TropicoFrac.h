/// @mainpage Tropic Math for Solve Decision Making
/// @author Oleynik Michael
/// @date 04.10.2023

/// @file
/// @brief Classes Tropic frac, Tropic matrix and solver.
/// @details File contains the definition of class of a degree fraction, class of a prime number in degree of fraction, 
/// class of a tropic number, class of a tropic matrix, class of solver of a single-criteria task and 
/// class of a solver multi-criteria task.
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <set>
#include <algorithm>
#include "BigInt.h"

/// @brief Class of degree fraction.
class DegreeFrac
{
private:
    int n, m;

    void Reduction();
public:
    /// @brief Constructor.
    /// @param _n Numerator.
    /// @param _m Denominator.
    DegreeFrac(int _n = 1, int _m = 1) : n(_n), m(_m) { Reduction(); }

    /// @brief Getter.
    /// @return Nominator.
    inline int GetN() const { return n; }
    /// @brief Getter.
    /// @return Denominator.
    inline int GetM() const { return m; }

    /// @brief Converts to double.
    /// @return Value of fraction.
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

/// @brief Class of prime number in degree of fraction.
class PrimeDegreeFrac
{
private:
    int n;
    DegreeFrac df;
public:
    /// @brief Constructor.
    /// @param _n Prime number.
    /// @param _df Degree fraction.
    PrimeDegreeFrac(int _n, DegreeFrac _df) : n(_n), df(_df) {};

    /// @brief Getter.
    /// @return Prime number.
    inline int GetPrime() const { return n; }
    /// @brief Getter.
    /// @return Degree fraction.
    inline const DegreeFrac& GetDegreeFrac() const { return df; }

    /// @brief Converts to double.
    /// @return Value of the number.
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

/// @brief Class of tropic number
class TropicoFrac
{
private:
    std::set<PrimeDegreeFrac> n;

    /// @brief Converts integer number to tropic number.
    /// @param _n Integer number.
    /// @param[in] delta Degree (optional).
    void Convert(int _n, const DegreeFrac& delta = 1);
    /// @brief Checks contain of a prime number.
    /// @param _n Prime number.
    /// @return 
    bool ContainPrime(int _n) const;
    /// @brief Multiplies a number by a prime number in degree.
    /// @param prime Prime number.
    /// @param df Degree.
    void AddDegree(int prime, DegreeFrac df);
    void PositiveDegree(TropicoFrac& tf1, TropicoFrac& tf2) const;

    /// @brief Getter.
    /// @param _n Prime number.
    /// @return Degree.
    DegreeFrac GetDegreePrime(int _n) const;
public:
    /// @brief Constructor.
    /// @param _n Numerator.
    /// @param _m Denominator.
    TropicoFrac(int _n = 0, int _m = 1);
    TropicoFrac(const PrimeDegreeFrac _prime) { n.insert(_prime); }
    TropicoFrac(const std::set<PrimeDegreeFrac>& _n) : n(_n) {}

    /// @brief Converts to integer if all degree is natural.
    /// @return Value of the number.
    bigint ValueIntPositive() const;

    /// @brief Converts to double.
    /// @return Value of the number.
    double ToDouble() const;

    bool operator<=(const TropicoFrac& tf) const;
    inline bool operator>=(const TropicoFrac& tf) const { return tf <= *this; }
    inline bool operator==(const TropicoFrac& tf) const { return tf >= *this && tf <= *this; }
    inline bool operator!=(const TropicoFrac& tf) const { return !(tf == *this); }
    inline bool operator<(const TropicoFrac& tf) const { return *this <= tf && tf != *this; };
    inline bool operator>(const TropicoFrac& tf) const { return *this >= tf && tf != *this; };

    /// @brief Maximum of two numbers.
    /// @param tf 
    /// @return 
    TropicoFrac& operator+=(const TropicoFrac& tf);
    /// @brief Maximum of two numbers.
    /// @param tf 
    /// @return 
    TropicoFrac operator+(const TropicoFrac& tf) const;
    TropicoFrac& operator*=(const TropicoFrac& tf);
    TropicoFrac operator*(const TropicoFrac& tf) const;
    TropicoFrac& operator^=(const DegreeFrac& df);
    TropicoFrac operator^(const DegreeFrac& df) const;

    friend std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tf);
    friend std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf);
    friend std::istream& operator>>(std::istream& in, TropicoFrac& tf);
};

/// @brief Outputs a number in TeX format.
/// @param out Output stream.
/// @param tm Number.
/// @return Output stream.
std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tm);
std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf);
std::istream& operator>>(std::istream& in, TropicoFrac& tf);

/// @brief Class of tropic matrix
class TropicoMatrix
{
private:
    std::vector<std::vector<TropicoFrac>> matr;
public:
    TropicoMatrix(size_t height = 0, size_t width = 0, TropicoFrac elem = TropicoFrac());

    /// @brief Getter.
    /// @return Height of the matrix.
    inline size_t GetHeight() const { return matr.size(); }
    /// @brief Getter.
    /// @return Width of the matrix.
    size_t GetWidth() const;

    /// @brief Gets the unit matrix.
    /// @param size Height and width of the unit matrix.
    /// @return
    static TropicoMatrix GetI(size_t size);

    void Resize(size_t height, size_t width);

    /// @brief Gets a submatrix of the matrix.
    /// @param x First coordinate of the left up angle of a submatrix.
    /// @param y Second coordinate of the left up angle of a submatrix.
    /// @param height Height of a submatrix.
    /// @param width Width of a submatrix.
    /// @return Submatrix.
    TropicoMatrix SubMatrix(size_t x, size_t y, size_t height, size_t width) const;
    /// @brief Multiplicative-conjugate transposition.
    /// @return 
    TropicoMatrix MultiConjTrans() const;
    /// @brief Standardizes matrix columns.
    /// @return 
    TropicoMatrix& Standardization();

    /// @brief Converts to double.
    /// @return Double matrix.
    std::vector<std::vector<double>> ToDouble() const;

    const std::vector<TropicoFrac>& operator[](size_t i) const { return matr[i]; }
    std::vector<TropicoFrac>& operator[](size_t i) { return matr[i]; }

    TropicoMatrix& operator+=(const TropicoMatrix& tm);
    TropicoMatrix operator+(const TropicoMatrix& tm) const;
    TropicoMatrix& operator*=(const TropicoMatrix& tm);
    TropicoMatrix operator*(const TropicoMatrix& tm) const;
    TropicoMatrix& operator*=(const TropicoFrac& tf);
    TropicoMatrix operator*(const TropicoFrac& tf) const;

    /// @brief Trace of the matrix.
    /// @return 
    TropicoFrac Trace() const;
    /// @brief Spectral radius of the matrix.
    /// @return 
    TropicoFrac SpectralRadius(std::vector<TropicoMatrix>* degreeMatrix) const;
    /// @brief Kleene matrix of the matrix.
    /// @return 
    TropicoMatrix Kleene() const;

    /// @brief Removes correlated columns.
    /// @return 
    TropicoMatrix RemoveCorrel() const;

    friend bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2);

    friend std::ostream& ToDoubleTeXOut(std::ostream& out, const TropicoMatrix& tm);
    friend std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm);
    friend std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm);
    friend std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm);
    friend std::istream& operator>>(std::istream& in, TropicoMatrix& tm);
};

/// @brief Checks two vectors for correlation.
/// @param tm1 First vector.
/// @param tm2 Second vector.
/// @return true, if two vectors is correlated, false, else.
bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2);

inline TropicoMatrix operator*(const TropicoFrac& tf, const TropicoMatrix& tm) { return tm * tf; }

/// @brief Outputs a double matrix in TeX format.
/// @param out Output stream.
/// @param tm Matrix.
/// @return Output stream.
std::ostream& ToDoubleTeXOut(std::ostream& out, const TropicoMatrix& tm);
/// @brief Outputs a double matrix.
/// @param out Output stream.
/// @param tm Matrix.
/// @return Output stream.
std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm);
/// @brief Outputs a matrix in TeX format.
/// @param out Output stream.
/// @param tm Matrix.
/// @return Output stream.
std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm);
std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm);
std::istream& operator>>(std::istream& in, TropicoMatrix& tm);

/// @brief Class of solver of a single-criteria task.
class TropicoSolve
{
private:
    TropicoMatrix tm;
    std::vector<TropicoMatrix> degreeTM;
    TropicoFrac spectral;
    TropicoMatrix kleene;
    TropicoMatrix worstMatr;
    std::vector<TropicoMatrix> bestMatr;
    TropicoMatrix worst;
    std::vector<TropicoMatrix> best;
    TropicoFrac delta;
    TropicoMatrix P;
    std::vector<TropicoMatrix> Plk;

    std::string name;
public:
    /// @brief Constructor.
    /// @param _tm Matrix of a criteria.
    TropicoSolve(TropicoMatrix _tm = TropicoMatrix());

    /// @brief Finds the worst solution.
    /// @return Worst solution.
    TropicoMatrix& WorstSolve();
    /// @brief Makes P-matrix.
    /// @return P-matrix.
    TropicoMatrix MakeP() const;
    /// @brief Makes P_{lk}-matrix.
    /// @return P_{lk}-matrix.
    std::vector<TropicoMatrix> MakePlk() const;
    /// @brief Finds the best solution.
    /// @return Best solution.
    std::vector<TropicoMatrix>& BestSolve();

    /// @brief Gives the matrix a name.
    /// @param _name Name.
    inline void GiveName(std::string _name) { name = _name; }

    /// @brief Getter.
    /// @return Matrix of a criteria.
    inline const TropicoMatrix& GetMatrix() const { return tm; }
    /// @brief Getter.
    /// @return Degrees of matrix.
    inline const std::vector<TropicoMatrix>& GetDegreeMatrix() const { return degreeTM; }
    /// @brief Getter.
    /// @return Spectral radius of the matrix.
    inline const TropicoFrac& GetSpectral() const { return spectral; }
    /// @brief Getter.
    /// @return Kleene matrix of the matrix.
    inline const TropicoMatrix& GetKleene() const { return kleene; }
    /// @brief Getter.
    /// @return Matrix of the worst solution.
    inline const TropicoMatrix& GetWorstMatrix() const { return worstMatr; }
    /// @brief Getter.
    /// @return Vector of matrixes of the best solution.
    inline const std::vector<TropicoMatrix>& GetBestMatrix() const { return bestMatr; }
    /// @brief Getter.
    /// @return Worst solution.
    inline const TropicoMatrix& GetWorst() const { return worst; }
    /// @brief Getter.
    /// @return Vector of the best solution.
    inline const std::vector<TropicoMatrix>& GetBest() const { return best; }
    /// @brief Getter.
    /// @return Delta.
    inline const TropicoFrac& GetDelta() const { return delta; }
    /// @brief Getter.
    /// @return P-matrix.
    inline const TropicoMatrix& GetP() const { return P; }
    /// @brief Getter.
    /// @return P_{lk}-matrix.
    inline const std::vector<TropicoMatrix>& GetPlk() const { return Plk; }
    /// @brief Getter.
    /// @return Matrix name.
    inline const std::string& GetName() const { return name; }

    /// @brief Solves a single-criteria problem.
    void Solve();

    friend std::ostream& MakeTeX(std::ostream& out, const TropicoSolve& ts);
    friend std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts);
};

/// @brief Outputs a solve in TeX format.
/// @param out Output stream.
/// @param ts Solve.
/// @return Output stream.
std::ostream& MakeTeX(std::ostream& out, const TropicoSolve& ts);
std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts);

/// @brief Class of a solver multi-criteria task.
class TropicoMultiSolve
{
private:
    TropicoMatrix C;
    TropicoMatrix BWorst;
    std::vector<TropicoMatrix> BBest;
    TropicoSolve tsC, tsBW;
    std::vector<TropicoSolve> tsBB;
    std::vector<TropicoMatrix> A;
public:
    TropicoMultiSolve() = default;
    TropicoMultiSolve(TropicoMatrix& _C, std::vector<TropicoMatrix>& _A) : C(_C), A(_A) {}
    
    /// @brief Solves a multi-criteria problem.
    void Solve();

    friend std::ostream& MakeTeX(std::ostream& out, const TropicoMultiSolve& tms);
    friend std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms);
    friend std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms);
};

/// @brief Outputs a matrix of solve in TeX format.
/// @param out Output stream.
/// @param str Name of a matrix.
/// @param tm Matrix.
/// @return Output stream.
std::ostream& MakeTeXMatrix(std::ostream& out, std::string str, const TropicoMatrix& tm);
/// @brief Outputs a number with a name in TeX format.
/// @param out Output stream.
/// @param str Name of a number.
/// @param tf Number.
/// @return Output stream.
std::ostream& MakeTeXFrac(std::ostream& out, std::string str, const TropicoFrac& tf);

std::ostream& MakeTeX(std::ostream& out, const TropicoMultiSolve& tms);
std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms);
std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms);