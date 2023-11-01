#include <exception>
#include <sstream>
#include <map>
#include <iomanip>
#include <numeric>
#include "TropicoFrac.h"

void GcdFrac(int& n, int& m)
{
    int minNM = std::min(std::abs(n), m), maxDel = 1;

    for (int i = 2; i <= minNM; ++i)
        if (n % i == 0 && m % i == 0)
            maxDel = i;

    n /= maxDel;
    m /= maxDel;

    if (n == 0)
        m = 1;
}

void DegreeFrac::Reduction()
{
    if (m < 0)
    {
        m *= -1;
        n *= -1;
    }

    GcdFrac(n, m);
}

DegreeFrac& DegreeFrac::operator+=(const DegreeFrac& df)
{
    n = n * df.m + df.n * m;
    m *= df.m;

    Reduction();

    return *this;
}

DegreeFrac DegreeFrac::operator+(const DegreeFrac& df) const
{
    DegreeFrac res(*this);

    return res += df;
}

DegreeFrac& DegreeFrac::operator*=(const DegreeFrac& df)
{
    n *= df.n;
    m *= df.m;

    Reduction();

    return *this;
}

DegreeFrac DegreeFrac::operator*(const DegreeFrac& df) const
{
    DegreeFrac res(*this);

    return res *= df;
}

DegreeFrac DegreeFrac::operator-() const
{
    DegreeFrac res(*this);

    res.n *= -1;

    return res;
}

bool DegreeFrac::operator<=(const DegreeFrac& df) const
{
    return n * df.m <= df.n * m;
}

std::ostream& operator<<(std::ostream& out, const DegreeFrac& df)
{
    out << "(" << df.n;

    if (df.m != 1)
        out << "/" << df.m;
        
    out << ")";

    return out;
}

std::istream& operator>>(std::istream& in, DegreeFrac& df)
{
    size_t pos;
    std::string str;

    in >> str;

    try
    {
        pos = str.find('/');
        df.n = std::stoi(str.substr(0, pos));
        str = str.substr(pos + 1);

        if (pos != std::string::npos)
            df.m = std::stoi(str);
    }
    catch(const std::exception& e)
    {
        throw std::invalid_argument("Invalid format. Try: int[/int]. For example: 1/2.");
    }

    df.Reduction();

    return in;
}

PrimeDegreeFrac& PrimeDegreeFrac::operator*=(const PrimeDegreeFrac& pdf)
{
    df += pdf.df;

    return *this;
}

PrimeDegreeFrac PrimeDegreeFrac::operator*(const PrimeDegreeFrac& pdf) const
{
    PrimeDegreeFrac res(*this);

    return res *= pdf;
}

PrimeDegreeFrac& PrimeDegreeFrac::operator^=(const DegreeFrac& _df)
{
    df *= _df;

    return *this;
}

PrimeDegreeFrac PrimeDegreeFrac::operator^(const DegreeFrac& df) const
{
    PrimeDegreeFrac res(*this);

    return res^df;
}

bool PrimeDegreeFrac::operator<=(const PrimeDegreeFrac& pdf) const
{
    return n <= pdf.n;
}

TropicoFrac::TropicoFrac(int _n, int _m)
{
    if (_m == 0)
        throw std::invalid_argument("The denominator cannot be 0.");

    if (_n < 0 || _m < 0)
        throw std::invalid_argument("The number cannot be less 0.");

    GcdFrac(_n, _m);
    Convert(_n, 1);
    Convert(_m, -1);
}

bool TropicoFrac::ContainPrime(int _n) const
{
    for (const PrimeDegreeFrac& p : n)
    {
        if (p.GetPrime() == _n)
            return true;
    }

    return false;
}

DegreeFrac TropicoFrac::GetDegreePrime(int _n) const
{
    DegreeFrac df;

    for (const PrimeDegreeFrac& p : n)
    {
        if (p.GetPrime() == _n)
            return p.GetDegreeFrac();
    }

    throw std::out_of_range("It not have this prime number.");
}

void TropicoFrac::AddDegree(int prime, DegreeFrac df)
{
    if (prime == 1)
        return;

    if (ContainPrime(prime))
    {
        DegreeFrac EraseDf = GetDegreePrime(prime);
        df += EraseDf;
        n.erase(PrimeDegreeFrac(prime, EraseDf));
    }

    if (df != DegreeFrac(0, 1))
        n.insert(PrimeDegreeFrac(prime, df));
}

void TropicoFrac::Convert(int _n, const DegreeFrac& delta)
{
    if (_n == 0)
        AddDegree(0, 1);

    for (int i = 2; i <= _n; ++i)
        while (_n % i == 0)
        {
            AddDegree(i, delta);
            _n /= i;
        }
}

void TropicoFrac::PositiveDegree(TropicoFrac& tf1, TropicoFrac& tf2) const
{
    DegreeFrac df;

    for (const PrimeDegreeFrac& p : n)
    {
        df = p.GetDegreeFrac();

        if (df < DegreeFrac(0, 1))
        {
            tf1.AddDegree(p.GetPrime(), -df);
            tf2.AddDegree(p.GetPrime(), -df);
        }
    }
}

unsigned long long TropicoFrac::ValueIntPositive() const
{
    unsigned long long res = 1;

    if (ContainPrime(0))
        return 0;

    for (const PrimeDegreeFrac& p : n)
        for (int i{}; i < p.GetDegreeFrac().GetN(); ++i)
            res *= p.GetPrime();

    return res;
}

void SqrtProd(std::ostream& out, const std::pair<DegreeFrac, std::vector<int>> t, int sign = 1)
{
    if (t.first.GetM() != 1)
    {
        if (t.first.GetM() == 2)
            out << "\\sqrt";
        else
            out << "\\sqrt[" << t.first.GetM() << "]";
    }

    int prod = 1;

    for (const int& i : t.second)
        prod *= i;
    
    out << "{" << prod;

    if (t.first.GetN() != sign)
        out << "^{" << sign * t.first.GetN() << "}";

    out << "}";
}

double TropicoFrac::ToDouble() const
{
    double res = 1;

    for (const PrimeDegreeFrac& p : n)
        res *= p.ToDouble();

    return res;
}

std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tf)
{
    bool isFrac = false, positFrac = false;
    std::map<DegreeFrac, std::vector<int>> TeXPDF;

    for (const PrimeDegreeFrac& p : tf.n)
    {
        if (p.GetDegreeFrac() < 0)
            isFrac = true;
        
        TeXPDF[p.GetDegreeFrac()].push_back(p.GetPrime());
    }

    if (isFrac)
        out << "\\frac";

    out << "{";

    for (const auto& t : TeXPDF)
        if (t.first > 0)
        {
            positFrac = true;
            SqrtProd(out, t);
        }

    if (!positFrac)
        out << 1;

    out << "} {";

    for (const auto& t : TeXPDF)
        if (t.first < 0)
            SqrtProd(out, t, -1);

    out << "}";

    return out;
}

bool TropicoFrac::operator<=(const TropicoFrac& tf) const
{
    DegreeFrac df;
    TropicoFrac tf1(*this), tf2(tf), tf3;

    // std::cout << tf1 << " " << tf2 << "\n";

    this->PositiveDegree(tf1, tf2);
    tf.PositiveDegree(tf1, tf2);

    int m = 1;

    for (const PrimeDegreeFrac& p : n)
        m = std::lcm(m, p.GetDegreeFrac().GetM());

    for (const PrimeDegreeFrac& p : tf.n)
        m = std::lcm(m, p.GetDegreeFrac().GetM());

    tf1 ^= DegreeFrac(m);
    tf2 ^= DegreeFrac(m);
    tf3 = tf1;

    for (const PrimeDegreeFrac& p : tf3.n)
    {
        if (tf2.ContainPrime(p.GetPrime()))
        {
            df = tf1.GetDegreePrime(p.GetPrime()) <= tf2.GetDegreePrime(p.GetPrime()) ?
                 tf1.GetDegreePrime(p.GetPrime()) : tf2.GetDegreePrime(p.GetPrime());
            tf1 *= PrimeDegreeFrac(p.GetPrime(), -df);
            tf2 *= PrimeDegreeFrac(p.GetPrime(), -df);
        }
    }

    if ((tf1.ValueIntPositive() <= tf2.ValueIntPositive()) != (this->ToDouble() <= tf.ToDouble()))
    {
        std::cout << "ERROR!!!\n";
        std::cout << m << "\n";
        std::cout << tf1 << " " << tf2 << "\n";
        std::cout << tf1.ValueIntPositive() << " " << tf2.ValueIntPositive() << "\n";
        std::cout << tf1.ToDouble() << " " << tf2.ToDouble() << "\n";
        std::cout << *this << " " << tf << "\n";
        std::cout << this->ToDouble() << " " << tf.ToDouble() << "\n\n";
    }

    return tf1.ValueIntPositive() <= tf2.ValueIntPositive();
}

TropicoFrac& TropicoFrac::operator+=(const TropicoFrac& tf)
{
    *this = *this >= tf ? *this : tf;

    return *this;
}

TropicoFrac TropicoFrac::operator+(const TropicoFrac& tf) const
{
    TropicoFrac res(*this);

    return res += tf;
}

TropicoFrac& TropicoFrac::operator*=(const TropicoFrac& tf)
{
    TropicoFrac res(*this);

    if (tf.ContainPrime(0))
    {
        n.clear();
        n.insert(PrimeDegreeFrac(0, DegreeFrac()));

        return *this;
    }

    if (ContainPrime(0))
        return *this;

    for (const PrimeDegreeFrac& p : tf.n)
        res.AddDegree(p.GetPrime(), p.GetDegreeFrac());

    this->n = res.n;

    return *this;
}

TropicoFrac TropicoFrac::operator*(const TropicoFrac& tf) const
{
    TropicoFrac res(*this);

    return res *= tf;
}

TropicoFrac& TropicoFrac::operator^=(const DegreeFrac& df)
{
    std::set<PrimeDegreeFrac> res;

    for (const PrimeDegreeFrac& p : n)
        res.insert(PrimeDegreeFrac(p.GetPrime(), p.GetDegreeFrac() * df));
    
    n = res;

    return *this;
}

TropicoFrac TropicoFrac::operator^(const DegreeFrac& df) const
{
    TropicoFrac res(*this);

    return res ^= df;
}

std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf)
{
    bool first = true;

    if (!tf.n.size())
        out << 1;

    for (const PrimeDegreeFrac& p : tf.n)
    {
        if (first)
            first = false;
        else
            out << "*";

        out << p.GetPrime();

        if (p.GetDegreeFrac() != DegreeFrac(1, 1))
            out << "^" << p.GetDegreeFrac();
    }

    return out;
}

std::istream& operator>>(std::istream& in, TropicoFrac& tf)
{
    size_t pos;
    int n, m;
    DegreeFrac df;
    std::stringstream ss;
    std::string str, token;

    tf.n.clear();

    in >> str;

    try
    {
        if ((pos = str.find('/')) == std::string::npos && str.find('^') == std::string::npos)
        {
            n = std::stoi(str.substr(0));

            tf = TropicoFrac(n, 1);
        }
        else if (str.find('^') == std::string::npos)
        {
            n = std::stoi(str.substr(0, pos));
            m = std::stoi(str.substr(pos + 1));

            tf = TropicoFrac(n, m);
        }
        else
        {
        do
            {
                df = DegreeFrac();

                pos = str.find('*');
                token = str.substr(0, pos);

                if (pos != std::string::npos)
                    str = str.substr(pos + 1);
                else
                    str = "";

                pos = token.find('^');
                n = std::stoi(token.substr(0, pos));
                token = token.substr(pos + 1);
                ss.flush();

                if (pos != std::string::npos)
                {
                    ss.clear();
                    ss << token.substr(0);
                    ss >> df;
                }

                tf.AddDegree(n, df);
            } while (!str.empty());
        }
    }
    catch(const std::exception& e)
    {
        throw std::invalid_argument("Invalid format. Try: int[^int[/int]][*...*][int[^int[/int]]] or int/int. For example: 2^3/4*5^6/7.");
    }

    return in;
}

TropicoMatrix::TropicoMatrix(size_t height, size_t width, TropicoFrac elem)
{
    matr.resize(height, std::vector<TropicoFrac>(width, elem));
}

size_t TropicoMatrix::GetWidth() const
{
    if (GetHeight() == 0)
        return 0;
    else
        return matr[0].size();
}

TropicoMatrix TropicoMatrix::GetI(size_t size)
{
    TropicoMatrix res(size, size);

    for (size_t i{}; i < size; ++i)
        res[i][i] = 1;

    return res;
}

void TropicoMatrix::Resize(size_t height, size_t width)
{
    matr.resize(height);

    std::for_each(matr.begin(), matr.end(), [width](std::vector<TropicoFrac>& v)
    {
        v.resize(width, TropicoFrac());
    });
}

TropicoMatrix TropicoMatrix::SubMatrix(size_t x, size_t y, size_t height, size_t width) const
{
    if (x + height > GetHeight() || y + width > GetWidth())
        throw std::invalid_argument("SubMatrix must be less than Matrix.");

    TropicoMatrix sub(height, width);

    for (size_t i{}; i < height; ++i)
        for (size_t j{}; j < width; ++j)
            sub[i][j] = matr[i + x][j + y];

    return sub;
}

TropicoMatrix TropicoMatrix::MultiConjTrans() const
{
    TropicoMatrix MCT(GetWidth(), GetHeight());

    for (size_t i{}; i < MCT.GetHeight(); ++i)
        for (size_t j{}; j < MCT.GetWidth(); ++j)
            if (matr[j][i] != 0)
                MCT[i][j] = matr[j][i] ^ (-1);
            else
                MCT[i][j] = 0;

    return MCT;
}

void TropicoMatrix::Standardization()
{
    if (GetWidth() != 1)
        throw std::invalid_argument("It is not vector.");

    if (GetHeight() < 1)
        throw std::invalid_argument("Vector is empty.");

    TropicoFrac sum;

    for (size_t i{}; i < GetHeight(); ++i)
        sum += matr[i][0];
    
    *this *= sum ^ (-1);
}

std::vector<std::vector<double>> TropicoMatrix::ToDouble() const
{
    std::vector<std::vector<double>> res(GetHeight(), std::vector<double>(GetWidth()));

    for (size_t i{}; i < GetHeight(); ++i)
        for(size_t j{}; j < GetWidth(); ++j)
            res[i][j] = matr[i][j].ToDouble();

    return res;
}

std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm)
{
    std::vector<std::vector<double>> outMatrix = tm.ToDouble();

    for (size_t i{}; i < outMatrix.size(); ++i)
    {
        for(size_t j{}; j < outMatrix[0].size(); ++j)
            out << std::setprecision(4) << std::setw(8) << std::left << double(outMatrix[i][j]) << " ";

        out << "\n";
    }

    return out;
}

std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm)
{
    out << "\\begin{pmatrix}\n";

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        for(size_t j{}; j < tm.GetWidth(); ++j)
            MakeTeX(out,tm.matr[i][j]) << " & ";

        out << "\\\\\n";
    }

    out << "\\end{pmatrix}";

    return out;
}

bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2)
{
    if (tm1.GetWidth() != 1 || tm2.GetWidth() != 1)
        throw std::invalid_argument("It is not vector.");

    if (tm1.GetHeight() < 1 || tm2.GetHeight() < 1)
        throw std::invalid_argument("Vector is empty.");

    if (tm1.GetHeight() != tm2.GetHeight())
        throw std::invalid_argument("Size of vectors are not equal.");

    TropicoFrac coefProp(tm1[0][0] * (tm2[0][0] ^ (-1)));

    for (size_t i{}; i < tm1.GetHeight(); ++i)
        if (tm1[i][0] * (tm2[i][0] ^ (-1)) != coefProp)
            return false;

    return true;
}

TropicoMatrix& TropicoMatrix::operator+=(const TropicoMatrix& tm)
{
    if (GetHeight() != tm.GetHeight() || GetWidth() != tm.GetWidth())
        throw std::invalid_argument("Size of matrixes not equal.");

    for (size_t i{}; i < tm.GetHeight(); ++i)
        for(size_t j{}; j < tm.GetWidth(); ++j)
            matr[i][j] += tm[i][j];

    return *this;
}

TropicoMatrix TropicoMatrix::operator+(const TropicoMatrix& tm) const
{
    TropicoMatrix res(*this);

    return res += tm;
}

TropicoMatrix& TropicoMatrix::operator*=(const TropicoMatrix& tm)
{
    *this = *this * tm;

    return *this;
}

TropicoMatrix TropicoMatrix::operator*(const TropicoMatrix& tm) const
{
    if (GetWidth() != tm.GetHeight())
        throw std::invalid_argument("Width of first matrix and height of second matrix not equal.");

    TropicoMatrix res(GetHeight(), tm.GetWidth());

    for (size_t i{}; i < res.GetHeight(); ++i)
        for(size_t j{}; j < res.GetWidth(); ++j)
            for (size_t l{}; l < GetWidth(); ++l)
                res[i][j] += matr[i][l] * tm[l][j];
        
    return res;
}

TropicoMatrix& TropicoMatrix::operator*=(const TropicoFrac& tf)
{
    for (size_t i{}; i < GetHeight(); ++i)
        for(size_t j{}; j < GetWidth(); ++j)
            matr[i][j] *= tf;

    return *this;
}

TropicoMatrix TropicoMatrix::operator*(const TropicoFrac& tf) const
{
    TropicoMatrix res(*this);

    return res *= tf;
}

TropicoFrac TropicoMatrix::Trace() const
{
    if (GetWidth() != GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");

    TropicoFrac res;

    for (size_t i{}; i < GetHeight(); ++i)
        res += matr[i][i];

    return res;
}

TropicoFrac TropicoMatrix::SpectralRadius() const
{
    if (GetWidth() != GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");

    TropicoFrac res;
    TropicoMatrix C(*this);

    for (size_t i{}; i < GetHeight(); ++i, C *= *this)
        res += C.Trace() ^ DegreeFrac(1, i + 1);

    return res;
}

TropicoMatrix TropicoMatrix::Kleene() const
{
    if (GetWidth() != GetHeight())
        throw std::out_of_range("Width of matrix and height of matrix not equal.");

    TropicoMatrix res = TropicoMatrix::GetI(GetHeight());
    TropicoMatrix C(*this);

    for (size_t i{1}; i < GetHeight(); ++i, C *= *this)
        res += C;

    return res;
}

std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm)
{
    size_t maxLength = 0;
    std::string str;
    std::stringstream ss;
    std::vector<std::string> strMatrix;

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        for(size_t j{}; j < tm.GetWidth(); ++j)
        {
            ss.clear();
            ss << tm[i][j];
            std::getline(ss, str);

            maxLength = std::max(maxLength, str.length());
            
            strMatrix.push_back(str);
        }
    }

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        for(size_t j{}; j < tm.GetWidth(); ++j)
        {
            strMatrix[j + i * tm.GetWidth()].resize(maxLength, ' ');

            out << strMatrix[j + i * tm.GetWidth()] << " ";
        }
        
        out << std::endl;
    }

    return out;
}

std::istream& operator>>(std::istream& in, TropicoMatrix& tm)
{
    size_t height, width;

    in >> height >> width;

    tm.Resize(height, width);

    for (size_t i{}; i < tm.GetHeight(); ++i)
        for(size_t j{}; j < tm.GetWidth(); ++j)
            in >> tm[i][j];

    return in;
}

TropicoSolve::TropicoSolve(TropicoMatrix _tm) : tm(_tm)
{
    if (tm.GetWidth() != tm.GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");
}

TropicoMatrix TropicoSolve::WorstSolve()
{
    delta = (TropicoMatrix(1, tm.GetWidth(), 1) * kleene * TropicoMatrix(tm.GetHeight(), 1, 1))[0][0];
    worst = ((delta ^ (-1)) * TropicoMatrix(tm.GetHeight(), tm.GetWidth(), 1) + (spectral ^ (-1)) * tm).Kleene() * TropicoMatrix(tm.GetHeight(), 1, 1);

    worst.Standardization();

    return worst;
}

TropicoMatrix TropicoSolve::MakeP() const
{
    bool isCorrel = false;
    TropicoMatrix _P(kleene.SubMatrix(0, 0, kleene.GetHeight(), 1));

    for (size_t i{1}; i < kleene.GetWidth(); ++i)
    {
        isCorrel = false;

        for (size_t j{}; j < _P.GetWidth(); ++j)
            if (CorrelVector(kleene.SubMatrix(0, i, kleene.GetHeight(), 1), _P.SubMatrix(0, j, _P.GetHeight(), 1)))
                isCorrel = true;

        if (!isCorrel)
        {
            _P.Resize(_P.GetHeight(), _P.GetWidth() + 1);

            for (size_t j{}; j < _P.GetHeight(); ++j)
                _P[j][_P.GetWidth() - 1] = kleene[j][i];
        }
    }

    return _P;
}

TropicoMatrix TropicoSolve::MakePlk() const
{
    size_t k, l;
    TropicoFrac maxPj, value;
    TropicoMatrix _Plk(P.GetHeight(), P.GetWidth());

    for (size_t i{}; i < P.GetWidth(); ++i)
    {
        if ((value = (TropicoMatrix(1, tm.GetWidth(), 1) * (P.SubMatrix(0, i, P.GetHeight(), 1) * P.SubMatrix(0, i, P.GetHeight(), 1).MultiConjTrans()) *
                      TropicoMatrix(tm.GetHeight(), 1, 1))[0][0]) > maxPj)
        {
            maxPj = value;
            k = i;      
        }
    }

    maxPj = 0;

    for (size_t i{}; i < P.GetHeight(); ++i)
    {
        if ((value = P[i][k] ^ (-1)) > maxPj)
        {
            maxPj = value;
            l = i;
        }
    }

    _Plk[l][k] = P[l][k];

    return _Plk;
}

TropicoMatrix TropicoSolve::BestSolve()
{
    P = MakeP();
    Plk = MakePlk();

    best = P * (TropicoMatrix::GetI(P.GetWidth()) + Plk.MultiConjTrans() * P) * TropicoMatrix(P.GetWidth(), 1, 1);
    best.Standardization();

    return best;
}

void TropicoSolve::Solve()
{
    spectral = tm.SpectralRadius();
    kleene = ((spectral ^ (-1)) * tm).Kleene();

    WorstSolve();
    BestSolve();
}

std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts)
{
    out << "Matrix " << ts.GetName() << ":\n";
    out << ts.GetMatrix() << "~\n";
    ToDoubleOut(out, ts.GetMatrix()) << "\n";
    out << "Spectral Radius: " << ts.GetSpectral() << " ~ " << ts.GetSpectral().ToDouble() << "\n\n";
    out << "Kleene:\n" << ts.GetKleene() << "~\n";
    ToDoubleOut(out, ts.GetKleene()) << "\n";
    out << "Delta: \n" << ts.GetDelta() << " ~ " << ts.GetDelta().ToDouble() << "\n\n";
    out << "Worst Solve:\n" << ts.GetWorst() << "~\n";
    ToDoubleOut(out, ts.GetWorst()) << "\n";
    out << "P:\n" << ts.GetP() << "~\n";
    ToDoubleOut(out, ts.GetP()) << "\n";
    out << "P_{lk}:\n" << ts.GetPlk() << "~\n";
    ToDoubleOut(out, ts.GetPlk()) << "\n";
    out << "Best Solve:\n" << ts.GetBest() << "~\n";
    ToDoubleOut(out, ts.GetBest()) << "\n";

    return out;
}

void TropicoMultiSolve::Solve()
{
    tsC = C;

    tsC.Solve();

    for (size_t i{}; i < tsC.GetWorst().GetHeight(); ++i)
        BWorst += tsC.GetWorst()[i][0] * A[i];

    for (size_t i{}; i < tsC.GetBest().GetHeight(); ++i)
        BBest += tsC.GetBest()[i][0] * A[i];

    tsBW = BWorst;
    tsBB = BBest;

    tsBW.Solve();
    tsBB.Solve();

    tsC.GiveName("C");
    tsBW.GiveName("B_1");
    tsBB.GiveName("B_2");
}

std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms)
{
    out << tms.tsC << "\n";
    out << tms.tsBW << "\n";
    out << tms.tsBB << "\n";

    return out;
}

std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms)
{
    TropicoMatrix tm;
    in >> tms.C;

    for (size_t i{}; i < tms.C.GetHeight(); ++i)
    {
        in >> tm;
        tms.A.push_back(tm);
    }

    tms.BWorst.Resize(tms.A[0].GetHeight(), tms.A[0].GetWidth());
    tms.BBest.Resize(tms.BWorst.GetHeight(), tms.BWorst.GetWidth());

    return in;
}