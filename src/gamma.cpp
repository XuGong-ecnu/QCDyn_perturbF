#include "gamma.h"
#include <cmath>
double Gamma::GetPrecision() const
{
    return m_precision;
}
Gamma Gamma::SetPrecision(double precision)
{
    m_precision = precision;
    return *this;
}
double Gamma::GetMaxNum() const
{
    return m_maxNum;
}
Gamma Gamma::SetMaxNum(double maxNum)
{
    m_maxNum = maxNum;
    return *this;
}
double Gamma::GetS() const
{
    return m_s;
}
Gamma Gamma::SetS(double s)
{
    m_s = s;
    return *this;
}
double Gamma::GetX() const
{
    return m_x;
}
Gamma Gamma::SetX(double x)
{
    m_x = x;
    return *this;
}
double Gamma::Exec() const
{
    double ans = 0;
    for (double ptr = m_x + m_precision/2; ptr < m_maxNum; ptr += m_precision) {
        ans += static_cast<double>(std::pow<double, double>(ptr, m_s - 1) * std::exp(-ptr));
    }
    return ans * m_precision;
}
Gamma::Gamma(double s, double x) : m_s(s), m_x(x) {}
