#ifndef GAMMAFUNC_GAMMA_H
#define GAMMAFUNC_GAMMA_H
class Gamma {
public:
    Gamma() = default;
    Gamma(double s, double x);
    ~Gamma() = default;
    [[nodiscard]] double GetPrecision() const;
    Gamma SetPrecision(double precision);
    [[nodiscard]] double GetMaxNum() const;
    Gamma SetMaxNum(double maxNum);
    [[nodiscard]] double GetS() const;
    Gamma SetS(double s);
    [[nodiscard]] double GetX() const;
    Gamma SetX(double x);
    double Exec() const;
private:
    double m_precision = 0.01;
    double m_maxNum = 1e5;
    double m_s = 0;
    double m_x = 0;
};
#endif //GAMMAFUNC_GAMMA_H

