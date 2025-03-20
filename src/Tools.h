/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Oct. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include <map>
#include <tuple>
#include <cmath>
#include <ctime>
#include <regex>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <typeinfo>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include "Eigen/Eigen"

/**
 * Definition of data types.
 */
typedef std::complex<double> Complex;
typedef std::vector<std::vector<Complex> > Complex_Matrix;
typedef std::vector<std::vector<std::vector<Complex> > > Complex_3D_Array;
typedef std::vector<std::vector<double> > Real_Matrix;
typedef std::vector<std::vector<std::vector<double> > > Real_3D_Array;

/**
 * Definition of constants.
 * Note that hbar is reduced Planck constant, hbar = h/(2pi) = 1.0545718e-34 J·s
 * If using atomic units, hbar = 1.
 */
const double hbar = 1.0;
const double pi = 3.14159265358979324;
const double pi_rad = 180;
const Complex I(0.0, 1.0);
// constants from  Python-OpenMM unit/constants.py (CODATA 2018)
const double NA   = 6.02214076e23; // Avogadro constant in 1/mol
const double kB   = 1.380649e-23;  // Boltzmann constant in J/K
const double Rbar = 8.31446261815324e-3; // gas constant = NA*kB in kJ/K/mol
// CustomForceField for Ewald Sum 
#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))
const double sqrpi = 1.7724538509055159;
const double e = 1.602176634e-19;// Elementary charge (C)
const double e0 = 8.8541878128e-12;// Permittivity of Vacuum (C^2/Jm)
const double nm = 1.0E-9;// nm   
const double ONE_4PI_EPS0 = e * e * NA / (4 * pi * e0 * nm * 1000 );// unit:KJ/mol
// For perturb
const double Fm2 = 4 * pi * e0 * 1e-30;
const double mol = 1 /NA;
// For RPMD
const double hbar_SI = 1.054571800e-34; // J.s
/**
 * Definition of constants for unit conversions.
 * Ref: https://en.wikipedia.org/wiki/Hartree
 * Ref: https://en.wikipedia.org/wiki/Hartree_atomic_units
 */
const double au2eV    = 27.211386246;      // 1 hartree (a.u.) = 27.211386246 eV
const double eV2au    = 0.0367493221;      // 1 eV = 0.0367493221 hartree (a.u.)
const double au2kcal  = 627.5094740631;    // 1 hartree (a.u.) = 627.5094740631 kcal/mol
const double kcal2au  = 0.00159360144;     // 1 kcal/mol = 0.00159360144 hartree (a.u.)
const double au2cm    = 219474.6313632;    // 1 hartree (a.u.) = 219474.6313632 cm-1 (energy)
const double cm2au    = 4.55633525e-6;     // 1 cm-1 (energy) = hc/λ = 4.55633525e-6 hartree (a.u.)
const double au2kj    = 2625.49963948;     // 1 hartree (a.u.) = 2625.49963948 kj/mol
const double kj2au    = 3.808798847e-4;    // 1 kj/mol = 3.808798847e-4 hartree (a.u.)
const double kj2eV    = 0.0103642697;      // 1 kj/mol = 0.0103642697 eV
const double kcal2kj  = 4.18400;           // 1 kcal/mol = 4.184 kj/mol
const double kcal2eV  = 0.0433641043;      // 1 kcal/mol = 0.0433641043 eV
const double eV2kcal  = 23.0605478;        // 1 eV = 23.0605478 kcal/mol
const double eV2kj    = 96.485332;         // 1 eV = 96.485332 kj/mol
const double au2kT    = 315775.024804;     // 1 au = 315775.024804 K (Kelvin)
const double kT2au    = 3.16681156e-6;     // kT (1 Kelvin) = 3.16681156e-6 a.u.
const double au2Hz    = 6.5796839205e12;   // 1 a.u. = 6.5796839205e12 Hz (s^-1)
const double Hz2au    = 1.5198298463-13;   // 1 Hz (s^-1) = 1.5198298463-13 eV
const double au2s     = 2.4188843266e-17;  // 1 a.u. = 2.4188843266e-17 s (time)
const double au2ps    = 2.4188843266e-5;   // 1 a.u. = 2.4188843266e-5 ps
const double au2fs    = 2.4188843266e-2;   // 1 a.u. = 2.4188843266e-2 fs
const double ps2au    = 41341.373335;      // 1 ps = 41.34137 a.u.
const double fs2au    = 41.341373335;      // 1 fs = 41.34137 a.u.
const double amu2kg   = 1.660539040e-27;   // 1 atomic mass unit = 1.660539040e-27 kg
const double amu2g    = 1.660539040e-24;   // 1 atomic mass unit = 1.660539040e-24 g
const double au2mdmass = amu2kg * NA; // // 1 atomic mass unit to  kg/mol
const double perturbunit = 1.11265e-40; // This is the Raman perturb force unit
/**
 * Get current system date and time (accurate to the seconds).
 *
 * @return    a string records current system date and time
 */
inline std::string CurrentTime() {
    time_t now = time(NULL);
    std::string time = ctime(&now);
    return time.substr(0, time.size()-1); // delete the "\n"
}
/**
 * Check a string is blank (including space, "\t", "\n") or not.
 *
 * @param str  the string to be checked
 * @return     return true, if it is empty or consisted by space (including "\t", "\n")
 */
inline bool IsBlankString(const std::string& str) {
    return str.length() ==
        std::count_if(str.begin(), str.end(), [](unsigned char ch){ return std::isspace(ch); });
}
/**
 * Converts the given string to lowercase.
 * Uppercase letters ABCDEFGHIJKLMNOPQRSTUVWXYZ -->
 * Lowercase letters abcdefghijklmnopqrstuvwxyz
 * Or unmodified chracters if no lowercase version is listed.
 *
 * @param str the string to be converted
 * @return    Lowercase version of the string without modification of the original one
 */
inline std::string StringToLower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char ch){ return tolower(ch); });
    return str;
}
/**
 * Converts the given string to uppercase.
 * Lowercase letters abcdefghijklmnopqrstuvwxyz -->
 * Uppercase letters ABCDEFGHIJKLMNOPQRSTUVWXYZ
 * Or unmodified chracters if no uppercase version is listed.
 *
 * @param str the string to be converted
 * @return    Uppercase version of the string without modification of the original one
 */
inline std::string StringToUpper(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char ch){ return toupper(ch); });
    return str;
}
/**
 * Separate the fields for a line to a vector of string (field delimiter is space).
 *
 * For example, a line "name    id    age", after using Splitline(fields, line),
 * the strings "name", "id" and "age" will be stored in fields[0], fields[1]
 * and fields[2], respectively.
 *
 * @param fields the vector of string to store the value after separated
 * @param line   the line want to be separated
 */
inline void SplitLine(std::vector<std::string>& fields, const std::string& line) {
    std::istringstream stream(line);
    for(std::string str; stream >> str; fields.emplace_back(str));
}
/**
 * Separate the fields for a line to a vector of double (field delimiter is space).
 *
 * The content of this line is a batch of numbers.
 *
 * @param values the vector of double to store the value after separated
 * @param line   the line want to be separated
 */
inline void SplitLine(std::vector<double>& values, const std::string& line) {
    std::istringstream stream(line);
    for(double value; stream >> value; values.emplace_back(value));
}
/**
 * Separate the fields for a line to a vector of int (field delimiter is space).
 *
 * The content of this line is a batch of integers.
 *
 * @param integers the vector of int to store the integers after separated
 * @param line     the line want to be separated
 */
inline void SplitLine(std::vector<int>& integers, const std::string& line) {
    std::istringstream stream(line);
    for(int integer; stream >> integer; integers.emplace_back(integer));
}
/**
 * Split a string with a given delimiter to a vector of string.
 *
 * For example, a string "id=123", after using SplitString(list, str, '='),
 * the strings "id" and "123" will be stored in list[0] and list[1], respectively.
 *
 * @param list   the vector of string to store the value after splitted
 * @param str    the string want to be splitted
 * @param dlim   the char used as delimiter
 */
void SplitString(std::vector<std::string>& list, const std::string& str, char dlim = ',');
/**
 * Split a string which inlcudes a serious of point numbers with a given delimiter
 * to a vector of double.
 *
 * @param values the vector of double to store the value after splitted
 * @param str    the string want to be splitted
 * @param dlim   the char used as delimiter
 */
void SplitString(std::vector<double>& values, const std::string& str, char dlim = ',');
/**
 * Split a string which inlcudes a serious of integers with a given delimiter
 * to a vector of int.
 *
 * @param integers the vector of int to store the integers after splitted
 * @param str      the string want to be splitted
 * @param dlim     the char used as delimiter
 */
void SplitString(std::vector<int>& integers, const std::string& str, char dlim = ',');
/**
 * Split a string which inlcudes a serious of indices (integers) with a delimiter
 * '-' and/or ',', for example "1-28, 30".
 *
 * The splitted indices list is in assending order. And, if a index is
 * given more than once, only one index is kept.
 *
 * @param indices  the vector of int to store the indices after splitted
 * @param str      the string want to be splitted
 */
void GetIndexList(std::vector<int>& indices, const std::string& str);
/**
 * Get the extension name of file.
 *
 * @param file   the file name
 * @return       the extension name of file (without .)
 */
inline std::string GetFileSuffix(const std::string& file) {
    std::vector<std::string> list;
    SplitString(list, file, '.');
    return list.back();
}
/**
 * Get the name of file without extension name.
 *
 * @param file   the file name
 * @return       the name of file (without .extension)
 */
inline std::string GetFilePrefix(const std::string& file) {
    std::vector<std::string> list;
    SplitString(list, file, '.');
    return file.substr(0, file.size()-list.back().size()-1);
}
/**
 * Get the number of lines of file.
 *
 * @param file   the file name
 * @return       number of lines
 */
double GetLineNumber(const std::string& file);
/**
 * Check the status of a file and return FILE* used as output by fprintf().
 * Remember to close it after output using fclose(FILE*);
 *
 * If the output is first time, a new file will be created (If an old file with
 * same name already exists, its contents are discarded).
 * If the output is not first time, new contents will be appended.
 *
 * @param file   file name
 * @param count  the time of output
 * @return       a C-style FILE* used to output by fprintf()
 */
FILE* CheckFile(const std::string& file, int count = 0);
/**
 * Save data (vector<double>) to .csv or .dat file.
 *
 * In .dat file, columns are separated by at least one space with fixed width.
 * In .csv file, columns are separated by one comma (,).
 *
 * If the output is first time, a new file will be created (If an old file with
 * same name already exists, its contents are discarded).
 * If the output is not first time, new contents will be appended (if multi=true,
 * which means this is multi-frame file).
 *
 * This function will add one column (index) to output file (1st column).
 *
 * The default values of parameters are used for the single-frame file.
 *
 * @param file          the file name to save data
 * @param headers       the headers of columns
 * @param data          the data for output, a vector of pointers to data
 * @param type          the type of file, cvs or dat. By default it is empty, and
 *                      type is decided according to the extension name of file.
 * @param multi         Multi-frame, true or false, (used for multi-frame file)
 *                      if true, new contents will be appended to old file.
 * @param prefix        the prefix to of multi-frame, e.g., "traj=". (for multi-frame)
 *                      if it is empty (default), then no prefix will be output
 * @param count         the time/frame of output (used for multi-frame file)
 * @param multiHeaders  if true, the headers will be appended every time for
 *                      multi-frame output, by default is false.
 * @param index         whether to add one column (index) to output file (1st column)
 *                      This is enables by default.
 * @param precision     the precision (number of significant figures), by default is 16
 */
void SaveDataFile(const std::string& file, const std::vector<std::string>& headers,
                  const std::vector<std::vector<double>*>& data, std::string type = "",
                  bool multi = false, std::string prefix = "", int count = 0,
                  bool multiHeaders = false, bool addIndex = true, int precision = 16);
/**
 * Load data from .csv or .dat file (e.g., created by SaveDataFile()) to vector<double>.
 *
 * Note the vector of data used to store the loading data should be resized before
 * call this function, which means that this function don't accept an empty vector.
 *
 * A multi-frame file is supported, but load one frame only.
 *
 * In .dat file, columns are separated by at least one space with fixed width.
 * In .csv file, columns are separated by one comma (,).
 *
 * The default values of parameters are used for the single-frame file.
 *
 * @param file        the file name
 * @param data        the data used to store, a vector of pointers to data
 * @param type        the type of file, cvs or dat. By default it is empty, and
 *                    type is decided according to the extension name of file.
 * @param skip_line   the number of lines to be skipped, e.g., the headers
 * @param skip_column the number of columns to be skipped, e.g., the index
 * @param count       the time/frame to be loaded (used for multi-frame file only)
 *                    by default, it is 0, which means only one frame.
 */
void LoadDataFile(const std::string& file, std::vector<std::vector<double>*>& data,
                  std::string type = "", int skip_line = 1, int skip_column = 1, int count = 0);
/**
 * Compute integral \int f(x)dx using Trapezoidal Rule:
 * int f(x) dx = dx/2 * (f[0] + 2 f[1] + 2 f[2] + ... + 2 f[n-2] + f[n-1]).
 *
 * @param fx    the data of function f(x)
 * @param dx    dx
 */
double Integrate(const std::vector<double>& fx, double dx);
// # The following functions are used to manipulate of matrix (vector of vector in C++).
/**
 * Provide an matrix A_trans to store the Transpose matrix square Matrix A.
 * The matrix A_trans will be resized and initialized, so it can be any or empty matrix.
 *
 * @param A          the square matrix to be transported
 * @param A_trans    the Transpose matrix of A
 */
void Matrix_Transpose(const Real_Matrix& A, Real_Matrix& A_trans);
/**
 * Provide an matrix C to store the results from matrix A multiply matrix B.
 * The matrix C will be resized and initialized, so it can be any or empty matrix.
 *
 * @param A    the matrix A
 * @param B    the matrix B
 * @param C    the matrix C to store result of A * B
 */
void Matrix_Multiply(const Real_Matrix& A, const Real_Matrix& B, Real_Matrix& C);
void Matrix_Multiply(const Complex_Matrix& A, const Complex_Matrix& B, Complex_Matrix& C);
/**
 * Provide a vector C to store the results from matrix A multiply vector B.
 * The vector C will be resized and initialized, so it can be any or empty matrix.
 *
 * @param A    the matrix A
 * @param B    the vector B
 * @param C    the vector C to store result of A * B
 */
void Matrix_Multiply(const Complex_Matrix& A, const std::vector<Complex>& B, std::vector<Complex>& C);
/**
 * Provide an matrix C to store the results from matrix A + factor*B.
 * The matrix C can be any or empty matrix.
 *
 * @param A      the matrix A
 * @param B      the matrix B
 * @param C      the matrix C to store result of A + factor*B
 * @param factor if factor = -1, do A - B, defalut is 1.
 */
void Matrix_Add(const Real_Matrix& A, const Real_Matrix& B, Real_Matrix& C, double factor = 1.0);
void Matrix_Add(const Complex_Matrix& A, const Complex_Matrix& B, Complex_Matrix& C, double factor = 1.0);
double Delta(int a, int b);
void MatrixMultiplyMu(double A[][3], double B[3], double C[3]);
void MatrixMultiplyPi(double A[][3], double B[][3], double C[][3]);