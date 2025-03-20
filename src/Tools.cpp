/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Oct. 08, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Tools.h"

void Matrix_Transpose(const Real_Matrix& A, Real_Matrix& A_trans) {
    const int size = A.size();
    A_trans.resize(size, std::vector<double>(size, 0));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            A_trans[i][j] = A[j][i];
}

void Matrix_Multiply(const Real_Matrix& A, const Real_Matrix& B, Real_Matrix& C) {
    const int row_A = A.size();
    const int col_A = A[0].size();
    const int row_B = B.size();
    const int col_B = B[0].size();
    if (col_A != row_B)
        throw std::runtime_error("ERROR: Can't do A*B, since col_A != row_B.");
    C.resize(row_A, std::vector<double>(col_B, 0));
    for (int i = 0; i < row_A; ++i)
        for (int j = 0; j < col_B; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < col_A; ++k)
                C[i][j] += A[i][k] * B[k][j];
        }
}

void Matrix_Multiply(const Complex_Matrix& A, const Complex_Matrix& B, Complex_Matrix& C) {
    const int row_A = A.size();
    const int col_A = A[0].size();
    const int row_B = B.size();
    const int col_B = B[0].size();
    if (col_A != row_B)
        throw std::runtime_error("ERROR: Can't do A*B, since col_A != row_B.");
    C.resize(row_A, std::vector<Complex>(col_B, 0));
    for (int i = 0; i < row_A; ++i)
        for (int j = 0; j < col_B; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < col_A; ++k)
                C[i][j] += A[i][k] * B[k][j];
        }
}

void Matrix_Multiply(const Complex_Matrix& A, const std::vector<Complex>& B, std::vector<Complex>& C) {
    const int row_A = A.size();
    const int col_A = A[0].size();
    const int row_B = B.size();
    if (col_A != row_B)
        throw std::runtime_error("ERROR: Can't do A*B, since col_A != row_B.");
    C.resize(row_A, 0);
    for (int i = 0; i < row_A; ++i) {
        C[i] = 0;
        for (int j = 0; j < col_A; ++j)
            C[i] += A[i][j] * B[j];
    }
}

void Matrix_Add(const Real_Matrix& A, const Real_Matrix& B, Real_Matrix& C, double factor) {
    const int row_A = A.size();
    const int col_A = A[0].size();
    const int row_B = B.size();
    const int col_B = B[0].size();
    if (row_A != row_B || col_A != col_B)
        throw std::runtime_error("ERROR: Can't do A+B, since the number rows or columns of them are not same.");
    C = A;
    for (int i = 0; i < row_A; ++i)
        for (int j = 0; j < col_A; ++j)
            C[i][j] += factor * B[i][j];
}

void Matrix_Add(const Complex_Matrix& A, const Complex_Matrix& B, Complex_Matrix& C, double factor) {
    const int row_A = A.size();
    const int col_A = A[0].size();
    const int row_B = B.size();
    const int col_B = B[0].size();
    if (row_A != row_B || col_A != col_B)
        throw std::runtime_error("ERROR: Can't do A+B, since the number rows or columns of them are not same.");
    C = A;
    for (int i = 0; i < row_A; ++i)
        for (int j = 0; j < col_A; ++j)
            C[i][j] += factor * B[i][j];
}

void SplitString(std::vector<std::string>& list, const std::string& str, char dlim) {
    std::string field;
    std::istringstream s(str);
    while (getline(s, field, dlim))
        list.push_back(field);
}

void SplitString(std::vector<double>& values, const std::string& str, char dlim) {
    std::string field;
    std::istringstream s(str);
    while (getline(s, field, dlim))
        values.push_back(std::stod(field));
}

void SplitString(std::vector<int>& intergers, const std::string& str, char dlim) {
    std::string field;
    std::istringstream s(str);
    while (getline(s, field, dlim))
        intergers.push_back(std::stoi(field));
}

void GetIndexList(std::vector<int>& indices, const std::string& str) {
    if (str.empty())
        throw std::runtime_error("ERROR: The input string is empty.");
    // Check if the str contains unrecoginzed char.
    // ASCII number of char '0-9' is 48-57, and char ',' is 44, char '-' is 45
    // Ref: https://en.cppreference.com/w/cpp/language/ascii
    for (int i = 0; i < str.size(); ++i)
        if (str[i] < ',' ||  (str[i] > '-' && str[i] < '0') || str[i] > '9')
            throw std::runtime_error("ERROR: Unrecoginzed character in string.");
    // Split the string with delimiter ',' to a vector of string firstly.
    std::vector<std::string> list;
    SplitString(list, str, ',');
    // Get intergers of indices from string
    for (int i = 0; i < list.size(); ++i) {
        // Find char '-' in string list
        // Ref: https://en.cppreference.com/w/cpp/string/basic_string/find
        std::string::size_type pos = list[i].find('-');
        if (pos == std::string::npos) // no '-'
            indices.push_back(std::stoi(list[i]));
        else { // with '-'
            // Get the start and end index
            std::vector<int> start_end;
            SplitString(start_end, list[i], '-');
            if (start_end.size() != 2) // onle one '-' is allowed
                throw std::runtime_error("ERROR: More than one '-' in: " + list[i]);
            for (int j = start_end[0]; j <= start_end[1]; ++j)
                indices.push_back(j);
        }
    }
    // Sort indices from samll to large and remove the duplicated indices
    // Ref: https://blog.csdn.net/xiangxianghehe/article/details/90637998
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

double GetLineNumber(const std::string& file) {
    std::ifstream inpfile(file);
    if (!inpfile)
        throw std::runtime_error("ERROR: File not found: " + file);
    int lineNumber = 0;
    for (std::string line; getline(inpfile, line); lineNumber++);
    inpfile.close();
    return lineNumber;
}

FILE* CheckFile(const std::string& file , int count) {
    FILE* f = nullptr;
    if (count == 0) {
        // w: write only, the file is created if it does not exist. If a file with
        // the same name already exists, its contents are discarded and the file
        // is treated as a new empty file.
        f = fopen(file.c_str(), "w");
        if (f == nullptr)
            throw std::runtime_error("ERROR: Can't create the file: " + file);
    }
    else {
        // r: read only, the file must exist.
        // If this is not the first time to write, new data will be appended.
        // Check if file already exists firstly, if not exist, throw an exception.
        f = fopen(file.c_str(), "r");
        if(f == nullptr)
            throw std::runtime_error("ERROR: Can't open the file: " + file);
        fclose(f);
        // a: append only, the file is created if it does not exist.
        f = fopen(file.c_str(), "a");
    }
    return f;
}

void SaveDataFile(const std::string& file, const std::vector<std::string>& headers,
                  const std::vector<std::vector<double>*>& data, std::string type,
                  bool multi, std::string prefix, int count, bool multiHeaders,
                  bool addIndex, int precision) {
    // if count = 0, create a new file, else, append to this file.
    FILE* outfile = CheckFile(file, count);
    const int cols = headers.size(); // number of columns
    if (data.size() != cols)
        throw std::runtime_error("ERROR: The number of headers doesn't equal to the number of data.");
    const int size = data[0]->size(); // length of records
    if (type.empty())
        type = GetFileSuffix(file);
    if (multi && !prefix.empty()) // print a sinle line with prefix
        fprintf(outfile, "%s%-d\n", prefix.c_str(), count+1);
    // fields are separated by at least one space (fixed width)
    if (type == "dat") {
        if ((count == 0) || (multi && multiHeaders)) { // Print headers
            if (addIndex) // add index (start from 1) as 1st column
                fprintf(outfile, "%5s", "index");
            for (int i = 0; i < cols; ++i)
                fprintf(outfile, "%*s", precision+8, headers[i].c_str());
            fprintf(outfile, "\n");
        }
        for (int j = 0; j < size; ++j) { // print data
            fprintf(outfile, "%5d", j+1);
            for (int i = 0; i < cols; ++i)
                // You can't use data[i][j], since it is a pointer.
                fprintf(outfile, "%*.*g",precision+8, precision, data[i]->at(j));
            fprintf(outfile, "\n");
        }
    }
    // fields are separated by one comma (,).
    else if (type == "csv") {
        if ((count == 0) || (multi && multiHeaders)) { // Print headers
            if (addIndex)
                fprintf(outfile, "%s", "index");
            for (int i = 0; i < cols; ++i)
                fprintf(outfile, ",%s", headers[i].c_str());
            fprintf(outfile, "\n");
        }
        for (int j = 0; j < size; ++j) { // print data
            fprintf(outfile, "%d", j+1);
            for (int i = 0; i < cols; ++i)
                fprintf(outfile, ",%.*g", precision, data[i]->at(j));
            fprintf(outfile, "\n");
        }
    }
    else
        throw std::runtime_error("ERROR: Unknown file type for output: " + type);
    fclose(outfile);
}

void LoadDataFile(const std::string& file, std::vector<std::vector<double>*>& data, std::string type, int skip_line, int skip_column, int count) {
    std::ifstream inpfile(file);
    if (!inpfile)
        throw std::runtime_error("ERROR: Can't open input file: " + file);
    // According to the file extension name to decide the function used
    // to process each line of file
    std::function<void(std::vector<std::string>&, const std::string&)> split;
    if (type.empty())
        type = GetFileSuffix(file);
    if (type == "dat") // fields are separated by at least one space.
        // For overloading functions, we must specify the type of function pointer
        // e.g., int(*)(int) is the type of function pointer, which can be as parameter
        // for function, and int(int) is the type of function.
        split = (void(*)(std::vector<std::string>&, const std::string&))SplitLine;
    else if (type == "csv") // fields are separated by one comma (,).
        split = std::bind((void(*)(std::vector<std::string>&, const std::string&, char))SplitString, std::placeholders::_1, std::placeholders::_2, ',');
    else
        throw std::runtime_error("ERROR: Unknown file type for input: " + type);
    const int size = data[0]->size(); // length of records
    const int cols = data.size(); // number of columns
    for (int i = 0; i < cols; ++i)
        if (data[i]->size() != size || data[i]->size() == 0)
            throw std::runtime_error("ERROR: The size of vectors must be same in data and can't be 0.");
    const int shift = (size+skip_line) * count; // shift for the start line (multi-frame)
    int lineNumber = 0, index = 0;
    for (std::string line; getline(inpfile, line); ) {
        lineNumber++; // Let line number start from 1
        if (lineNumber <= (shift+skip_line)) continue;
        else if (lineNumber >= (shift+skip_line+1) && lineNumber <= (shift+size+skip_line)) {
            std::vector<std::string> fields;
            split(fields, line);
            if (fields.size() < cols+skip_column)
                throw std::runtime_error("ERROR: Too few columns at line " + std::to_string(lineNumber));
            for (int i = 0; i < cols; ++i)
                data[i]->at(index) = std::stod(fields[i+skip_column]);
            index++;
        }
        if (index == size) break;
    }
    if (index != size)
        throw std::runtime_error("ERROR: Too few records of data in " + file);
    inpfile.close();
}

double Integrate(const std::vector<double>& fx, double dx) {
    // Compute integral \int f(x)dx using Trapezoidal Rule:
    // int f(x) dx = dx/2 * (f[0] + 2f[1] + 2f[2] + ... + 2f[n-2] + f[n-1])
    const int size = fx.size();
    double integral = 0.5 * (fx[0] + fx[size-1]);
    for(int x = 1; x < size-1; ++x)
        integral += fx[x];
    return integral*dx;
}

double Delta(int a, int b) { //kronecker delta
    double d;
    if (a==b) d = 1;
    else d = 0;
    return d;
}

 void MatrixMultiplyPi(double A[][3], double B[][3], double C[][3]) { // 3*3 matrix multiplication

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            C[i][j]=0;
    }
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++)
            for (int j = 0; j < 3; j++) {      
                C[i][j] += A[i][k]*B[k][j];
            }
    }
}

void MatrixMultiplyMu(double A[][3], double B[3], double C[3]) {
    int i,j,k;
    for (i=0; i<3; i++) {
        C[i]=0;
    }
    for (i=0; i<3; i++) {
        for (k=0; k<3; k++)
            C[i] += A[i][k]*B[k];
    }
    return;
}
