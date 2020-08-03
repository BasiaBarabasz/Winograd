#include <vector>

using namespace std;

class Fraction {
  public:
    long long numerator;
    long long denominator;

    Fraction();
    Fraction(long long numer, long long denomin);

    void print();
    Fraction reduce();
    Fraction mul(Fraction a);
    Fraction add(Fraction a);
    int equal(Fraction a);
    int low_zero();
    Fraction minus();
    Fraction divide(Fraction a);
};

class symb_polynomial {
  public:
    int rows;
    int cols;
    vector<vector<Fraction> > elements;

    symb_polynomial(int rows, int columns);
    void print();
    void set(int row, int col, Fraction element);
    Fraction get(int row, int col);
    int get_rows();
    int get_cols();
    int degree();
    symb_polynomial copy(symb_polynomial A);
    symb_polynomial cut_rows();
};

class Matrix {
  public:
    int rows;
    int cols;
    vector<vector<Fraction> > elements;

    Matrix(int rows, int columns);
    void print();
    void set(int row, int col, Fraction element);
    Fraction get(int row, int col);
    int get_rows();
    int get_cols();
    Matrix add(Matrix B);
    Matrix multiplication(Matrix A);
    Matrix multiplication_pol(symb_polynomial A);
    Matrix multiplication_wise(Matrix A);
    Matrix put_pol(int a, int b, symb_polynomial A);
    Matrix put_matrix(int a, int b, Matrix A); //wstawia macierz A do macierzy lewy gorny rog w
    Matrix transpose(); //zwraca macierz transponowana
};

typedef vector<Fraction> polynomial;
typedef vector<Fraction> row;
typedef vector<Fraction> column;
typedef vector<row> matrix;
typedef vector<Fraction> polynomialh;// coefficients h_0,..,h_n
typedef vector<Fraction> rowh;
typedef vector<rowh> matrixh;
typedef vector<Fraction> row_int;
typedef vector<row_int> matrix_int; //macierz liczbowa faktoryzacji (1)
typedef vector<polynomial> set_of_polynomial; //wektor wielomianow

int degree(polynomial a);
void print(polynomial a);
polynomial multiplication_pol(polynomial a, polynomial b);
