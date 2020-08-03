#include <iostream>
#include "utils_win.h"

using namespace std;

Fraction::Fraction() {
  numerator = 0;
  denominator = 1;
}

Fraction::Fraction(long long numer, long long denomin) {
  numerator = numer;
  denominator = denomin;
  reduce();
}

void Fraction::print() {
  cout << numerator << " "<< denominator;
}

Fraction Fraction::minus() {
  Fraction s;
  s.numerator = -numerator;
  if(-s.numerator != numerator) {
    throw "Overflow error occurred in Fraction::minus()";
  }
  s.denominator = denominator;
  return s;
}

int Fraction::low_zero() {
   long long t;
   t = numerator * denominator;
   if(t / denominator != numerator) {
     throw "Overflow error occurred in Fraction::low_zero()";
   }
   if(numerator*denominator < 0) return 1;
   else return 0;
}

Fraction Fraction::reduce() {
  Fraction s;
  long long a = numerator;
  long long b = denominator;
  long long c;

  if(numerator == 0) return s;
  if(a < 0) a = -a;
  if(b < 0) b = -b;
  if(b > a) {
    long long t = a;
    a = b;
    b = t;
  }
  while(b != 0) {
    c = a % b;
    a = b;
    b = c;
  }
  s.numerator = numerator / a;
  s.denominator = denominator / a;
  return s;
}

Fraction Fraction::mul(Fraction a) { //multiplication
  Fraction s;
  s.numerator = a.numerator * numerator;
  if(a.numerator != 0) {
    if(s.numerator / a.numerator != numerator) {
      throw "Overflow error occurred in Fraction::mul()";
    }
  }
  if(numerator != 0) {
    if(s.numerator / numerator != a.numerator) {
      throw "Overflow error occurred in Fraction::mul()";
    }
  }
  s.denominator = a.denominator * denominator;
  if (s.denominator / a.denominator != denominator) {
    throw "Overflow error occurred in Fraction::mul()";
  }
  s = s.reduce();
  return s;
}

int Fraction::equal(Fraction a) { //return 1 if equal a, 0 if not
  long long t1;
  long long t2;

  t1 = a.numerator * denominator;
  if(t1 / denominator != a.numerator){
    throw "Overflow error occurred in Fraction::equal()";
  }
  t2 = a.denominator * numerator;
  if(t2 / a.denominator != numerator) {
    throw "Overflow error occurred in Fraction::equal()";
  }
  return a.numerator * denominator == a.denominator * numerator;
}

Fraction Fraction::divide(Fraction a) {
   Fraction s;
  s.numerator = numerator * a.denominator;
  if (s.numerator / a.denominator != numerator) {
    throw "Overflow error occurred in Fraction::divide()";
  }
  s.denominator = denominator * a.numerator;
  if (s.denominator / denominator != a.numerator) {
    throw "Overflow error occurred in Fraction::divide()";
  }
  s = s.reduce();
  return s;
}

Fraction Fraction::add(Fraction a) { //addition
  Fraction s;

  long long left_numerator = a.numerator * denominator;
  if (left_numerator / denominator != a.numerator) {
    throw "Overflow error occurred in Fraction::add()";
  }
  long long right_numerator = a.denominator * numerator;
  if (right_numerator / a.denominator != numerator) {
    throw "Overflow error occurred in Fraction::add()";
  }
  s.numerator = left_numerator + right_numerator;
  if (s.numerator - left_numerator != right_numerator) {
    throw "Overflow error occurred in Fraction::add()";
  }
  s.denominator = a.denominator * denominator;
  if (s.denominator / a.denominator != denominator) {
    throw "Overflow error occurred in Fraction::add()";
  }
  s = s.reduce();
  return s;
}

symb_polynomial::symb_polynomial(int rows, int cols) {
  this->rows = rows;
  this->cols = cols;

  this->elements = vector<vector<Fraction> >(rows, vector<Fraction>(cols));

  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      elements[i][j] = Fraction(0,1);
    }
  }
}

void symb_polynomial::print() {
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      elements[i][j].print();
      cout << "  ";
    }
    cout << "\n";
  }
}

void symb_polynomial::set(int row, int col, Fraction element) {
  elements[row][col] = element;
}

Fraction symb_polynomial::get(int row, int col) {
  return elements[row][col];
}

int symb_polynomial::get_rows() {
  int r = rows;
  return r;
}

int symb_polynomial::get_cols() {
  int c = cols;
  return c;
}

int symb_polynomial::degree() {//zwraca stopien wielomianu symbolicznego w macierzy
 int d = get_rows()-1;
 Fraction z = Fraction(0,1);
 int x = 0; //znacznik czy wiersz jest zerowy, gdy nie bedzie zmieniamy na 1
 while(d >= 0) {
     for(int i = 0; i < get_cols(); i++) if(!get(d,i).equal(z)) x = 1;
     if(x == 1) return d;
     d--;
 }
 return d;
}

symb_polynomial symb_polynomial::copy(symb_polynomial A) { //wstawia macierz A od
//podanych wspolrzednych a,b
//  if((a+A.rows >= rows) || (b+A.cols >= cols))
  symb_polynomial result(rows,cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++)
      result.set(i,j,A.get(i,j));
  return result;
}

symb_polynomial symb_polynomial::cut_rows() {
  symb_polynomial result(degree()+1,cols);
  for(int i = 0; i < result.get_rows(); i++)
    for(int j = 0; j < result.get_cols(); j++)
      result.set(i,j,get(i,j));
  return result;
}

Matrix::Matrix(int rows, int cols) {
  this->rows = rows;
  this->cols = cols;

  this->elements = vector<vector<Fraction> >(rows, vector<Fraction>(cols));

  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      elements[i][j] = Fraction(0,1);
    }
  }
}

void Matrix::print() {
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      elements[i][j].print();
      cout << "  ";
    }
    cout << "\n";
  }
}

void Matrix::set(int row, int col, Fraction element) {
  elements[row][col] = element;
}

Fraction Matrix::get(int row, int col) {
  return elements[row][col];
}

int Matrix::get_rows() {
  int r = rows;
  return r;
}

int Matrix::get_cols() {
  int c = cols;
  return c;
}

Matrix Matrix::add(Matrix B) {
  Matrix result(rows, cols);
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      Fraction tmp = get(i, j).add(B.get(i, j));
      result.set(i, j, tmp);
    }
  }
  return result;
}

Matrix Matrix::multiplication(Matrix A) {
  if(cols != A.rows) exit(1);
  Matrix result(rows, A.cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < A.cols; j++) {
      Fraction tmp = Fraction(0,1);
      for(int k = 0; k < cols; k++)
        tmp = tmp.add(get(i,k).mul(A.get(k,j)));
      result.set(i, j, tmp);
  }
  return result;
}

Matrix Matrix::multiplication_pol(symb_polynomial A) {
  if(cols != A.rows) exit(1);
  Matrix result(rows, A.cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < A.cols; j++) {
      Fraction tmp = Fraction(0,1);
      for(int k = 0; k < cols; k++)
        tmp = tmp.add(get(i,k).mul(A.get(k,j)));
      result.set(i, j, tmp);
  }
  return result;
}

Matrix Matrix::multiplication_wise(Matrix A) {
  if((cols != A.cols) || (rows != A.rows)) exit(1);
  Matrix result(rows, cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++) {
      Fraction tmp = get(i,j).mul(A.get(i,j));
      result.set(i,j, tmp);
  }
    return result;
}

Matrix Matrix::put_pol(int a, int b, symb_polynomial A) { //wstawia macierz A od podanych wspolrzednych
//a,b
//  if((a+A.rows >= rows) || (b+A.cols >= cols))
  Matrix result(rows,cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++)
      result.set(i,j,get(i,j));
  for(int i = 0; i < A.get_rows(); i++)
    for(int j = 0; j < A.get_cols(); j++) {
      result.set(a+i,b+j, A.get(i,j));
  }
  return result;
}

Matrix Matrix::put_matrix(int a, int b, Matrix A) { //wstawia macierz A od podanych wspolrzednych
//a,b
//  if((a+A.rows >= rows) || (b+A.cols >= cols))
  Matrix result(rows,cols);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++)
      result.set(i,j,get(i,j));
  for(int i = 0; i < A.get_rows(); i++)
    for(int j = 0; j < A.get_cols(); j++) {
      result.set(a+i,b+j, A.get(i,j));
  }
  return result;
}

Matrix Matrix::transpose() {
  Matrix result(cols,rows);

  for(int i = 0; i < cols; i++)
    for(int j = 0; j < rows; j++)
      result.set(i,j,get(j,i));
  return result;
}

// zwraca stopien wielomianu
int degree(polynomial a) {
  int d = a.size() - 1;
  if(d == -1) return -1;
  Fraction z = Fraction(0, 1);
  while((a[d].equal(z)) && (d >= 0)) d--;
  return d;
}
// wypisuje wielomian a
void print(polynomial a) {
  for(int i = 0; i <= degree(a); i++) {
    a[i].print();
    cout << " ";
  }
  cout << std::endl;
}

// zwraca wielomian bedacy wynikiem
// mnozenia wielomianow a(p) i b(p)
polynomial multiplication_pol(polynomial a, polynomial b) {
  polynomial result(degree(a) + degree(b) + 1);
  for(int i = 0; i <= degree(b); i++) {
    for(int j = 0; j <= degree(a); j++) {
      Fraction tmp = result[i + j];
      tmp = tmp.add(a[j].mul(b[i]));
      result[i + j] = tmp;
    }
  }
  return result;
}
