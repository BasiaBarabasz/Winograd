#include <iostream>
#include <ctime>
#include <cstdlib>
#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "utils_win.h"

using namespace std;

//return division remainder
//Horner method
polynomial modulo_horner(polynomial a, polynomial b) {
  polynomial s(degree(a));
  polynomial r(1);

  s[degree(a) - 1] = a[degree(a)];
  for(int i = degree(a) - 2; i >= 0; i--) {
    Fraction t1;
    t1 = s[i + 1].minus();
    Fraction t2;
    t2 = b[0].mul(t1);
    s[i] = a[i + 1].add(t2);
  }
  Fraction temp1 = s[0].minus();
  Fraction temp2 = b[0].mul(temp1);
  r[0] = a[0].add(temp2);
  return r;
}

//return result of division a/b
//Horner method
polynomial divide_horner(polynomial a, polynomial b) {
  polynomial s(degree(a));
  
  s[degree(a) - 1] = a[degree(a)];
  for(int i = degree(a) - 2; i >= 0; i--) {
    Fraction t1;
    t1 = s[i + 1].minus();
    Fraction t2;
    t2 = b[0].mul(t1);
    s[i] = a[i + 1].add(t2);
  }
  return s;
}

//return coefficients a_0,..,a_n for a(q)
polynomialh value_pol(polynomial a, Fraction q) {
  polynomialh b(degree(a) + 1);
  Fraction temp = Fraction(1, 1);

  b[0] = a[0];
  for(int i = 1; i <= degree(a); i++) {
    temp = temp.mul(q);
    b[i] = temp.mul(a[i]);
  }
  return b;
} 

// create polynomial h_0,...,h_n of the degree equal to d (coefficients equal to one)
polynomialh create_pol(int d) {
  polynomialh a(d + 1);
  for(int i = 0; i <= d; i++)
    a[i] = Fraction(1, 1);
  return a;
}

// return root point of the linear polynomial
Fraction root_pol(polynomial a){
  Fraction b;

  b = a[0].minus();
  return b;
}

/*
void initialize() {
  //srand((unsigned)time(0));
  std::random_device rd;
  std::mt19937 gen(rd());
}

long long get_number(long long range1, long long range2) { //losowanie z przedzialu range1 - range2
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(range1,range2);
  return dis(gen);
}

Fraction get_fraction() { //ustalanie przedzialu losowania
  long long denominator = get_number(1, 10);
  long long numerator = get_number(1,4 * denominator) - 2 * denominator;
  return Fraction(numerator, denominator); 
}

int fraction_in_vector(Fraction f, vector<Fraction> v) { //zwraca 1 gdy punkt sie powtarza
  for(int i = 0; i < v.size(); i++)
    if(f.equal(v[i]))
      return 1;
  return 0;
}

Fraction get_fraction_no_repeat(vector<Fraction> except) { //losowanie roznych punktow
  Fraction f = get_fraction();
  while(fraction_in_vector(f, except))
    f = get_fraction();
  return f;
}
*/
int main(int argc, char* argv[]) {
  try {

    if(argc < 3) {
      cout << "Not enough arguments \n";
      return 1;
    }

    int Dh = atoi(argv[1]) - 1; //size of kernel (polynomial h)
    int Dx = atoi(argv[2]) - 1; // size of output (polynomial x)
    int K = Dh + Dx + 1; //degree of the polynomial m(p)

    if(argc < 2 * (K-1) + 3) {
      cout << "Not enough arguments";
      return 1;
    }

    polynomialh h(Dh); //polynomial h (kernel)
    polynomialh x(Dx); //polynomial x (output)
    matrix m_i(K, row(K)); //matrix of the linear polynomials m_0(p),..., m_k(p), m_i(p) = p + a_i
    polynomial m(K); // product of m_0(p),..., m_k(p)
    matrix M_i(K, row(K)); // matrix of the polynomials M_i(p),... M_i(p) (M_i(p)=m(p)/m_i(p))
    matrix r_i(K, row(K)); // matrix of thepolynomials M_i(p) mod m_i(p)
    matrixh Vh_i(K, rowh(K)); //matrix of the polynomials with coefficients h_0,...,h_n computed
    //as h(p) mod m_i(p) from Bezou theorem = h(a_i) where a_i is a root point of m_i(p)
    matrixh Vx_i(K, rowh(K)); // matrix of the polynomials with coefficients x_0,..., x_n computed
    //as x(p) mod m_i(p) from Bezou theorem = x(a_i) where a_i is a root point of m_i(p)
    //matrix_int S(SIZE, row_int(SIZE));
    matrixh B(K, rowh(K));

    h = create_pol(Dh);  //polynomial of the degree Dh
    x = create_pol(Dx);  //polynomial of the degree Dx

    for(int i = 0; i < K - 1; i++) {//load polynomials m_i(p) - coefficient as Fraction
      long long m = atoi(argv[3 + 2 * i]);
      long long n = atoi(argv[3 + 2 * i + 1]);
      m_i[i][0] = Fraction(-m, n);
      m_i[i][1] = Fraction(1, 1);
    }

    //product of the polynomials m_i
    m = m_i[0];
    for(int i = 1; i < K - 1; i++)
      m = multiplication_pol(m, m_i[i]);
    for(int i = 0; i < K - 1; i++){
      M_i[i] = divide_horner(m, m_i[i]);  // M_i(p)

      r_i[i] = modulo_horner(M_i[i], m_i[i]); // r_i for - M_i[i] mod m_i[i]
      Vh_i[i] = value_pol(h, root_pol(m_i[i])); //h_i(p) = h(p)modm_i(p) for m_i(p)=p+a
      Vx_i[i] = value_pol(x, root_pol(m_i[i])); //x_i(p) = x(p)modm_i(p) for m_i(p)=p+a
    }

    for(int i = 0; i < K - 1; i++) {
      for(int j = 0; j < K - 1; j++) {
        B[i][j] = M_i[i][j];
      }
    }

    for(int i = 0; i < K; i++) {
      B[i][K - 1] = Fraction(0, 1);
    }
    for(int j = 0; j < K; j++) {
      B[K - 1][j] = m[j];
    }

    ofstream myfile ("matrix_output/BTMat.txt");
    if (myfile.is_open()) {
      for(int i = 0; i < K; i++) {
        for(int j = 0; j < K; j++) {
          myfile << B[i][j].numerator << " " << B[i][j].denominator << "  ";
        }
        myfile << '\n';
      }
      myfile.close();
    }

    for(int i = 0; i < K - 1; i++) {
      for(int j = 0; j < Dh + 1; j++) {
        Vh_i[i][j] = Vh_i[i][j].divide(r_i[i][0]); //dzielenie przez r
      }
    }

    for(int j = 0; j < K; j++) {
      Vh_i[K - 1][j] = Fraction(0, 1);
    }
    Vh_i[K - 1][Dh] = Fraction(1, 1);

    ofstream myfile1 ("matrix_output/GMat.txt");
    if(myfile1.is_open()) {
      for(int i = 0; i < K; i++) {
        for(int j = 0; j < Dh + 1; j++) {
          myfile1 << Vh_i[i][j].numerator << " " << Vh_i[i][j].denominator << "  ";
        }
        myfile1 << '\n';
      }
      myfile1.close();
    }

    for(int j = 0; j < Dx + 1; j++) {
      Vx_i[K - 1][j] = Fraction(0, 1);
    }
    Vx_i[K - 1][Dx] = Fraction(1, 1);

    ofstream myfile2 ("matrix_output/ATMat.txt");
    if(myfile2.is_open()) {
      for(int j = 0; j < Dx + 1; j++) {
        for(int i = 0; i < K; i++) {
          myfile2 << Vx_i[i][j].numerator << " " << Vx_i[i][j].denominator << "  ";
        }
        myfile2 << '\n';
      }
      myfile2.close();
    }
  } catch (const char* msg) {
    std::cerr << msg << std::endl;
  }
  return 0;
}
