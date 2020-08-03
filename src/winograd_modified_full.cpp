#include <iostream>
#include <ctime>
#include <cstdlib>
#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <random>
#include <sstream>
#include <string>
#include <stdio.h>
#include "utils_win.h"

using namespace std;

//funkcje dla wielomianow

// zwraca wynik dzielenia calkowitego a(p)/b(p)
// dla wektorow
polynomial divide_pol(polynomial a, polynomial b) {
  polynomial sol(degree(a) - degree(b) + 1);
  polynomial r(degree(a) + 1);
  r = a;
  int i;

  if(degree(b) == 0) {
    for(i = 0; i <= degree(a); i++) {
      sol[i] = a[i].divide(b[0]); //sol[i] = a[i]/b[0];
      r[i] = Fraction(0, 1);
    }
  }  else {
    while(degree(r) >= degree(b)) {
      int x = degree(r) - degree(b);
      sol[x] = r[degree(r)].divide(b[degree(b)]);
      for(i = 0; i <= degree(b); i++) {
        r[i + x] = r[i + x].add(b[i].mul(sol[x].minus()));
      }
    }
  }
  return sol;
}

// zwraca reszte dzielenia a(p)/b(p)
// dla wektorow
polynomial remain_pol(polynomial a, polynomial b) {
  polynomial sol(degree(a) - degree(b) + 1);
  polynomial r(degree(a) + 1);
  r = a;
  int i;

  if(degree(b) == 0) {
    for(i = 0; i <= degree(a); i++) {
      sol[i] = a[i].divide(b[0]);
      r[i] = Fraction(0, 1);
    }
  } else {
    while(degree(r) >= degree(b)) {
      int x = degree(r) - degree(b);
      sol[x] = r[degree(r)].divide(b[degree(b)]);
      for(i = 0; i <= degree(b); i++) {
        r[i + x] = r[i + x].add(b[i].mul(sol[x].minus()));
      }
    }
  }
  return r;
}

polynomial copy(polynomial x) {
  polynomial result(degree(x) + 1);

  for(int i = 0; i <= degree(x); i++) {
    result[i].numerator = x[i].numerator;
    result[i].denominator = x[i].denominator;
  }
  return result;
}
//koniec funkcji dla wielomianow


//funkcje dla wielomianow symbolicznych

// zwraca reszte z dzielenia
//wielomianu symbolicznego A(p) przez wielomianb(p) 
symb_polynomial divide_rem(symb_polynomial A, polynomial b) {
  if(A.get_rows() - 1 < degree(b)) return A;

  symb_polynomial solution(A.degree() - degree(b) + 1, A.get_cols());
  symb_polynomial remainder(A.get_rows(), A.get_cols());
  remainder = remainder.copy(A);

  if(degree(b) == 0) {
    for(int i = 0; i <= A.degree() - degree(b); i++) 
      for(int j = 0; j < A.get_cols(); j++) {
        solution.set(i,j, A.get(i,j).divide(b[0])); //sol[i] = a[i]/b[0];
      }
  } else {
    while(int degree_r = remainder.degree() >= degree(b)) {
       int x = remainder.degree() - degree(b); 
       for(int i = 0; i < solution.get_cols(); i++) {
         solution.set(x, i, remainder.get(remainder.degree(), i).divide(b[degree(b)]));
       }
       for(int i = 0; i <= degree(b); i++) {
         for(int j = 0; j < remainder.get_cols(); j++) {
           remainder.set(i + x, j, remainder.get(i + x, j).add(b[i].mul(solution.get(x, j).minus())));
         }
       }
    }
  }
  remainder = remainder.cut_rows();
  return remainder;  
}

polynomial minus_pol(polynomial a, polynomial b) {
  int d = degree(a) + 1;
  int low = degree(b) + 1;
  if(degree(b) > degree(a)) {
    d = degree(b) + 1;
    low = degree(a) + 1;
  }
  polynomial c(d);
  int i = 0;
  while(i < low) {
    // cout << b[i].denominator << endl;
    c[i] = a[i].add(b[i].minus());
    i++;
  }
  if(d == degree(a) + 1) {
    for(int i = low; i < d; i++) c[i] = a[i];
  } else {
    for(int i = low; i < d; i++) c[i] = b[i].minus();
  }
  return c;
}

symb_polynomial form_column(polynomial a) {
  symb_polynomial solution(degree(a) + 1, 1);

  for(int i = 0; i < degree(a)+1; i++) {
    solution.set(i, 0, a[i]);
  }
  return solution;
}

symb_polynomial mat(Matrix A) {
  symb_polynomial solution(A.get_rows(),A.get_cols());
  
  for(int i = 0; i < solution.get_rows(); i++) {
    for(int j = 0; j < solution.get_cols(); j++) {
      solution.set(i, j, A.get(i, j));
    }
  }
  return solution;
}

// zwraca N^i(p) - wspolczynnik dla
// kombinacji liniowej dla M_i(p)
polynomial euclidean_ext(polynomial a, polynomial b) {
  polynomial alpha1 = {Fraction(0, 1)};
  polynomial alpha2 = {Fraction(1, 1)};
  polynomial beta1 = {Fraction(1, 1)};
  polynomial beta2 = {Fraction(0, 1)};
  polynomial r(degree(b) - 1);
  polynomial alpha = {Fraction(0, 1)};
  polynomial beta = {Fraction(0, 1)};
  polynomial sol(abs(degree(a) - degree(b)) + 1);
  Fraction zero;
  zero.numerator = 0;
  zero.denominator = 1;

  int change = 0;
  if(degree(a) < degree(b)) {
    polynomial temp = a;
    a = b;
    b = temp;
    change = 1;
  }
  while((degree(b) > 0) || (!b[0].equal(zero))) {
    sol = divide_pol(a, b);
    r = remain_pol(a, b);
    alpha = minus_pol(alpha2, multiplication_pol(sol, alpha1));
    beta = minus_pol(beta2, multiplication_pol(sol, beta1));
    a = b;
    b = r;
    alpha2 = alpha1;
    alpha1 = alpha;
    beta2 = beta1;
    beta1 = beta;
  }
  alpha = alpha2;
  beta = beta2; 
  // print(alpha);
  // cout << '\n';
  // print(a);
  if(change == 1) alpha = beta;
  alpha = divide_pol(alpha, a);
  return alpha;
}

Matrix read_matrix(char* path) {
  ifstream file;
  file.open(path);
  string line;
  int rows = 0;
  int cols = 0;
  while(getline(file, line)) {
    if(cols == 0) {
      istringstream iss(line);
      while (iss) {
        string sub;
        iss >> sub;
        if(!sub.empty()) cols++;
      }
    }
    rows++;
  }
  file.clear();
  file.seekg(0, ios::beg);
  Matrix R(rows, cols / 2);
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols / 2; j++) {
      long long num;
      long long denom;
      file >> num;
      file >> denom;
      R.set(i, j, Fraction(num, denom));
    }
  }
  file.close();
  return R;
}

void save_matrix(Matrix M, std::string filename) {
  ofstream myfile1 (filename);
  if (myfile1.is_open()) {
    for(int i = 0; i < M.get_rows(); i++) {
      for(int j = 0; j < M.get_cols(); j++) {
        Fraction element = M.get(i, j);
        myfile1 << element.numerator << " " << element.denominator << "  ";
      }
      myfile1 << '\n';
    }
    myfile1.close();
  }
}


int main() {

  //macierze dla wielomianow p, p - 1, p + 1
  //Matrix TC_h(3, 2);
  //Matrix TC_x(3, 2);
  //symb_polynomial TC_inv(3, 3);

  char *g_mat = (char*) "matrix_input/GMat_%d.txt";
  char *at_mat = (char*) "matrix_input/ATMat_%d.txt";
  char *bt_mat = (char*) "matrix_input/BTMat_%d.txt";
  char path[100];

  int h_size = 3;
  int x_size = 6; //input size (output size)
  int h_deg = h_size - 1;
  int x_deg = x_size - 1;
  int K = h_deg + x_deg + 1;
  int W = 6; //liczba wielomianow
  int NL2 = 1; //liczba wielomianow nie liniowych
  int NL4 = 0;
  int NL6 = 0;
  int NL8 = 0;
  int NL10 = 0;
  //int Maxd = 5;//maksymalny stopien wielomianu faktoryzacji +1

  cout << "zmienne" << std::endl;
  // wielomian dla wyjscia
  symb_polynomial X(x_size, x_size);
  // wielomian dla filtru
  symb_polynomial H(h_size, h_size);

  for(int i = 0; i < h_size; i++) {
    H.set(i, i, Fraction(1, 1));
  }

  for(int i = 0; i < x_size; i++) {
    X.set(i, i, Fraction(1, 1));
  }

  // macierz wielomianow faktoryzacji m_i(p)
  set_of_polynomial m_i(W, polynomial(K + 1));
  // macierz wielomianow M(p)/m_i(p)
  set_of_polynomial M_i(W, polynomial(K + 1));
  // wielomian M(p) (iloczyn m_i(p))
  polynomial M(K + 1);

  // przykladowa faktoryzacja
  // f(x) = x
  polynomial t1 = {Fraction(0, 1), Fraction(1, 1)};

  // f(x) = x + 1/4
  // polynomial t2 = {Fraction(1, 4), Fraction(1, 1)};

  // f(x) = x + 1/2
  polynomial t3 = {Fraction(1, 2), Fraction(1, 1)};

  // f(x) = x - 1/2
  polynomial t4 = {Fraction(-1, 2), Fraction(1, 1)};

  // f(x) = x - 1
  // polynomial t5 = {Fraction(-1, 1), Fraction(1, 1)};
  // f(x) = x + 1
  // polynomial t6 = {Fraction(1, 1), Fraction(1, 1)};

  // f(x) = x + 2
  polynomial t7 = {Fraction(2, 1), Fraction(1, 1)};

  // f(x) = x - 2
  polynomial t8 = {Fraction(-2, 1), Fraction(1, 1)};

  // f(x) = x - 3
  // polynomial t2 = {Fraction(-3, 1), Fraction(1, 1)};
  // f(x) = x - 1/4
  // polynomial t9 = {Fraction(-1, 4), Fraction(1, 1)};
  // f(x) = x
  // polynomial t10 = {Fraction(0, 1), Fraction(1, 1)};
  // f(x) = x^2
  // polynomial t10 = {Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)};
  // f(x) = x^2 - 1
  // polynomial t12 = {Fraction(-1, 1), Fraction(0, 1), Fraction(1, 1)};

  // f(x) = x^2 + 1
  polynomial t11 = {Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)};

  // f(x) = x^2 + x + 1
  // polynomial t13 = {Fraction(1, 1), Fraction(1, 1), Fraction(1, 1)};
  // f(x) = x^2 - x + 1
  // polynomial t14 = {Fraction(1, 1), Fraction(-1, 1), Fraction(1, 1)};
  // f(x) = x^2 + x
  // polynomial t14 = {Fraction(0, 1), Fraction(1, 1), Fraction(1, 1)};
  // f(x) = x^2 - x - 1
  // polynomial t15 = {Fraction(-1, 1), Fraction(-1, 1), Fraction(1, 1)};
  // f(x) = x^2 + x - 1
  // polynomial t16 = {Fraction(-1, 1), Fraction(1, 1), Fraction(1, 1)};
  // f(x) = x^2 - x
  // polynomial t15 = {Fraction(0, 1), Fraction(-1, 1), Fraction(1, 1)};
  // f(x) = x^4 + 1
  // polynomial t14 = {Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(1, 1)};
  // f(x) = x^4 - x^2 + 1
  // polynomial t15 = {Fraction(1, 1), Fraction(0, 1), Fraction(-1, 1),
  //                   Fraction(0, 1), Fraction(1, 1)};
  // f(x) = x^6 + 1
  // polynomial t16 = {Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(1, 1)};
  // f(x) = x^8 + 1
  // polynomial t12 = {Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)};
  // f(x) = x^10 + 1
  // polynomial t13 = {Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0, 1), Fraction(0, 1), Fraction(0, 1),
  //                   Fraction(0,1), Fraction(1, 1)};
  // f(x) = x + 1/4
  // polynomial t17 = {Fraction(1, 4), Fraction(1, 1)};
  // f(x) = x + 3/4
  // polynomial t18 = {Fraction(3, 4), Fraction(1, 1)};
  // f(x) = x - 4/3
  // polynomial t19 = {Fraction(-4, 3), Fraction(1, 1)};

  m_i[0]= copy(t1);
  m_i[1]= copy(t3);
  m_i[2]= copy(t4);
  m_i[3]= copy(t7);
  m_i[4]= copy(t8);
  m_i[5]= copy(t11);
  // m_i[6]= copy(t11);
  // m_i[7] = copy(t10);
  // m_i[8] = copy(t10);
  // m_i[9] = copy(t1);
  // m_i[10] = copy(t1);
  //koniec przykladowej faktoryzacji

  M = copy(m_i[0]);
  for(int i = 1; i < W; i++) {
    M = multiplication_pol(M, m_i[i]);
  }

  //print(M);
  //cout << endl;
  for(int i = 0; i < W; i++){
    M_i[i] = divide_pol(M, m_i[i]);
  }

  // rozmiary macierzy uzaleznic od roznych stopninie tylko kwadratowego
  // int add = 0;
  // for(int i = 0; i < W; i++) {

  Matrix E(W + NL2 + (3 * NL4) + (5 * NL6) + (7 * NL8) + (9 * NL10) + 1,
           W + NL2 + (3 * NL4) + (5 * NL6) + (7 * NL8) + (9 * NL10)); //
  Matrix Mod_h(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (13 * NL10) + 1,
               h_size);
  Matrix Mod_x(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (11 * NL10) + 1,
               x_size);
  // Matrix Mod_x(W + (2 * NL2) + (5 * NL4) + 2, x_size);
  Matrix C(W + NL2 + (3 * NL4) + (5 * NL6) + (7 * NL8) + (9 * NL10),
           W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (11 * NL10) + 1);//dla kwadr ok

  int id = 2; //numer zestawu macierzy (w katalogu)
  //macierz G
  int it1 = 0;
  int place1 = 0;
  while(it1 < W) {
    symb_polynomial Th = divide_rem(H, m_i[it1]);
    if(degree(m_i[it1]) == 1) {
      Mod_h = Mod_h.put_pol(place1, 0, Th);
      place1++;
    } else {
      id = degree(m_i[it1]);
      sprintf(path, g_mat, id);
      Matrix TC_h = read_matrix(path);
      Mod_h = Mod_h.put_matrix(place1, 0, TC_h.multiplication_pol(Th));
      place1 = place1 + TC_h.get_rows();
      // id++;
    }
    it1++;
  }
  // cout << Mod_h.get_rows();
  // cout << endl;

  for(int i = 0; i < h_size; i++) {
    Mod_h.set(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (13 * NL10),
              i,
              Fraction(0, 1));
  }
  Mod_h.set(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (13 * NL10),
            h_size - 1,
            Fraction(1, 1));
  // Mod_h.print();
  // cout << endl;

  id = 2;
  //macierz A
  int it2 = 0;
  place1 = 0;
  while(it2 < W) {
    symb_polynomial Tx = divide_rem(X, m_i[it2]);
    if(degree(m_i[it2]) == 1) {
      Mod_x = Mod_x.put_pol(place1, 0, Tx);
      place1++;
    } else {
      id = degree(m_i[it2]);
      sprintf(path, at_mat, id);
      Matrix TC = read_matrix(path);
      Matrix TC_x = TC.transpose();
      // TC_x.print();
      // cout << endl;
      // Tx.print();
      // cout << endl;
      Mod_x = Mod_x.put_matrix(place1, 0, TC_x.multiplication_pol(Tx));
      place1 = place1 + TC_x.get_rows();
      // id++;
    }
    it2++;
  }

  for(int i = 0; i < x_size; i++) {
    Mod_x.set(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (13 * NL10),
              i,
              Fraction(0, 1));
  }
  Mod_x.set(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (13 * NL10),
            x_size - 1,
            Fraction(1, 1));

  // Mod_x.print();
  // cout << endl;

  id = 2;
  //macierz C
  int it3 = 0;
  place1 = 0;
  int place2 = 0;
  while(it3 < W) {
    if(degree(m_i[it3]) == 1) {
      C.set(place1, place2, Fraction(1, 1));
      place1++;
      place2++;
    } else {
      id = degree(m_i[it3]);
      sprintf(path, bt_mat, id);
      Matrix TC_in = read_matrix(path);
      Matrix TC = TC_in.transpose();
      symb_polynomial TC_inv = mat(TC);
      C = C.put_pol(place1, place2, divide_rem(TC_inv, m_i[it3]));
      place1 = place1 + divide_rem(TC_inv, m_i[it3]).get_rows();
      place2 = place2 + divide_rem(TC_inv, m_i[it3]).get_cols();
      // id++;
    }
    it3++;
  }
  // C.print();
  // cout << endl;

  int it4 = 0;
  place2 = 0;
  while(it4 < W) {
    for(int i = 0; i < degree(m_i[it4]); i++) {
      polynomial v = euclidean_ext(M_i[it4], m_i[it4]);
      // print(v);
      // cout << endl;
      polynomial temp(i + 1);
      temp[i] = Fraction(1, 1);
      v = multiplication_pol(temp, v);
      polynomial w = multiplication_pol(v, M_i[it4]);
      if(degree(w) >= degree(M)) w = remain_pol(w, M);
      symb_polynomial z = form_column(w);

      // z.print();
      // cout << endl << endl;
      E = E.put_pol(0, place2, z);
      place2++;
    }
    it4++;
  } 
  // E.print();
  // cout << endl;

  Matrix B = E.multiplication(C);
  // B.print();
  // cout << endl;
  Matrix BT = B.transpose();
  // cout << BT.get_rows() << " " << BT.get_cols() << endl;
  for(int i = 0; i < W + NL2 + (3 * NL4) + (5 * NL6) + (7 * NL8) + (9 * NL10) + 1; i++) {
    BT.set(W + (2 * NL2) + (5 * NL4) + (7 * NL6) + (9 * NL8) + (11 * NL10), i, M[i]);
  }
  // BT.print();
  // cout << endl;
  Matrix G = Mod_h;
  // G.print();
  // cout << endl;
  Matrix AT = Mod_x.transpose();
  // AT.print();
  // cout << endl;

  cout << "  Saving transformation matrices" << std::endl;
  save_matrix(BT, "matrix_output/BTMat_Win_6_sym.txt");
  save_matrix(G, "matrix_output/GMat_Win_6_sym.txt");
  save_matrix(AT, "matrix_output/ATMat_Win_6_sym.txt");

  std::cout << "Run test:" << std::endl;
  // test
  symb_polynomial HT(h_size, 1);
  HT.set(0, 0, Fraction(1, 1));
  HT.set(1, 0, Fraction(2, 1));
  HT.set(2, 0, Fraction(1, 1));
  // HT.print();
  // cout << endl;

  symb_polynomial I(K, 1);
  I.set(0, 0, Fraction(1, 1));
  I.set(1, 0, Fraction(1, 1));
  I.set(2, 0, Fraction(1, 1));
  I.set(3, 0, Fraction(2, 1));
  I.set(4, 0, Fraction(2, 1));
  // I.print();
  // cout << endl;
  I.set(5, 0, Fraction(1, 1));
  I.set(6, 0, Fraction(2, 1));
  I.set(7, 0, Fraction(1, 1));
  // I.set(8, 0, Fraction(1, 1));
  // I.set(9, 0, Fraction(2, 1));
  // I.set(10, 0, Fraction(1, 1));
  // I.set(11, 0, Fraction(2, 1));
  // I.set(12, 0, Fraction(2, 1));
  // I.set(13, 0, Fraction(1, 1));
  // I.set(14, 0, Fraction(1, 1));

  Matrix test1 = G.multiplication_pol(HT);
  // test1.print();
  // cout << endl;
  Matrix test2 = BT.multiplication_pol(I);
  // test2.print();
  // cout << endl;
  Matrix test3 = test1.multiplication_wise(test2);
  // test3.print();
  // cout << endl;
  Matrix test4 = AT.multiplication(test3);

  test4.print();

  return 0;
}
