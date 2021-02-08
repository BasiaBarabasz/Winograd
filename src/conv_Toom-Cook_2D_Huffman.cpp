#include <iostream>
#include <fstream>
#include <random>
#include "utils.h"

using namespace std;

/*
 * Performs summation of left row using Huffman tree
 * template function
 *
 * @oaram a row
 * @param n size of a row
 * @param b row
 * @return result of a summation
 */
template<typename T>
T multiplication_left_huffman(Row<T> a, int n, Row<T> b) {
  Queue<T> Q;
  for(int j = 0; j < n; j++) {
    double e = log2(abs(a[j]));
    Tree<T>* t = new Tree<T>(e, a[j] * b[j]);
    Q.insert(t);
  }

  QueueNode<T>* p = Q.first();
  while(p->next != Q.last()) {
    Tree<T>* t = new Tree<T>(0, 0);
    t = merge(p->next->root, p->next->next->root);
    Q.qdelete();
    Q.qdelete();
    Q.insert(t);
  }
  return p->next->root->root->result();
}

/*
 * Performs summation of right row using Huffman tree
 * template function
 *
 * @oaram a row
 * @param n size of a row
 * @param b row
 * @return result of a summation
 */
template<typename T>
T multiplication_right_huffman(Row<T> a, int n, Row<T> b) {
  Queue<T> Q;
  for(int j = 0; j < n; j++) {
    double e = log2(abs(b[j]));
    Tree<T>* t = new Tree<T>(e, a[j] * b[j]);
    Q.insert(t);
  }

  QueueNode<T>* p = Q.first();
  while(p->next != Q.last()) {
    Tree<T>* t = new Tree<T>(0, 0);
    t = merge(p->next->root, p->next->next->root);
    Q.qdelete();
    Q.qdelete();
    Q.insert(t);
  }
  return p->next->root->root->result();
}

template<typename T>
Matrix<T> multiplication_matrix_left_huffman(Matrix<T> A, Matrix<T> B, int n, int m, int l) {
  Matrix<T> result(n, Row<T>(m));

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
       Row<T> a1(l);
       Row<T> a2(l);
       for(int it = 0; it < l; it++) {
         a1[it] = A[i][it];
         a2[it] = B[it][j];
       }
       result[i][j] = multiplication_left_huffman<T>(a1, l, a2);
    }
  }
  return result;  
}

template<typename T>
Matrix<T> multiplication_matrix_right_huffman(Matrix<T> A, Matrix<T> B, int n, int m, int l) {
  Matrix<T> result(n, Row<T>(m));

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      Row<T> a1(l);
      Row<T> a2(l);
      for(int it = 0; it < l; it++) {
        a1[it] = A[i][it];
        a2[it] = B[it][j];
      }
      result[i][j] = multiplication_right_huffman<T>(a1, l, a2);
    }
  }
  return result;
}

template<typename T>
Matrix<T> multiplication_matrix(Matrix<T> A, Matrix<T> B, int n, int m, int l) {
  Matrix<T> result(n, Row<T>(m));

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      for(int k = 0; k < l; k++) {
        result[i][j] = result[i][j] + (A[i][k] * B[k][j]);
      }
    }
  }
  return result;
}

/*
 * Performs element wise summation of matrices
 * template function
 *
 * @param A matrix
 * @param B matrix
 * @param n number of rows
 * @param m number of columns
 * @return matrix with sums
 */
template<typename T>
Matrix<T> multiplication_wise(Matrix<T> A, Matrix<T> B, int n, int m) {
  Matrix<T> result(n, Row<T>(m));

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      result[i][j] = (T)((float)A[i][j] * (float)B[i][j]);
    }
  }
  return result;
}

/*
 * Performs matrix transposition
 * template function
 *
 * @param A matrix
 * @param n number of rows
 * @param m number of columns
 * @return transposed matrix
 */
template<typename T>
Matrix<T> transp(Matrix<T> A, int n, int m) {
  int m_size = std::max(n, m);
  Matrix<T> AT(m_size, Row<T>(m_size));

  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      AT[i][j] = A[j][i];
    }
  }
  return AT;
}

void print_usage() {
  cout << "Incorrect number of arguments!" << std::endl;
  cout << "usage: conv_Toom-Cook_2D_Huffman "
       << "<kernel_size> <output_size> [no_experiments]"
       << std::endl;
}

/*
 * Prints matrix for debug purposes
 * template function
 *
 * @param m_name name of matrix
 * @param M matrix
 * @param m_rows number of rows
 * @param m_cols number of columns
 */
template<typename T>
void print_matrix_debug(string m_name, Matrix<T> M, int m_rows, int m_cols) {
  if(DEBUG) {
    cout << "  matrix " << m_name << ":" << std::endl;
    print_matrix<T>(M, m_rows, m_cols);
    cout << std::endl;
  }
}


int main(int argc, char* argv[]) {
  if(argc < 3){
    print_usage();
    return 1;
  }

  // size of the kernel
  int n_h = atoi(argv[1]);
  // size of the output
  int n_o = atoi(argv[2]);
  // size of the input
  int n_i = n_h + n_o - 1;
  //int m_h = n_h;
  //int m_o = n_o; //liczba kolumn wyjscia
  //int m_i = m_h + m_o - 1; //liczba kolumn wejscia
  cout << "Experiment parameters:" << std::endl;
  cout << "  Kernel size: " << n_h << std::endl;
  cout << "  Output size: " << n_o << std::endl;
  cout << "  Input size:  " << n_i << std::endl;

  // number of experiments
  int no_exp = 5000;
  if(argc == 4) {
    no_exp = atoi(argv[3]);
    cout << "Performing custom number of experiments: "
         << no_exp
         << std::endl;
  } else {
    cout << "Performing default number of experiments: "
         << no_exp
         << std::endl;
  }

  // Row<double> sd(no_exp);

  // matrix G in single precision
  Matrix<float> G1(n_i, Row<float>(n_h));
  // matrix G^T in single precision
  Matrix<float> G1T(n_h, Row<float>(n_i));
  // matrix G in double precision
  Matrix<double> G1d(n_i, Row<double>(n_h));
  // matrix G^T in double precision
  Matrix<double> G1Td(n_h, Row<double>(n_i));
  // matrix B^T in single precision
  Matrix<float> B1T(n_i, Row<float>(n_i));
  // matrix B in single precision
  Matrix<float> B1(n_i, Row<float>(n_i));
  // matrix B^T in double precision
  Matrix<double> B1Td(n_i, Row<double>(n_i));
  // matrix B in double precision
  Matrix<double> B1d(n_i, Row<double>(n_i));
  // matrix A^T in single precision
  Matrix<float> A1T(n_o, Row<float>(n_i));
  // matrix A in single precision
  Matrix<float> A1(n_i, Row<float>(n_o));
  // matrix A^T in double precision
  Matrix<double> A1Td(n_o, Row<double>(n_i));
  // matrix A in double precision
  Matrix<double> A1d(n_i, Row<double>(n_o));

  // Winograd convolution
  Matrix<float> G(n_i, Row<float>(n_h));
  Matrix<float> GG(n_i, Row<float>(n_h));
  Matrix<float> B(n_i, Row<float>(n_i));
  Matrix<float> BB(n_i, Row<float>(n_i));
  Matrix<float> A(n_o, Row<float>(n_i));
  Matrix<float> AA(n_o, Row<float>(n_i));

  // Winograd convolution with Huffman in mixed prec
  Matrix<float> W(n_h, Row<float>(n_h));
  Matrix<float> Wh(n_h, Row<float>(n_h));
  Matrix<double> Gd(n_i, Row<double>(n_h));
  Matrix<double> GGd(n_i, Row<double>(n_h));
  Matrix<double> Bd(n_i, Row<double>(n_i));
  Matrix<double> BBd(n_i, Row<double>(n_i));
  Matrix<double> Ad(n_o, Row<double>(n_i));
  Matrix<double> AAd(n_o, Row<double>(n_i));

  // Winograd convolution with Huffman
  Matrix<double> Wd(n_h, Row<double>(n_h));
  Matrix<float> Gh(n_i, Row<float>(n_h));
  Matrix<float> GGh(n_i, Row<float>(n_h));
  Matrix<float> Bh(n_i, Row<float>(n_i));
  Matrix<float> BBh(n_i, Row<float>(n_i));
  Matrix<float> Ah(n_o, Row<float>(n_i));
  Matrix<float> AAh(n_o, Row<float>(n_i));

  // kernel matrix in single precision
  Matrix<float> H(n_h, Row<float>(n_h));
  // input matrix in single precision
  Matrix<float> X(n_i, Row<float>(n_i));
  // kernel matrix in double precision
  Matrix<double> Hd(n_h, Row<double>(n_h));
  // input matrix in double precision
  Matrix<double> Xd(n_i, Row<double>(n_i));
  // direct convolution in single precision
  Matrix<float> S(n_o, Row<float>(n_o));
  // direct convolution in double precision
  Matrix<double> D(n_o, Row<double>(n_o));

  // double error; //blad rozwiazania
  // double L2_error = 0;
  double error_W_all = 0;
  double error_Wh_all = 0;
  double error_W_d_all = 0;
  double error_d_all = 0;

  long long num;
  long long den;
  ifstream fB1T("matrix_input/BTMat.txt");
  for(int i = 0; i < n_i; i++) {
    for(int j = 0; j < n_i; j++) {
      fB1T >> num >> den;
      B1T[i][j] = (float)num / (float)den;
      B1Td[i][j] = (double)B1T[i][j];
    }
  }
  fB1T.close();

  ifstream fA1T("matrix_input/ATMat.txt");
  for(int i = 0; i < n_o; i++) {
    for(int j = 0; j < n_i; j++) {
      fA1T >> num >> den;
      A1T[i][j] = (float)num / (float)den;
      A1Td[i][j] = (double)A1T[i][j];
    }
  }
  fA1T.close();

  ifstream fG1("matrix_input/GMat.txt");
  for(int i = 0; i < n_i; i++) {
    for(int j = 0; j < n_h; j++) {
      fG1 >> num >> den;
      G1[i][j] = (float)num / (float)den;
      G1d[i][j] = (double)G1[i][j];
    }
  }
  fG1.close();

  G1T = transp<float>(G1, n_i, n_h);
  G1Td = transp<double>(G1d, n_i, n_h);
  
  B1 = transp<float>(B1T, n_i, n_i);
  B1d = transp<double>(B1Td, n_i, n_i);

  A1 = transp<float>(A1T, n_o, n_i);
  A1d = transp<double>(A1Td, n_o, n_i);

  for(int it = 0; it < no_exp; it++) {
    for(int i = 0; i < n_o; i++) {
      for(int j = 0; j < n_o; j++) {
        S[i][j] = D[i][j] = 0;
      }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1, 1);

    for(int i = 0; i < n_h; i++) {
      for(int j = 0; j < n_h; j++) {
        H[i][j] =  dis(gen);
        Hd[i][j] = (double)H[i][j];
      }
    }

    for(int i = 0; i < n_i; i++) {
      for(int j = 0; j < n_i; j++) {
        X[i][j] = dis(gen);
        Xd[i][j] = (double)X[i][j];
      }
    }

     G = multiplication_matrix<float>(G1, H, n_i, n_h, n_h);
    Gh = multiplication_matrix_left_huffman<float>(G1, H, n_i, n_h, n_h);
    Gd = multiplication_matrix_left_huffman<double>(G1d, Hd, n_i, n_h, n_h);

     GG = multiplication_matrix<float>(G, G1T, n_i, n_i, n_h);
    GGh = multiplication_matrix_right_huffman<float>(Gh, G1T, n_i, n_i, n_h);
    GGd = multiplication_matrix_right_huffman<double>(Gd, G1Td, n_i, n_i, n_h);

     B = multiplication_matrix<float>(B1T, X, n_i, n_i, n_i);
    Bh = multiplication_matrix_left_huffman<float>(B1T, X, n_i, n_i, n_i);
    Bd = multiplication_matrix_left_huffman<double>(B1Td, Xd, n_i, n_i, n_i);

     BB = multiplication_matrix<float>(B, B1, n_i, n_i, n_i);
    BBh = multiplication_matrix_right_huffman<float>(Bh, B1, n_i, n_i, n_i);
    BBd = multiplication_matrix_right_huffman<double>(Bd, B1d, n_i, n_i, n_i);

     W = multiplication_wise<float>(GG, BB, n_i, n_i);
    Wh = multiplication_wise<float>(GGh, BBh, n_i, n_i);
    Wd = multiplication_wise<double>(GGd, BBd, n_i, n_i);

     A = multiplication_matrix<float>(A1T, W, n_o, n_i, n_i);
    Ah = multiplication_matrix_left_huffman<float>(A1T, Wh, n_o, n_i, n_i);
    Ad = multiplication_matrix_left_huffman<double>(A1Td, Wd, n_o, n_i, n_i);

     AA = multiplication_matrix<float>(A, A1, n_o, n_o, n_i);
    AAh = multiplication_matrix_right_huffman<float>(Ah, A1, n_o, n_o, n_i);
    AAd = multiplication_matrix_right_huffman<double>(Ad, A1d, n_o, n_o, n_i);
    print_matrix_debug<float>("AA", AA, n_o, n_o);

    for(int i = 0; i < n_o; i++) {
      for(int j = 0; j < n_o; j++) {
        for(int k = 0; k < n_h; k++) {
          for(int l = 0; l < n_h; l++) {
            S[i][j] = S[i][j] + (H[k][l] * X[k + i][l + j]);
          }
        }
      }
    }
    print_matrix_debug<float>("S", S, n_o, n_o);

    for(int i = 0; i < n_o; i++) {
      for(int j = 0; j < n_o; j++) {
        for(int k = 0; k < n_h; k++) {
          for(int l = 0; l < n_h; l++) {
            D[i][j] = D[i][j] + (Hd[k][l] * Xd[k + i][l + j]);
          }
        }
      }
    }
    print_matrix_debug<double>("D", D, n_o, n_o);

    // Winograd
    double error_W_g = 0;
    // Winograd Huffman
    double error_Wh_g = 0;
    // Winograd mixed Huffman
    double error_d_g = 0;
    // direct
    double error_W_d = 0;

    for(int i = 0; i < n_o; i++) {
      for(int j = 0; j < n_o; j++) {
        error_W_g = error_W_g + abs(AA[i][j] - D[i][j]);
        error_Wh_g = error_Wh_g + abs(AAh[i][j] - D[i][j]);
        error_W_d = error_W_d + abs(AAd[i][j] - D[i][j]);
        error_d_g = error_d_g + abs(S[i][j] - D[i][j]);
      }
    }
    // sd[it] = error_W_d;
    // error_W_g = sqrt(error_W_g);
    // error_Wh_g = sqrt(error_Wh_g);
    // error_W_d = sqrt(error_W_d);
    // error_d_g = sqrt(error_d_g);

    error_W_all = error_W_all + error_W_g;
    error_Wh_all = error_Wh_all + error_Wh_g;
    error_W_d_all = error_W_d_all + error_W_d;
    error_d_all = error_d_all + error_d_g;
  }

  error_W_all = error_W_all / no_exp;
  error_Wh_all = error_Wh_all / no_exp;
  error_W_d_all = error_W_d_all / no_exp;
  error_d_all = error_d_all / no_exp;
  // double Ss = 0;
  // for(int i = 0; i < no_exp; i++)
  //   Ss = Ss + pow(sd[i] - error_W_d_all, 2);
  // Ss = Ss / (no_exp - 1);
  // Ss = sqrt(Ss);
  // double error_W_output = error_W_all / (double)n_o;
  cout << "Experiment results:" << std::endl;
  print_result("Winograd mixed Huffman:\t", error_W_d_all);
  print_result("Winograd:\t\t\t", error_W_all);
  print_result("Winograd huffman:\t\t", error_Wh_all);
  print_result("Direct:\t\t\t", error_d_all);

  // print_result("Double+float", error_W_d_all);
  // print_result("Winograd vs direct:", error_W_d);

 return 0;
}
