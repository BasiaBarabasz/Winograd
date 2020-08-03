#include <iostream>
#include <fstream>
#include <random>
#include "utils.h"

using namespace std;

/*
 * Performs multiplication with Huffman tree of matrix M (n x m) and
 * vector b (m)
 * template function
 *
 * @param M matrix n x m
 * @param n number of rows in matrix
 * @param m number of columns in matrix
 * @param b vector
 * @return vector of size n
 */
template<typename T>
Row<T> multiplication_huffman(Matrix<T> M, int n, int m, Row<T> b) {
  Row<T> result(n);

  // for each row
  for(int i = 0; i < n; i++) {
    Queue<T> Q;
    for(int j = 0; j < m; j++) {
      double e = log2(abs(M[i][j]));
      Tree<T>* t = new Tree<T>(e, M[i][j] * b[j]);
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
    result[i] = p->next->root->root->result();
  }
  return result;
}

/*
 * Performs vector element wise multiplication
 * template function
 *
 * @param a vector of size n
 * @param b vector of size n
 * @param n size of vectors
 * @return vector of size n
 */
template<typename T>
Row<T> multiplication_wise(Row<T> a, Row<T> b, int n) {
  Row<float> c(n);
  Row<T> result(n);

  for(int i = 0; i < n; i++) {
      c[i] = (float)a[i] * (float)b[i];
      result[i] = (T)c[i];
  }
  return result;
}

/*
 * Performs multiplication of matrix M (n x m) and vector (m)
 * template function
 *
 * @param M matrix of size n x m
 * @param n number of rows in matrix M
 * @param m number of columns in matrix M
 * @param b vector of size m
 * @return vector of size n
 */
template<typename T>
Row<T> multiplication(Matrix<T> M, int n, int m, Row<T> b) {
  Row<T> result(n);

  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++) {
      result[i] = result[i] + M[i][j] * b[j];
  }
  return result;
}

// prints program usage
void print_usage() {
  cout << "Incorrect number of arguments!" << std::endl;
  cout << "usage: conv_Toom-Cook_1D "
       << "<kernel_size> <output_size> [no_experiments]"
       << std::endl;
}

/*
 * Prints vector when built for debug
 * template function
 *
 * @param v_name vector name
 * @param v vector
 * @param v_size size of a vector
 */
template<typename T>
void print_vector_debug(string v_name, Row<T> v, int v_size) {
  if(DEBUG) {
    cout << "  vector " << v_name << ":" << std::endl;
    print_vector<T>(v, v_size);
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

  // transformation matrix G
  Matrix<float> G(n_i, Row<float>(n_h));
  // transformation matrix B^T
  Matrix<float> BT(n_i, Row<float>(n_i));
  // transformation matrix A^T
  Matrix<float> AT(n_o, Row<float>(n_i));

  // transformation matrix G in double precision
  Matrix<double> Gd(n_i, Row<double>(n_h));
  // transformation matrix B^T in double precision
  Matrix<double> BTd(n_i, Row<double>(n_i));
  // transformation matrix A^T in double precision
  Matrix<double> ATd(n_o, Row<double>(n_i));

  long long num, den;
  print_debug("Reading matrix G from matrix_input/GMat.txt...");
  ifstream fG("matrix_input/GMat.txt");
  for(int i = 0; i < n_i; i++) {
    for(int j = 0; j < n_h; j++) {
      fG >> num >> den;
      G[i][j] = (float)num / (float)den;
      Gd[i][j] = (double)G[i][j];
    }
  }
  fG.close();

  print_debug("Reading matrix B^T from matrix_input/BTMat.txt...");
  ifstream fBT("matrix_input/BTMat.txt");
  for(int i = 0; i < n_i; i++) {
    for(int j = 0; j < n_i; j++) {
      fBT >> num >> den;
      BT[i][j] = (float)num / (float)den;
      BTd[i][j] = (double)BT[i][j];
    }
  }
  fBT.close();

  print_debug("Reading matrix A^T from matrix_input/ATMat.txt...");
  ifstream fAT("matrix_input/ATMat.txt");
  for(int i = 0; i < n_o; i++) {
    for(int j = 0; j < n_i; j++) {
      fAT >> num >> den;
      AT[i][j] = (float)num / (float)den;
      ATd[i][j] = (double)AT[i][j];
    }
  }
  fAT.close();

  Row<float> g(n_i);
  Row<float> b(n_i);
  Row<float> a(n_o);
  Row<float> w(n_i);
  Row<float> gh(n_i);
  Row<float> bh(n_i);
  Row<float> ah(n_o);
  Row<float> wh(n_i);
  Row<double> gd(n_i);
  Row<double> bd(n_i);
  Row<double> ad(n_o);
  Row<double> wd(n_i);

  // kernel
  Row<float> h(n_h);
  // input
  Row<float> x(n_i);
  // kernel in double precision
  Row<double> hd(n_h);
  // input in double precision
  Row<double> xd(n_i);

  // direct convolution
  Row<float> s(n_o);
  // direct convolution in double precision
  Row<double> d(n_o);

  // winograd huffman in float vs direct in double
  double error_Wh_all = 0;
  // winograd in float vs direct in double
  double error_W_all = 0;
  // direct in float vs direct in double
  double error_d_all = 0;
  // winograd in mixed vs direct in double
  double error_Wd_all = 0;

  for(int it = 0; it < no_exp; it++) {
    print_debug("Running experiment - " + std::to_string(it));
    for(int i = 0; i < n_o; i++)
      s[i] = d[i] = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1, 1);

    for(int i = 0; i < n_h; i++) {
      h[i] =  dis(gen);
      hd[i] = (double)h[i];
    }

    for(int i = 0; i < n_i; i++) {
      x[i] = dis(gen);
      xd[i] = (double)x[i];
    }

/*
    h[0] = hd[0] = 1;
    h[1] = hd[1] = 2;
    h[2] = hd[2] = 1/2;

    x[0] = xd[0] = 1;
    x[1] = xd[1] = 1;
    x[2] = xd[2] = 1;
    x[3] = xd[3] = 1;
    x[4] = xd[4] = 1;
    x[5] = xd[5] = 1;
    x[6] = xd[6] = 1;
    x[7] = xd[7] = 1;
    x[8] = xd[8] = 1;
    x[9] = xd[9] = 1;
    x[10] = xd[10] = 1;
    x[11] = xd[11] = 1;
*/

    gh = multiplication_huffman<float>(G, n_i, n_h, h);
     g = multiplication<float>(G, n_i, n_h, h);
    gd = multiplication_huffman<double>(Gd, n_i, n_h, hd);
    print_vector_debug<float>("g", g, n_i);

    bh = multiplication_huffman<float>(BT, n_i, n_i, x);
     b = multiplication<float>(BT, n_i, n_i, x);
    bd = multiplication_huffman<double>(BTd, n_i, n_i, xd);
    print_vector_debug<float>("b", b, n_i);
  
    wh = multiplication_wise<float>(gh, bh, n_i);
     w = multiplication_wise<float>(g, b, n_i);
    wd = multiplication_wise<double>(gd, bd, n_i);
    print_vector_debug<float>("w", w, n_i);

    ah = multiplication_huffman<float>(AT, n_o, n_i, wh);
     a = multiplication<float>(AT, n_o, n_i, w);
    ad = multiplication_huffman<double>(ATd, n_o, n_i, wd);
    print_vector_debug<float>("a", a, n_o);
  
    for(int i = 0; i < n_o; i++)
      for(int k = 0; k < n_h; k++)
        s[i] = s[i] + (h[k] * x[k + i]);
    print_vector_debug<float>("s", s, n_o);

    for(int i = 0; i < n_o; i++)
      for(int k = 0; k < n_h; k++)
        d[i] = d[i] + (hd[k] * xd[k + i]);
    print_vector_debug<double>("d", d, n_o);

    // Winograd with Huffman
    double error_Wh_g = 0;
    // direct
    double error_d_g = 0;
    // Winograd
    double error_W_g = 0;
    // mixed with Huffman
    double error_Wd = 0;
    // double error_W_d = 0;

    for(int i = 0; i < n_o; i++) {
      error_Wh_g = error_Wh_g + abs(ah[i] - d[i]);
      error_d_g = error_d_g + abs(s[i] - d[i]);
      error_W_g = error_W_g + abs(a[i] - d[i]);
      error_Wd = error_Wd + abs(ad[i] - d[i]);
      // error_W_d = error_W_d + pow((a[i] - s[i]), 2);
    }
    // error_Wh_g = sqrt(error_Wh_g);
    // error_d_g = sqrt(error_d_g);
    // error_W_g = sqrt(error_W_g);
    // error_Wd = sqrt(error_Wd);
    // error_W_d = sqrt(error_W_d);

    error_Wh_all = error_Wh_all + error_Wh_g;
    error_d_all = error_d_all + error_d_g;
    error_W_all = error_W_all + error_W_g;
    error_Wd_all = error_Wd_all + error_Wd;
  }

  error_Wh_all = error_Wh_all / no_exp;
  error_d_all = error_d_all / no_exp;
  error_W_all = error_W_all / no_exp;
  error_Wd_all = error_Wd_all / no_exp;

  // double error_W_output = error_W_all / (double)n_o;
  cout << "Experiment results:" << std::endl;
  print_result("Winograd_huffman:\t", error_Wh_all);
  print_result("Mixed + Huffman:\t", error_Wd_all);
  print_result("Winograd:\t\t", error_W_all);
  print_result("Direct convolution:\t", error_d_all);

  // print_result("Winograd:\t\t", error_W_all);
  // print_result("Double+float:\t\t", error_Wd_all);
  // print_result("Direct:\t\t", error_d_all);
  // print_result("Winograd v direct:\t", error_W_d);

  return 0;
}
