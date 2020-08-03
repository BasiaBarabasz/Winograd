#include <iostream>
#include <iomanip>
#include "utils.h"

void print_debug(string msg) {
  if(DEBUG) {
    cout << msg << std::endl;
  }
}

void print_result(string label, double error) {
  cout << "  "
       << label
       << std::setprecision(12)
       << std::fixed
       << error
       << std::endl;
}
