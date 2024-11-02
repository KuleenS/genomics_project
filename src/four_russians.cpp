// cppttest.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <string>
#include <cstdio>
#include <sstream>
#include "four_russians/FourRussians.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<input_file> <t_value>" << std::endl;

        return 1;
    }

  std::ifstream input(argv[1]);

  if (!input.is_open()){
      std::cerr << "File " << argv[1] <<  " does not exist" << std::endl;
      return 1;
  }

  std::string x,y; 

  std::getline(input, x);
  std::getline(input, y);

  input.close();

  std::istringstream t_value_ss(argv[2]);

  int t_value;

  if (!(t_value_ss >> t_value)) {
      std::cerr << "Invalid number: " << argv[3] << '\n';
      return 1;
  } else if (!t_value_ss.eof()) {
      std::cerr << "Trailing characters after number: " << argv[3] << '\n';
      return 1;
  }

  FourRussians fourRussians(x, y, t_value);

  fourRussians.parallelPreProcessing();

  int** distanceMatrix = fourRussians.calculateEditMatrix();

  int minDistance = fourRussians.calculateMinDistance(distanceMatrix);

  std::cout << minDistance << std::endl;

  return 0;
}
