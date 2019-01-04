#include <igl/cotmatrix.h>
#include <igl/readSTL.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <cstring>
#include "stdlib.h"
#include "file_rw.hpp"
#include <fstream>
#include <sstream>
  
int main()
{


  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;

  std::string a = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/test_part1.STL";
  // std::string a = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/test_part1.STL";
  
  igl::readSTL(a, V, F, N);
  std::cout << V<< std::endl;
  std::cout << F<< std::endl << std::endl;
  std::cout << V.rows()<< std::endl;

  // file_rw::file_write_mat("vcpp.csv", V);
  // file_rw::file_write_mat("fcpp.csv", F);
  // file_rw::file_write_mat("ncpp.csv", N);

  return 0;
}



