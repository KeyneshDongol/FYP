//
//  main.cpp
//  Polynoimial Operator
//
//  Created by Keynesh on 15/2/23.
//

#include <iostream>
#include "Polynomial.hpp"
//#include <typeinfo>

int main(int argc, const char * argv[]) {
    // insert code here...
//    Polynomial<3> somePoly3(Eigen::Matrix<double, 35, 1>::Ones());
    
    //TODO: Figure out how to create polynomials and do basic operations on them
    
    double a=1;
    double b=1;
    double c=1;
    
    
    Eigen::Matrix<double, 1, 3> Row = Eigen::Matrix<double, 1,3>::Ones();
    Eigen::Matrix<double, 3, 1> Col = Eigen::Matrix<double, 3,1>::Ones();
    Eigen::Matrix<double, 3, 1> Col3 = Col.array()*Col.array()*Col.array();
    
    auto newMatrix = (Row*(Col3*a*b*c));
    auto newerMatrix = (Row*Col);
    auto newNewerMatrix = (Row*Col3*a);
    auto extremelyNewMatrix = newMatrix + newerMatrix + newNewerMatrix;
    
//    std::cout << extremelyNewMatrix << "\n";
//    std::cout << "Type of x : " << typeid(extremelyNewMatrix.value()).name()<< "\n";
//    std::cout << "Type of x : " << typeid(newNewerMatrix).name()<< ">>>\n";
//    std::cout << "Type of x : " << typeid(newerMatrix).name()<< "\n";
    std::cout << Row << "\n";
//    
    
    std::cout << "Hello, World!\n";
    return 0;
}




//#define EIGEN_USE_BLAS
//#define EIGEN_USE_MKL_ALL

//#include <iostream>
//#include <chrono>
//#include <Eigen/Dense>
//using namespace Eigen;
//using namespace std::chrono;
//int main()
//{
//
//  int n_a_rows = 10000;
//  int n_a_cols = 10000;
//  int n_b_rows = n_a_cols;
//  int n_b_cols = 1000;
//
//  MatrixXd a(n_a_rows, n_a_cols);
//
//  for (int i = 0; i < n_a_rows; ++ i)
//      for (int j = 0; j < n_a_cols; ++ j)
//        a (i, j) = n_a_cols * i + j;
//
//  MatrixXd b (n_b_rows, n_b_cols);
//  for (int i = 0; i < n_b_rows; ++ i)
//      for (int j = 0; j < n_b_cols; ++ j)
//        b (i, j) = n_b_cols * i + j;
//
//  MatrixXd d (n_a_rows, n_b_cols);
//
//  using wall_clock_t = std::chrono::high_resolution_clock;
//  auto const start = wall_clock_t::now();
//  clock_t begin = clock();
//
//  d = a * b;
//
//  clock_t end = clock();
//  auto const wall = std::chrono::duration<double>(wall_clock_t::now() - start);
//
//  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//  std::cout << "CPU time : " << elapsed_secs << std::endl;
//  std::cout << "Wall time : " << wall.count() << std::endl;
//  std::cout << "Speed up : " << elapsed_secs/wall.count() << std::endl;
//
//}
