#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <igl/readSTL.h>
#include "nabo/nabo.h"
#include "stdlib.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include "transformation_utilities.hpp"

template <typename t>
    std::vector<std::vector<t> > SortRows(std::vector<std::vector<t> > input, int col)
    {
    	for (long int i=0;i<input.size();++i)
    	{
    		t temp = input[i][0];
    		input[i][0] = input[i][col];
    		input[i][col] = temp; 
    	}
        std::sort(input.begin(), input.end());
        for (long int i=0;i<input.size();++i)
    	{
    		t temp = input[i][0];
    		input[i][0] = input[i][col];
    		input[i][col] = temp; 
    	}
        return input;   
    }

Eigen::VectorXi InPoly(Eigen::MatrixXd, Eigen::MatrixXd);
std::vector<int> InPoly2(Eigen::MatrixXd, Eigen::MatrixXd);
std::vector<int> find_idx(std::vector<int>);


int main()
{

    Eigen::MatrixXd v(3,2);
    Eigen::MatrixXd pts(5,2);
    v << 0.004,-12.001,4.82,7.42,-7.14,7.8;
    pts << 1,2,3,4,5,6,7,8,0.9,1;

    // Eigen::MatrixXi in = InPoly(v,pts);
    std::vector<int> in = InPoly2(v,pts);


    for (int i=0;i<in.size();++i)
    {
        std::cout << in[i] << std::endl;    
    }

    std::vector<int> idx = find_idx(in);    
    std::cout << std::endl << std::endl << std::endl;

    for (int i=0;i<idx.size();++i)
    {
        std::cout << idx[i] << std::endl;    
    }     

    std::cout << std::endl << std::endl << std::endl;
    
    std::cout << pts << std::endl;
    std::cout << std::endl << std::endl << std::endl;

    for (int i=0;i<idx.size();++i)
    {
        std::cout << pts.row(idx[i]) << std::endl;
    }
    // std::cout << pts(idx[0],0) << std::endl;


    return 0;
}

Eigen::VectorXi InPoly(Eigen::MatrixXd v, Eigen::MatrixXd pts)
{
    Eigen::VectorXi inp(pts.rows(),1);
    double A = 0.0;
    double A1 = 0.0;
    double A2 = 0.0;
    double A3 = 0.0;
    double tol = 1e-10;
    for (int i=0;i<pts.rows();++i)
    {
        A = fabs(0.5*(v(0,0)*(v(1,1)-v(2,1))-v(0,1)*(v(1,0)-v(2,0))+v(1,0)*v(2,1)-v(2,0)*v(1,1)));
        A1 = fabs((0.5*(pts(i,0)*(v(1,1)-v(2,1))-pts(i,1)*(v(1,0)-v(2,0))+v(1,0)*v(2,1)-v(2,0)*v(1,1))));
        A2 = fabs(0.5*(v(0,0)*(pts(i,1)-v(2,1))-v(0,1)*(pts(i,0)-v(2,0))+pts(i,0)*v(2,1)-v(2,0)*pts(i,1)));
        A3 = fabs(0.5*(v(0,0)*(v(1,1)-pts(i,1))-v(0,1)*(v(1,0)-pts(i,0))+v(1,0)*pts(i,1)-pts(i,0)*v(1,1)));
        if (fabs(A-A1-A2-A3)<tol)
        {
            inp(i,0) = 1;
        }
        else
        {
            inp(i,0) = 0;
        }
    }
    return inp;
}

std::vector<int> InPoly2(Eigen::MatrixXd v, Eigen::MatrixXd pts)
{
    std::vector<int> inp;
    double A = 0.0;
    double A1 = 0.0;
    double A2 = 0.0;
    double A3 = 0.0;
    double tol = 1e-10;
    for (int i=0;i<pts.rows();++i)
    {
        A = fabs(0.5*(v(0,0)*(v(1,1)-v(2,1))-v(0,1)*(v(1,0)-v(2,0))+v(1,0)*v(2,1)-v(2,0)*v(1,1)));
        A1 = fabs((0.5*(pts(i,0)*(v(1,1)-v(2,1))-pts(i,1)*(v(1,0)-v(2,0))+v(1,0)*v(2,1)-v(2,0)*v(1,1))));
        A2 = fabs(0.5*(v(0,0)*(pts(i,1)-v(2,1))-v(0,1)*(pts(i,0)-v(2,0))+pts(i,0)*v(2,1)-v(2,0)*pts(i,1)));
        A3 = fabs(0.5*(v(0,0)*(v(1,1)-pts(i,1))-v(0,1)*(v(1,0)-pts(i,0))+v(1,0)*pts(i,1)-pts(i,0)*v(1,1)));
        if (fabs(A-A1-A2-A3)<tol)
        {
            inp.push_back(1);
        }
        else
        {
            inp.push_back(0);
        }
    }
    return inp;
}

std::vector<int> find_idx(std::vector<int> vec)
{
    std::vector<int> idx;
    for (int i=0;i<vec.size();++i)
    {
        if (vec[i]!=0)
        {
            idx.push_back(i);   
        }
    }
    return idx;
}

