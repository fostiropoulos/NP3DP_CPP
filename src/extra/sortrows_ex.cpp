#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <igl/readSTL.h>
#include <igl/point_in_poly.h>
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


int main()
{
	std::vector<std::vector<double> > vec {{1,2},{2,1},{3,3},{1,1},{3,1},{1,3}};

	ut::disp_vec(vec);
	
	std::vector<std::vector<double> > vec2;
	std::vector<std::vector<double> > vec3;
	
	vec2 = SortRows(vec,1);
	// std::sort(vec.begin(), vec.end(), sortcol);
        
	std::cout << std::endl << std::endl;
	ut::disp_vec(vec2);
	// std::cout << vec2[0] << std::endl;

	vec3 = SortRows(vec,0);
	// std::sort(vec.begin(), vec.end(), sortcol);
        
	std::cout << std::endl << std::endl;
	ut::disp_vec(vec3);
}

