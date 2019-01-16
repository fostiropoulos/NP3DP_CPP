#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <igl/readSTL.h>
#include "stdlib.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include "transformation_utilities.hpp"
#include "NPAM_utilities.hpp"
#include <cstdlib>
#include <boost/tuple/tuple.hpp>
#include "plot_utilities.hpp"
#include <thread>

// A dummy function 
void fun1() 
{ 
    std::cout << "thread 1 activated" << std::endl; 
    while(true)
    {
    	int a = 1;
    }
} 

void fun2() 
{ 
    std::cout << "thread 2 activated" << std::endl; 
    while(true)
    {
    	int b = 1;
    }
} 

int main()
{
	std::cout << "start" << std::endl;
	std::thread thd1(fun1);
	std::thread thd2(fun2);
	thd1.join();
	thd2.join();
	// fun1();
	std::cout << "end" << std::endl;

}	