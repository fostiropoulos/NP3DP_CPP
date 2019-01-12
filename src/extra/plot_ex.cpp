#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <igl/readSTL.h>
#include <igl/point_in_poly.h>
#include "stdlib.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include "transformation_utilities.hpp"
#include "NPAM_utilities.hpp"
#include <vector>
#include <cmath>
#include <cstdlib>
#define GNUPLOT_ENABLE_PTY
#include <boost/tuple/tuple.hpp>
# include "gnuplot-iostream.h"

int main() {

	Gnuplot gp;
	
	Eigen::MatrixXd m(4,3);

	m << 1,2,3,4,5,6,7,8,9,10,11,12;
	std::vector<std::vector<double> > pts = ut::mat_to_vec(m);

	std::vector<std::vector<double> > pts2 = file_rw::file_read_vec("/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/complete_tool_path.csv");

	int a = 200; 
	// gp << "set xrange [0:200]\n";
	// gp << "set yrange [0:200]\n";
	// gp << "set zrange [0:200]\n";
	gp << "set xlabel 'X-axis'\n";
	gp << "set ylabel 'Y-axis'\n";
	gp << "set zlabel 'Z-axis'\n";
	
	// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
	gp << "splot '-' with lines title 'abc' \n";
	// gp << "show view \n";
	// gp << "using 1:1:1 \n";
	
	gp.send1d(pts2);
	double mx=0.5, my=0.5;
	int mb=1;
	
	while(mb != 3 && mb >= 0) 
	{
		gp.getMouse(mx, my, mb, "Left click to aim arrows, right click to exit.");
		printf("You pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
		if(mb < 0) 
			{
				printf("The gnuplot window was closed.\n");
			}
	}
	
}