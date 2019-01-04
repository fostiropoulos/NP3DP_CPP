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

int main()
{
    ///////////////////////////////////////////////////////
    ////////////////////// INPUTS /////////////////////////
	///////////////////////////////////////////////////////

	// Gap between 2 hatching lines -
	double pathgap_x = 2.0;
	double pathgap_y = 2.0;

	//Gap between 2 layers -
	double pathgap_z = 2.0;
	
	// grid aadition (NOTE: higher number for higher aspect ratio of part in xy plane)
	double grid_addition = 70.0;

	// hatching angle properties
	double start_hatch_angle = 3.14159/4;
	double hatch_angle_change = 0;	// consecutive layers will have this much change in hatching angle (degrees)
	
	// Path to generate -
	// 1 - Boundary
	// 2 - Hatching
	int Path_Number = 2;    // put number from above
	int FlipTravel = 0;     // 1 for yes, 0 for No
	int space = 2;          // number of points to skip for smoother path

    ///////////////////////////
    // Bottom layer generation
	///////////////////////////

	Eigen::MatrixXd v;
	Eigen::MatrixXd f;
	Eigen::MatrixXd n;
	std::string stl_file_name = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/test_part4.STL";
	igl::readSTL(stl_file_name, v, f, n);
	f = f.array()+1;	//NOTE : adding 1 make make start index as 1... this is compatible for rest of the code 
    Eigen::MatrixXd fillpts;
    
	//////////////////////////////////////
    // xmax and ymax to get the grid size
	//////////////////////////////////////

    double xmax = v.block(0,0,v.rows(),1).maxCoeff() + grid_addition;
	double ymax = v.block(0,1,v.rows(),1).maxCoeff() + grid_addition;
	double xmin = v.block(0,0,v.rows(),1).minCoeff() - grid_addition;
    double ymin = v.block(0,1,v.rows(),1).minCoeff() - grid_addition;

	//////////////////////////////////////////////////////////////////////////////
    // finding thickness of shell part and number of layers required for printing
	//////////////////////////////////////////////////////////////////////////////

	int num_of_layers = NPAM::number_of_layers(v,f,n,pathgap_z);

	//////////////////////////////////////////////////////
    // identifying the co-ordinates of bottom layer first
	//////////////////////////////////////////////////////
	
	Eigen::MatrixXd fnew = NPAM::identify_bottom_layer(v,f,n);

	////////////////////////////
    // generating the tool path
	////////////////////////////

	Eigen::MatrixXd	pts;
	double hatch_angle;
	double x_avg;
	double y_avg;
	for (int layer=1;layer<=num_of_layers;++layer)
	{
		hatch_angle = start_hatch_angle + (layer-1)*hatch_angle_change;

		pts = NPAM::generate_grid_points(pathgap_x,pathgap_y,xmin,ymin,xmax,ymax,hatch_angle);
		// apply rotation to points
		x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
		y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
		Eigen::MatrixXd rotated_pts = NPAM::rotate_pts(pts,hatch_angle,x_avg,y_avg);
		
		fillpts = NPAM::project_grid_points(fnew, v, rotated_pts, hatch_angle, x_avg, y_avg);
	}

	////////////////////////////
    // infill path
	////////////////////////////

	Eigen::MatrixXd	tool_path = NPAM::Infill_Path(fillpts, FlipTravel, space, hatch_angle, x_avg, y_avg);
	file_rw::file_write("tool_path.csv",tool_path);
	std::cout << "Tool Path Generated!" << std::endl;
    return 0;
}
