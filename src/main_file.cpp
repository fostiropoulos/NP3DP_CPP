#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include "nabo/nabo.h"
#include "stdlib.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include "transformation_utilities.hpp"

Eigen::MatrixXd rotate_pts(Eigen::MatrixXd, double, double, double);
Eigen::VectorXd linsp(double, double, double);
Eigen::MatrixXd identify_bottom_layer(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);
Eigen::MatrixXd generate_grid_points(double, double, double, double, double, double, double);

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
	double start_hatch_angle = 3.14159/2;
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

    // take data from file for now
    std::string data_dir = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/"; 
    std::string vert_file = data_dir + "vert.csv";
	std::string face_file = data_dir + "face.csv";
    std::string norm_file = data_dir + "norm.csv";

    Eigen::MatrixXd v;
    Eigen::MatrixXd f;
    Eigen::MatrixXd n;
    
    v = file_rw::file_read_mat(vert_file);
	f = file_rw::file_read_mat(face_file);
    n = file_rw::file_read_mat(norm_file);

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

	int num_of_layers = 1;

	//////////////////////////////////////////////////////
    // identifying the co-ordinates of bottom layer first
	//////////////////////////////////////////////////////
	
	Eigen::MatrixXd fnew = identify_bottom_layer(v,f,n);

	////////////////////////////
    // generating the tool path
	////////////////////////////

	for (int layer=1;layer<=num_of_layers;++layer)
	{
		double hatch_angle = start_hatch_angle + (layer-1)*hatch_angle_change;

		Eigen::MatrixXd	pts = generate_grid_points(pathgap_x,pathgap_y,xmin,ymin,xmax,ymax,hatch_angle);
		std::cout << pts << std::endl;
		

	}






    return 0;
}

Eigen::MatrixXd rotate_pts(Eigen::MatrixXd pts, double theta, double x_avg, double y_avg)
{
	// Function for applying rotation to points about their centroid
	// INPUT : xyz points, rotation angle, x-avg and y-avg
	// OUTPUT : rotated points about their centroid

	// translate points so that origin matches with the centroid
	pts.block(0,0,pts.rows(),1) = pts.block(0,0,pts.rows(),1).array() - x_avg;
	pts.block(0,1,pts.rows(),1) = pts.block(0,1,pts.rows(),1).array() - y_avg;

	// rotate points about centroid
	Eigen::Matrix3d r = rtf::rot_z(theta);
	Eigen::Vector3d t;
	t << 0,0,0;
	Eigen::Matrix4d T = rtf::hom_T(t, r);
	Eigen::MatrixXd pts_before_T(pts.rows(),pts.cols()+1);
	pts_before_T << pts,Eigen::MatrixXd::Constant(pts.rows(),1,0);
	Eigen::MatrixXd	pts_T = rtf::apply_transformation(pts_before_T, T); 

	// translate points back to the original position
	pts_T.block(0,0,pts_T.rows(),1) = pts_T.block(0,0,pts_T.rows(),1).array() + x_avg;
	pts_T.block(0,1,pts_T.rows(),1) = pts_T.block(0,1,pts_T.rows(),1).array() + y_avg;

	return pts_T.block(0,0,pts_T.rows(),2);
}

Eigen::MatrixXd generate_grid_points(double pathgap_x, double pathgap_y, double xmin, double ymin, double xmax, double ymax, double hatch_angle)
{
	// Function to generate the uniform mesh grid of points along the x-y plane
	// INPUT = gap between the adjacent points and maximum value in x and y direction
	// OUTPUT = All points consisting the uniform grid

	Eigen::VectorXd j = linsp(floor(ymin),ceil(ymax),pathgap_y);
	Eigen::VectorXd i_val = linsp(floor(xmin),ceil(xmax),pathgap_x);
	Eigen::MatrixXd	pts = Eigen::MatrixXd::Constant(j.rows()*i_val.rows(),2,0);
	long int st_pt = 0;
	for (long int i=0;i<i_val.rows();++i)
	{
		pts.block(st_pt,0,j.rows(),1) = i_val(i)*Eigen::MatrixXd::Constant(j.rows(),1,1);
		pts.block(st_pt,1,j.rows(),1) = j.block(0,0,j.rows(),j.cols());
		st_pt = st_pt + j.rows();
	}

	// apply rotation to points
	double x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
	double y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
	Eigen::MatrixXd rotated_pts = rotate_pts(pts,hatch_angle,x_avg,y_avg);
	return rotated_pts;	
}

Eigen::MatrixXd identify_bottom_layer(Eigen::MatrixXd v, Eigen::MatrixXd f, Eigen::MatrixXd n)
{
	// This function identifies the vertices of bottom layer first based on the normal direction
	// INPUT = vertices (3n-by-3 matrix), faces (n-by-3 matrix) and normals (n-by-3 matrix) generated by the STL file.
	// OUTPUT = vertices of the of the bottom surface (m-by-3 matrix)

	Eigen::VectorXd indexv = linsp(1,v.rows(),1);		// indexing vertices
	Eigen::VectorXd indexn = linsp(1,n.rows(),1);		// indexing normals...same as indexing faces
	Eigen::MatrixXd nnew(n.rows(),7);
	nnew << n,f,indexn;
	Eigen::MatrixXd vnew(v.rows(),4);
	vnew << v,indexv; 

	// storing index of only those normals whose z direction is negative,
	// i.e. normal pointing downwards... i.e normals from bottom layer
	Eigen::MatrixXd	store = Eigen::MatrixXd::Constant(nnew.rows(),nnew.cols(),0);
	long int idx_n = 0;
	for (long int i=0;i<nnew.rows();++i)
	{
		if (nnew(i,2)<0)
		{
			store.block(idx_n,0,1,store.cols()) = nnew.block(i,0,1,nnew.cols());
			++idx_n;
		}
	}

	// using index of those normals to get the faces and hence the points associated with those faces
	Eigen::MatrixXd fnew = store.block(0,3,idx_n,3);
	return fnew;
}


Eigen::VectorXd linsp(double strt, double end, double stp)
{
	int sz;
	if (strt<=end && stp>0 || strt>=end && stp<0)
	{
		sz = int((end-strt)/stp)+1;
	}
	else
	{
		if (strt>=end)
		{
			std::cerr << "start value is greater than the end value for incement!" << std::endl;
			std::terminate();	
		}
		else
		{
			std::cerr << "start value is less than the end value for decrement!" << std::endl;
			std::terminate();	
		}
	}
	return Eigen::VectorXd::LinSpaced(sz,strt,strt+stp*(sz-1));
}
