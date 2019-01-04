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

Eigen::MatrixXd rotate_pts(Eigen::MatrixXd, double, double, double);
Eigen::VectorXd linsp(double, double, double);
int number_of_layers(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, double);
Eigen::MatrixXd identify_bottom_layer(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);
Eigen::MatrixXd generate_grid_points(double, double, double, double, double, double, double);
std::vector<int> InPoly(Eigen::MatrixXd, Eigen::MatrixXd);
Eigen::MatrixXd Infill_Path(Eigen::MatrixXd, bool, double, double, double, double);
Eigen::MatrixXd	project_grid_points(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, double, double, double);
std::vector<int> find_idx(std::vector<int>);
std::vector<int> ismember(Eigen::MatrixXd, Eigen::MatrixXd);
double median(std::vector<double>);

template <typename t> std::vector<std::vector<t> > SortRows(std::vector<std::vector<t> > input, int col)
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
    ///////////////////////////////////////////////////////
    ////////////////////// INPUTS /////////////////////////
	///////////////////////////////////////////////////////

	// Gap between 2 hatching lines -
	double pathgap_x = 0.5;
	double pathgap_y = 0.5;

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

	// STLREAD- need to work this function ...generated data is not good
	Eigen::MatrixXd v;
	Eigen::MatrixXd f;
	Eigen::MatrixXd n;
	std::string stl_file_name = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/test_part1.STL";
	igl::readSTL(stl_file_name, v, f, n);
	f = f.array()+1;	//NOTE : adding 1 make make start index as 1... this is compatible for rest of the code 
    
    // take data from file for now
 	// std::string data_dir = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/"; 
	// std::string vert_file = data_dir + "vert.csv";
	// std::string face_file = data_dir + "face.csv";
	// std::string norm_file = data_dir + "norm.csv";
	// std::string fillpts_file = data_dir + "fillpts.csv";			////// DELETE THIS //////

    // Eigen::MatrixXd v;
    // Eigen::MatrixXd f;
    // Eigen::MatrixXd n;
    Eigen::MatrixXd fillpts;					/////// DELETE THIS //////
    
    // v = file_rw::file_read_mat(vert_file);
	// f = file_rw::file_read_mat(face_file);
    // n = file_rw::file_read_mat(norm_file);


    // fillpts = file_rw::file_read_mat(fillpts_file);    /////// DELETE THIS //////

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

	// int num_of_layers = 1;
	int num_of_layers = number_of_layers(v,f,n,pathgap_z);

	//////////////////////////////////////////////////////
    // identifying the co-ordinates of bottom layer first
	//////////////////////////////////////////////////////
	
	Eigen::MatrixXd fnew = identify_bottom_layer(v,f,n);

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

		pts = generate_grid_points(pathgap_x,pathgap_y,xmin,ymin,xmax,ymax,hatch_angle);
		// apply rotation to points
		x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
		y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
		Eigen::MatrixXd rotated_pts = rotate_pts(pts,hatch_angle,x_avg,y_avg);
		
		fillpts = project_grid_points(fnew, v, rotated_pts, hatch_angle, x_avg, y_avg);
	}

	////////////////////////////
    // infill path
	////////////////////////////

	Eigen::MatrixXd	tool_path = Infill_Path(fillpts, FlipTravel, space, hatch_angle, x_avg, y_avg);

	// std::cout << tool_path << std::endl;
	file_rw::file_write_mat("test_tool_path.csv",tool_path);


    return 0;
}


std::vector<int> ismember(Eigen::MatrixXd mat, Eigen::MatrixXd row_vec)
{
	////////MODIFY THIS>>>>THIS IS ONLY PRIMARY >>>CHANGES NEEDED//////
	std::vector<int> v;
	for (long int i=0;i<mat.rows();++i)
	{
		if (mat.row(i)==row_vec)
		{
			v.push_back(1);
		}
		else
		{
			v.push_back(0);
		}
	}
	return v;
}

double median(std::vector<double> vec)
{
	std::sort(vec.begin(),vec.end());
	if (vec.size()%2==1)
	{
		return vec[vec.size()/2];
	}
	else
	{
		return (vec[(vec.size()/2)-1]+vec[vec.size()/2])/2;	
	}
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

int number_of_layers(Eigen::MatrixXd v, Eigen::MatrixXd f, Eigen::MatrixXd n, double pathgap_z)
{
	Eigen::MatrixXd n_new(n.rows(),4);
    n_new.block(0,0,n_new.rows(),3) = (((n.array()*100.0).cast<int>()).cast<double>().array()/100.0);

    n_new.block(0,3,n_new.rows(),1) = linsp(0,n_new.rows()-1,1); 
    // std::cout << n_new << std::endl;

    // sorting normals for top and bottom face
    Eigen::MatrixXd nn_i = Eigen::MatrixXd::Constant(n_new.rows(),n_new.cols(),0);
    Eigen::MatrixXd np_i = Eigen::MatrixXd::Constant(n_new.rows(),n_new.cols(),0);
    long int nn_idx = 0;
    long int np_idx = 0;
    for (long int i=0;i<n_new.rows();++i)
    {
    	if (n_new(i,2)<0)
    	{
    		nn_i.row(nn_idx) = n_new.row(i);
    		++nn_idx;
    	}
    	else if (n_new(i,2)>0)
    	{
    		np_i.row(np_idx) = n_new.row(i);
    		++np_idx;
    	}
    }
    Eigen::MatrixXd nn = nn_i.block(0,0,nn_idx,nn_i.cols());
	Eigen::MatrixXd np = np_i.block(0,0,np_idx,np_i.cols());

    // getting normals which are opposite ot each other and forming pairs of parallel faces 
	Eigen::MatrixXd p11(1,3);
	Eigen::MatrixXd p12(1,3);
	Eigen::MatrixXd p13(1,3);
	Eigen::MatrixXd p21(1,3);
	double a = 0;
	double b = 0;
	double c = 0;
	double d = 0;
	
	std::vector<double> t;

	for (long int i=0;i<nn.rows();++i)
	{
		std::vector<int> idx = ismember(np.block(0,0,np.rows(),3),-nn.block(i,0,1,3));
		std::vector<int> store_idx = find_idx(idx);
		if (store_idx.size()!=0)
		{
			// calculating the gap between parallel faces by plane to plane distance formaula
			long int pair[2];
			pair[0] = np(store_idx[0],3);
			pair[1] = nn(i,3);
			p11 = v.row(f(pair[0],0)-1);
			p12 = v.row(f(pair[0],1)-1);
			p13 = v.row(f(pair[0],2)-1);
			p21 = v.row(f(pair[1],0)-1);

			a = ((p12(0,1)-p11(0,1))*(p13(0,2)-p11(0,2)))-((p13(0,1)-p11(0,1))*(p12(0,2)-p11(0,2)));
        	b = ((p12(0,2)-p11(0,2))*(p13(0,0)-p11(0,0)))-((p13(0,2)-p11(0,2))*(p12(0,0)-p11(0,0)));
        	c = ((p12(0,0)-p11(0,0))*(p13(0,1)-p11(0,1)))-((p13(0,0)-p11(0,0))*(p12(0,1)-p11(0,1)));
        	d = -(a*p11(0,0))-(b*p11(0,1))-(c*p11(0,2));

        	t.push_back(fabs((a*p21(0,0)+b*p21(0,1)+c*p21(0,2)+d)/sqrt((a*a)+(b*b)+(c*c))));
		}

	}

	double t_f = median(t);
	int num_of_layers = round(t_f/pathgap_z);
	if (num_of_layers==0)
	{
		num_of_layers=1;
	}

	return num_of_layers;
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

	// // apply rotation to points
	// double x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
	// double y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
	// Eigen::MatrixXd rotated_pts = rotate_pts(pts,hatch_angle,x_avg,y_avg);
	// return rotated_pts;

	return pts;	
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

std::vector<int> InPoly(Eigen::MatrixXd v, Eigen::MatrixXd pts)
{
    std::vector<int> inp;
    double A = 0.0;
    double A1 = 0.0;
    double A2 = 0.0;
    double A3 = 0.0;
    double tol = 1e-12;
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

Eigen::MatrixXd Infill_Path(Eigen::MatrixXd fillpts, bool FlipTravel, double space, double hatch_angle, double x_avg, double y_avg)
{
	// function for creating path from projected points for tcp travel along x direction
	// INPUT: projected points on the surface (equally spaced along x and y axes)
	// OUTPUT: points arranged along 0 degree path with their Rx and Ry value

	int dir1 = 1;
	int dir2 = 0;
	Eigen::MatrixXd allpts = fillpts;
	std::vector<std::vector<double> > allpts_vec;
	allpts_vec = ut::mat_to_vec(allpts);
	allpts_vec = SortRows(allpts_vec, dir1);
	allpts = ut::vec_to_mat(allpts_vec);
	int flip  = 0;

	Eigen::MatrixXd	storeset = Eigen::MatrixXd::Constant(allpts.rows(),allpts.cols(),0);
	Eigen::MatrixXd	storesort = Eigen::MatrixXd::Constant(allpts.rows(),allpts.cols(),0);
	long int storesort_strt = 0;
	long int storeset_idx = 0;

	for (long int i=0;i<allpts.rows()-1;++i)
	{
		if (std::abs(allpts(i,dir1)-allpts(i+1,dir1))<0.00001)
		{
			storeset.block(storeset_idx,0,1,allpts.cols()) = allpts.block(i,0,1,allpts.cols());
			++storeset_idx;
		}	
		else
		{
			storeset.block(storeset_idx,0,1,allpts.cols()) = allpts.block(i,0,1,allpts.cols());
			Eigen::MatrixXd storeset_nz(storeset_idx+1,storeset.cols());
			storeset_nz = storeset.block(0,0,storeset_idx+1,storeset.cols());
			std::vector<std::vector<double> > storeset_nz_vec;
			storeset_nz_vec = ut::mat_to_vec(storeset_nz);
			storeset_nz_vec = SortRows(storeset_nz_vec,dir2);
			storeset_nz = ut::vec_to_mat(storeset_nz_vec);
			if (double(flip)/2==round(double(flip)/2))
			{
				Eigen::MatrixXd storeset_nz_rev(storeset_nz.rows(),storeset_nz.cols());	
				storeset_nz_rev = storeset_nz.colwise().reverse();
				storeset_nz = storeset_nz_rev;
			}
			flip++;
			storesort.block(storesort_strt,0,storeset_nz.rows(),storeset_nz.cols()) = storeset_nz;
			storesort_strt = storesort_strt + storeset_idx+1;
			storeset_idx = 0;
		}		
	}

	storeset.block(storeset_idx,0,1,allpts.cols()) = allpts.block(allpts.rows()-1,0,1,allpts.cols());

	Eigen::MatrixXd storeset_nz(storeset_idx+1,storeset.cols());
	storeset_nz = storeset.block(0,0,storeset_idx+1,storeset.cols());
	std::vector<std::vector<double> > storeset_nz_vec;
	storeset_nz_vec = ut::mat_to_vec(storeset_nz);
	storeset_nz_vec = SortRows(storeset_nz_vec,dir2);
	storeset_nz = ut::vec_to_mat(storeset_nz_vec);

	// this is to get the direction of last line travel
	if (storesort(storesort_strt-1,dir2)==storeset_nz(0,dir2))
	{
		storesort.block(storesort_strt,0,storeset_nz.rows(),storeset_nz.cols()) = storeset_nz;	
	}
	else
	{
		Eigen::MatrixXd storeset_nz_rev(storeset_nz.rows(),storeset_nz.cols());
		storeset_nz_rev = storeset_nz.colwise().reverse();
		storeset_nz = storeset_nz_rev;	
		storesort.block(storesort_strt,0,storeset_nz.rows(),storeset_nz.cols()) = storeset_nz;
	}

	Eigen::VectorXd	Rx = -(storesort.block(0,4,storesort.rows(),1).array()/storesort.block(0,5,storesort.rows(),1).array())*(180/3.14159);
	Eigen::VectorXd	Ry = (storesort.block(0,3,storesort.rows(),1).array()/storesort.block(0,5,storesort.rows(),1).array())*(180/3.14159);


	Eigen::MatrixXd storesort0x1(storesort.rows(),5);
	storesort0x1.block(0,0,storesort.rows(),3) = storesort.block(0,0,storesort.rows(),3);
	storesort0x1.block(0,3,storesort.rows(),1) = Rx;
	storesort0x1.block(0,4,storesort.rows(),1) = Ry;

	// storing every n'th point to smoothen out the path 
	int count = 0;
	Eigen::MatrixXd store_spaced_pt(storesort0x1.rows(),storesort0x1.cols()); 
	store_spaced_pt.block(0,0,1,storesort0x1.cols()) = storesort0x1.block(0,0,1,storesort0x1.cols());
	long int store_spaced_pt_idx = 1;
	int flagg = 0;
	for (long int i=1;i<storesort0x1.rows()-1;++i)
	{
		if (flagg==1)
		{
			flagg=0;
			continue;
		}
		if (storesort0x1(i,1)==storesort0x1(i+1,1))
		{
			++count;
			if (double(count)/space == round(count/space))
			{
				store_spaced_pt.block(store_spaced_pt_idx,0,1,storesort0x1.cols()) = storesort0x1.block(i,0,1,storesort0x1.cols());
				++store_spaced_pt_idx;
				count = 0;
			}
		}
		else
		{
			store_spaced_pt.block(store_spaced_pt_idx,0,1,storesort0x1.cols()) = storesort0x1.block(i,0,1,storesort0x1.cols());
			++store_spaced_pt_idx;
			store_spaced_pt.block(store_spaced_pt_idx,0,1,storesort0x1.cols()) = storesort0x1.block(i+1,0,1,storesort0x1.cols());
			++store_spaced_pt_idx;
			flagg = 1;	
		}
	}

	// adding the last point
	store_spaced_pt.block(store_spaced_pt_idx,0,1,storesort0x1.cols()) = storesort0x1.block(storesort0x1.rows()-1,0,1,storesort0x1.cols());
	Eigen::MatrixXd storesort0x1tp = store_spaced_pt.block(0,0,store_spaced_pt_idx+1,store_spaced_pt.cols());

	// flip the direction of travel
	if (FlipTravel==1)
	{
		Eigen::MatrixXd storesort0x1tp_rev(storesort0x1tp.rows(),storesort0x1tp.cols());
		storesort0x1tp_rev = storesort0x1tp.colwise().reverse();
		storesort0x1tp = storesort0x1tp_rev;
	}

	//apply hatching angle
	
	Eigen::MatrixXd storesort0x1tp_new = rotate_pts(storesort0x1tp.block(0,0,storesort0x1tp.rows(),2),hatch_angle,x_avg,y_avg);
	storesort0x1tp.block(0,0,storesort0x1tp.rows(),2) = storesort0x1tp_new;
    
    return storesort0x1tp;
}

Eigen::MatrixXd	project_grid_points(Eigen::MatrixXd	fnew, Eigen::MatrixXd v, Eigen::MatrixXd pts, double hatch_angle, double x_avg, double y_avg)
{
	Eigen::MatrixXd p1(1,v.cols());
	Eigen::MatrixXd p2(1,v.cols());
	Eigen::MatrixXd p3(1,v.cols());
	Eigen::MatrixXd	tri(3,v.cols());
	Eigen::MatrixXd	fillpts = Eigen::MatrixXd::Constant(int(1.5*pts.rows()),6,0);
	double a;	
	double b;	
	double c;	
	double d;
	double zval;

	long int fillpts_idx = 0;

	for (long int i=0;i<fnew.rows();++i)
	{
		// vertices for each triangle
		p1 = v.row(fnew(i,0)-1);
		p2 = v.row(fnew(i,1)-1);
		p3 = v.row(fnew(i,2)-1);

		// forming the face with vertices
		tri.block(0,0,1,p1.cols()) = p1;
		tri.block(1,0,1,p1.cols()) = p2;
		tri.block(2,0,1,p1.cols()) = p3;

		// projecting triagles on to the xy plane and 
		// storing all grid points which are inside triangle

		std::vector<int> in = InPoly(tri, pts);

		std::vector<int> loc = find_idx(in);
	
		Eigen::MatrixXd store(loc.size(),pts.cols());
		Eigen::VectorXd storez(loc.size(),1);
		
		for (int loc_i=0;loc_i<loc.size();++loc_i)
		{
			// store.block(loc_i,0,1,pts.cols()) = pts.row(loc[loc_i]);
			store.row(loc_i) = pts.row(loc[loc_i]);
		}
			
		if (store.rows()!=0)
		{
			// creating plane eqaution to get z value of stored points

			a = ((p2(0,1)-p1(0,1))*(p3(0,2)-p1(0,2)))-((p3(0,1)-p1(0,1))*(p2(0,2)-p1(0,2)));
			b = ((p2(0,2)-p1(0,2))*(p3(0,0)-p1(0,0)))-((p3(0,2)-p1(0,2))*(p2(0,0)-p1(0,0)));
			c = ((p2(0,0)-p1(0,0))*(p3(0,1)-p1(0,1)))-((p3(0,0)-p1(0,0))*(p2(0,1)-p1(0,1)));
			d = -(a*p1(0,0))-(b*p1(0,1))-(c*p1(0,2));

			for (int store_i=0;store_i<store.rows();++store_i)
			{
				zval = (-d-(a*store(store_i,0))-(b*store(store_i,1)))/c;
				if (zval<0)
				{
					storez(store_i,0) = 0;
				}
				else
				{
					storez(store_i,0) = zval;	
				}
			}

			// along with points, storing unit normal 
			// associated with each point from which, tcp orientation is calculated.

			Eigen::MatrixXd n = Eigen::MatrixXd::Constant(store.rows(),3,1);
			n.col(0) = n.col(0).array()*(a/sqrt((a*a)+(b*b)+(c*c)));
			n.col(1) = n.col(1).array()*(b/sqrt((a*a)+(b*b)+(c*c)));
			n.col(2) = n.col(2).array()*(c/sqrt((a*a)+(b*b)+(c*c)));

			fillpts.block(fillpts_idx,0,store.rows(),2) = store;
			fillpts.block(fillpts_idx,2,store.rows(),1) = storez;
			fillpts.block(fillpts_idx,3,store.rows(),3) = -n;

			fillpts_idx = fillpts_idx + store.rows();			
		}
	}	

	// giving negative of hatching angle to make points aligned along reference axes

	Eigen::MatrixXd fillpts_final(fillpts_idx,fillpts.cols());	
	fillpts_final = fillpts.block(0,0,fillpts_idx,fillpts.cols());

	Eigen::MatrixXd fillpts_new = rotate_pts(fillpts_final.block(0,0,fillpts_final.rows(),2),-hatch_angle,x_avg,y_avg);
	fillpts_final.col(0) = fillpts_new.col(0);
	fillpts_final.col(1) = fillpts_new.col(1);
	
	return fillpts_final; 
}