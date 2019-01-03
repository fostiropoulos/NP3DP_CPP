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
	std::string data_dir = "/home/aniruddha/Desktop/Non_Planar_3D_Printing/CPP/data/"; 
    
	std::string fillpts_file = data_dir + "fillpts.csv";			////// DELETE THIS //////
	Eigen::MatrixXd fillpts;					/////// DELETE THIS //////
	fillpts = file_rw::file_read_mat(fillpts_file);

	double space = 3;
	double hatch_angle = 90*3.14159/180;
	bool FlipTravel = 0;

	////////////////////////////
    // infill path
    // MOVE THIS COMPLETE DATA TO A FUNCTION
	////////////////////////////

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

	Eigen:MatrixXd storesort0x1tp_new = rotate_pts(storesort0x1tp.block(0,0,storesort0x1tp.rows(),2),hatch_angle,x_avg,y_avg);
	storesort0x1tp.block(0,0,,storesort0x1tp.rows(),2) = storesort0x1tp_rev;
	return 0;
}	