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

Eigen::VectorXd linsp(double, double, double);
std::vector<int> ismember(Eigen::MatrixXd, Eigen::MatrixXd);
std::vector<int> find_idx(std::vector<int>);
double median(std::vector<double>);


int main()
{

	//Gap between 2 layers -
	double pathgap_z = 0.25;
	    
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

	std::cout << num_of_layers << std::endl;

	return 0;

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