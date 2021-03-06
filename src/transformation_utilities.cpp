#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include "transformation_utilities.hpp"

///////////////////////////////////////////////////////////

Eigen::Matrix4i rtf::hom_T(Eigen::Vector3i t, Eigen::Matrix3i r)
{
        Eigen::Matrix4i T;
        T.block(0,3,3,1) << t;
        T.block(0,0,3,3) << r;
        T.block(3,0,1,4) << 0,0,0,1;
        return T;
}

///////////////////////////////////////////////////////////

Eigen::Matrix4f rtf::hom_T(Eigen::Vector3f t, Eigen::Matrix3f r)
{
        Eigen::Matrix4f T;
        T.block(0,3,3,1) << t;
        T.block(0,0,3,3) << r;
        T.block(3,0,1,4) << 0,0,0,1;
        return T;
}

///////////////////////////////////////////////////////////

Eigen::Matrix4d rtf::hom_T(Eigen::Vector3d t, Eigen::Matrix3d r)
{
        Eigen::Matrix4d T;
        T.block(0,3,3,1) << t;
        T.block(0,0,3,3) << r;
        T.block(3,0,1,4) << 0,0,0,1;
        return T;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::apply_transformation(Eigen::MatrixXd data, Eigen::Matrix4d T_mat)
{
        //NOTE: Homogeneous Tranformation Matrix (4x4)

        // putting data in [x, y, z, 1]' format
        Eigen::MatrixXd data_with_fourth_row(data.cols()+1,data.rows());
        Eigen::VectorXd ones_vec = Eigen::VectorXd::Constant(data.rows(),1);
        data_with_fourth_row.block(0,0,data.cols(),data.rows()) = data.transpose();
        data_with_fourth_row.block(data.cols(),0,1,data.rows()) = ones_vec.transpose();
        Eigen::MatrixXd transformed_data = T_mat*data_with_fourth_row;
        Eigen::MatrixXd transformed_data_mat(transformed_data.rows()-1,transformed_data.cols());
        transformed_data_mat = transformed_data.block(0,0,transformed_data.rows()-1,transformed_data.cols());
        return transformed_data_mat.transpose();
}

////////////////////////////////////////////////////////////

std::string rtf::validate_seq(std::string seq)
{
        if(seq =="")
                seq = "ZYX";
        bool invalid_flag = false;
        if(seq.size()>0 && seq.size()<3)
        {
                invalid_flag = true;
        }
        for (int i =0; i<3; ++i)
                if(seq[i]!='X' && seq[i]!='Y' && seq[i]!='Z')
                {
                        invalid_flag = true;
                        break;
                }
        if(invalid_flag)
        {
                std::cerr << "ERROR: Invalid Rotations Sequence: " << seq << std::endl;
                std::terminate();
        }
        return seq;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::eul2rot(Eigen::MatrixXd eul_angles, std::string seq)
{
        seq = rtf::validate_seq(seq);
        Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
        for (int i=0; i<3; ++i)
        {
                if(seq[i]=='X')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
                else if(seq[i]=='Y')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());
                else if(seq[i]=='Z')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());
        }
        return rot_mat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::rot2eul(Eigen::Matrix3d rot_mat, std::string seq)
{
        seq = rtf::validate_seq(seq);
        int rot_idx[3];
        for (int i=0; i<3; ++i)
        {
                if(seq[i]=='X')
                        rot_idx[i] = 0;
                else if(seq[i]=='Y')
                        rot_idx[i] = 1;
                else if(seq[i]=='Z')
                        rot_idx[i] = 2;
        }
        Eigen::MatrixXd eul_angles(1,3);
        Eigen::Vector3d eul_angles_vec;
        eul_angles_vec = rot_mat.eulerAngles(rot_idx[0], rot_idx[1], rot_idx[2]);
        eul_angles(0,0) = eul_angles_vec[0];
        eul_angles(0,1) = eul_angles_vec[1];
        eul_angles(0,2) = eul_angles_vec[2];
        return eul_angles;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::qt2rot(Eigen::MatrixXd quat)
{
        Eigen::Quaterniond q;
        q.x() = quat(0,0);
        q.y() = quat(0,1);
        q.z() = quat(0,2);
        q.w() = quat(0,3);
        return q.normalized().toRotationMatrix();
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::rot2qt(Eigen::Matrix3d rot_mat)
{
        Eigen::MatrixXd quat(1,4);
        Eigen::Quaterniond q(rot_mat);
        quat(0,0) = q.x();
        quat(0,1) = q.y();
        quat(0,2) = q.z();
        quat(0,3) = q.w();
        return quat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::eul2qt(Eigen::MatrixXd eul_angles,std::string seq)
{
        seq = rtf::validate_seq(seq);
        Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
        for (int i=0; i<3; ++i)
        {
                if(seq[i]=='X')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
                else if(seq[i]=='Y')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());
                else if(seq[i]=='Z')
                        rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());
        }
        Eigen::MatrixXd quat(1,4);
        Eigen::Quaterniond q(rot_mat);
        quat(0,0) = q.x();
        quat(0,1) = q.y();
        quat(0,2) = q.z();
        quat(0,3) = q.w();
        return quat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::qt2eul(Eigen::MatrixXd quat, std::string seq)
{
        seq = rtf::validate_seq(seq);
        Eigen::Matrix3d rot_mat = rtf::qt2rot(quat);
        int rot_idx[3];
        for (int i=0; i<3; ++i)
        {
                if(seq[i]=='X')
                        rot_idx[i] = 0;
                else if(seq[i]=='Y')
                        rot_idx[i] = 1;
                else if(seq[i]=='Z')
                        rot_idx[i] = 2;
        }
        Eigen::MatrixXd eul_angles(1,3);
        Eigen::Vector3d eul_angles_vec;
        eul_angles_vec = rot_mat.eulerAngles(rot_idx[0], rot_idx[1], rot_idx[2]);
        eul_angles(0,0) = eul_angles_vec[0];
        eul_angles(0,1) = eul_angles_vec[1];
        eul_angles(0,2) = eul_angles_vec[2];
        return eul_angles;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::eul2bxbybz(Eigen::MatrixXd eul_angles)
{
        // input euler alpha,beta,gamma for ZYX (in radians)
        Eigen::MatrixXd bxbybz = Eigen::MatrixXd::Constant(eul_angles.rows(),9,0);
        for (unsigned int i=0; i<eul_angles.rows(); ++i)
        {
                Eigen::Matrix3d R = rtf::eul2rot(eul_angles.row(i));
                bxbybz.row(i) << R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(2,1),R(0,2),R(1,2),R(2,2);
        }
        return bxbybz;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::bxbybz2eul(Eigen::MatrixXd bxbybz)
{
        Eigen::MatrixXd Gx(1,3), Gy(1,3), Gz(1,3);
        Gx << 1,0,0;
        Gy << 0,1,0;
        Gz << 0,0,1;
        Eigen::Matrix3d R = Eigen::Matrix3d::Constant(0);
        Eigen::MatrixXd cba = Eigen::MatrixXd::Constant(bxbybz.rows(), 3, 0);
        for (unsigned int i=0; i<bxbybz.rows(); ++i)
        {
                R(0,0) = bxbybz(i,0)*Gx(0,0)+bxbybz(i,1)*Gx(0,1)+bxbybz(i,2)*Gx(0,2);
                R(0,1) = bxbybz(i,3)*Gx(0,0)+bxbybz(i,4)*Gx(0,1)+bxbybz(i,5)*Gx(0,2);
                R(0,2) = bxbybz(i,6)*Gx(0,0)+bxbybz(i,7)*Gx(0,1)+bxbybz(i,8)*Gx(0,2);
                R(1,0) = bxbybz(i,0)*Gy(0,0)+bxbybz(i,1)*Gy(0,1)+bxbybz(i,2)*Gy(0,2);
                R(1,1) = bxbybz(i,3)*Gy(0,0)+bxbybz(i,4)*Gy(0,1)+bxbybz(i,5)*Gy(0,2);
                R(1,2) = bxbybz(i,6)*Gy(0,0)+bxbybz(i,7)*Gy(0,1)+bxbybz(i,8)*Gy(0,2);
                R(2,0) = bxbybz(i,0)*Gz(0,0)+bxbybz(i,1)*Gz(0,1)+bxbybz(i,2)*Gz(0,2);
                R(2,1) = bxbybz(i,3)*Gz(0,0)+bxbybz(i,4)*Gz(0,1)+bxbybz(i,5)*Gz(0,2);
                R(2,2) = bxbybz(i,6)*Gz(0,0)+bxbybz(i,7)*Gz(0,1)+bxbybz(i,8)*Gz(0,2);
                cba.row(i) = rtf::rot2eul(R,"ZYX");
        }
        return cba;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::get_rob_T_part(Eigen::MatrixXd part_pts, Eigen::MatrixXd rob_pts)
{
// part_pts: co-ordinates with respect to CAD part frame
// rob_pts: points with respect to robot base frame (or world frame if it matches with robot base frame)
// input: part_pts = [x1, y1, z1;
//                             x2, y2, z2;
//                                :
//                                :
//                             xn, yn, zn]
// input: rob_pts = [x1, y1, z1;
//                            x2, y2, z2;
//                                :
//                                :
//                            xn, yn, zn]
        Eigen::MatrixXd centroid_part_pts(1,part_pts.cols());
        Eigen::MatrixXd centroid_rob_pts(1,rob_pts.cols());
        Eigen::MatrixXd shifted_part_pts(part_pts.rows(),part_pts.cols());
        Eigen::MatrixXd shifted_rob_pts(rob_pts.rows(),rob_pts.cols());
        Eigen::MatrixXd cros_cov_mat(part_pts.cols(),rob_pts.cols());
        Eigen::Matrix3d R;
        Eigen::Matrix3d U_T;
        Eigen::Matrix3d V;
        Eigen::Vector3d T;
        Eigen::Matrix3d M = Eigen::Matrix3d::Identity();
        Eigen::Matrix4d transform_mat = Eigen::Matrix4d::Constant(0);
        if (part_pts.rows()==rob_pts.rows())
        {
                centroid_part_pts = rtf::mean(part_pts);
                centroid_rob_pts = rtf::mean(rob_pts);
                shifted_part_pts.col(0) = part_pts.col(0).array() - centroid_part_pts(0,0);
                shifted_part_pts.col(1) = part_pts.col(1).array() - centroid_part_pts(0,1);
                shifted_part_pts.col(2) = part_pts.col(2).array() - centroid_part_pts(0,2);
                shifted_rob_pts.col(0) = rob_pts.col(0).array() - centroid_rob_pts(0,0);
                shifted_rob_pts.col(1) = rob_pts.col(1).array() - centroid_rob_pts(0,1);
                shifted_rob_pts.col(2) = rob_pts.col(2).array() - centroid_rob_pts(0,2);
                cros_cov_mat = shifted_part_pts.transpose()*shifted_rob_pts;

                // Singular Value Decomposition
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(cros_cov_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);

                // Take care of reflection case due to negative eigen vectors
                U_T = svd.matrixU().transpose();    V = svd.matrixV();
                M(2,2) = (V*U_T).determinant();
                R = V*M*U_T;
                if (R.determinant()>0)
                {
                        T = -R*centroid_part_pts.transpose() + centroid_rob_pts.transpose();
                        transform_mat.block(0,0,3,3) = R;
                        transform_mat.block(0,3,3,1) = T;
                        transform_mat(3,3) = 1;
                }
                else
                {
                        std::cerr << "ERROR: Determinant of rotation matrix is negative..." << std::endl;
                }
        }
        else
        {
                std::cerr << "ERROR: FUNCTION ERROR: For correspondance, number of rows of both matrices should be same..." << std::endl;
        }
        return transform_mat;
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::mean(Eigen::MatrixXi mat)
{
        Eigen::VectorXd vec(mat.cols());
        for (int i=0; i<mat.cols(); ++i)
        {
                vec(i) =  mat.block(0,i,mat.rows(),1).sum()/mat.rows();
        }
        return vec.transpose();
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::mean(Eigen::MatrixXf mat)
{
        Eigen::VectorXd vec(mat.cols());
        for (int i=0; i<mat.cols(); ++i)
        {
                vec(i) =  mat.block(0,i,mat.rows(),1).sum()/mat.rows();
        }
        return vec.transpose();
}

///////////////////////////////////////////////////////////

Eigen::MatrixXd rtf::mean(Eigen::MatrixXd mat)
{
        Eigen::VectorXd vec(mat.cols());
        for (int i=0; i<mat.cols(); ++i)
        {
                vec(i) =  mat.block(0,i,mat.rows(),1).sum()/mat.rows();
        }
        return vec.transpose();
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_x(double t)
{
        Eigen::Matrix3d rx;
        rx <<   1,    0,    0,
        0, cos(t),-sin(t),
        0, sin(t), cos(t);
        return rx;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_y(double t)
{
        Eigen::Matrix3d ry;
        ry << cos(t), 0, sin(t),
        0, 1,    0,
        -sin(t), 0, cos(t);
        return ry;
}

///////////////////////////////////////////////////////////

Eigen::Matrix3d rtf::rot_z(double t)
{
        Eigen::Matrix3d rz;
        rz << cos(t),-sin(t), 0,
        sin(t), cos(t), 0,
        0,      0, 1;
        return rz;
}

///////////////////////////////////////////////////////////
std::string rtf::generate_tool_path(std::string path){
///////////////////////////////////////////////////////
////////////////////// INPUTS /////////////////////////
////////////////////////////////////////////////////////

        // if true generate all layers, else generate only one layer
        bool generate_all_layers = false;
        bool generate_cont_toolpath = true; // add start and end point for each layer

        // Gap between 2 hatching lines -
        double pathgap_x = 1.0;
        double pathgap_y = 1.0;

        //Gap between 2 layers -
        double pathgap_z = 0.5;

        // grid aadition (NOTE: higher number for higher aspect ratio of part in xy plane)
        double grid_addition = 70.0;

        // hatching angle properties
        double start_hatch_angle = 0; //3.14159/4;
        double hatch_angle_change = 3.14159/4; // consecutive layers will have this much change in hatching angle (degrees)

        // Path to generate -
        // 1 - Boundary
        // 2 - Hatching
        int Path_Number = 2; // put number from above
        int FlipTravel = 0; // 1 for yes, 0 for No
        int space = 2; // number of points to skip for smoother path

        ///////////////////////////
        // Bottom layer generation
        ///////////////////////////

        Eigen::MatrixXd v;
        Eigen::MatrixXd f;
        Eigen::MatrixXd n;
        igl::readSTL(path, v, f, n);
        f = f.array()+1; //NOTE : adding 1 make make start index as 1... this is compatible for rest of the code
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

        int num_of_layers = 0;
        if (generate_all_layers)
        {
                num_of_layers = NPAM::number_of_layers(v,f,n,pathgap_z);
        }
        else
        {
                num_of_layers = 1;
        }

        //////////////////////////////////////////////////////
        // identifying the co-ordinates of bottom layer first
        //////////////////////////////////////////////////////

        Eigen::MatrixXd fnew = NPAM::identify_bottom_layer(v,f,n);

        ////////////////////////////
        // generating the tool path
        ////////////////////////////

        Eigen::MatrixXd pts;
        double hatch_angle;
        double x_avg;
        double y_avg;
        long int tool_path_strt_idx = 0;
        Eigen::MatrixXd tool_path(1,1);

        for (int layer=1; layer<=num_of_layers; ++layer)
        {
                hatch_angle = start_hatch_angle + (layer-1)*hatch_angle_change;
                pts = NPAM::generate_grid_points(pathgap_x,pathgap_y,xmin,ymin,xmax,ymax,hatch_angle);

                // apply rotation to points
                x_avg = pts.block(0,0,pts.rows(),1).sum()/pts.rows();
                y_avg = pts.block(0,1,pts.rows(),1).sum()/pts.rows();
                Eigen::MatrixXd rotated_pts = NPAM::rotate_pts(pts,hatch_angle,x_avg,y_avg);

                ////////////////////////////
                // project grid points
                ////////////////////////////
                fillpts = NPAM::project_grid_points(fnew, v, rotated_pts, hatch_angle, x_avg, y_avg);

                ////////////////////////////
                // infill path
                ////////////////////////////
                Eigen::MatrixXd layer_tool_path = NPAM::Infill_Path_with_euler(fillpts, FlipTravel, space, hatch_angle, x_avg, y_avg);
                layer_tool_path.col(2) = layer_tool_path.col(2).array() + ((layer-1)*pathgap_z);

                if (generate_cont_toolpath)
                {
                        ////////////////////////////////
                        // Making continunous tool_path:
                        // adding start and end points
                        ////////////////////////////////
                        tool_path.conservativeResize(layer_tool_path.rows()+tool_path_strt_idx+2,layer_tool_path.cols());
                        tool_path.block(tool_path_strt_idx,0,(layer_tool_path.rows()+2),layer_tool_path.cols()) = Eigen::MatrixXd::Constant((layer_tool_path.rows()+2),layer_tool_path.cols(),0);
                        Eigen::MatrixXd start_pt(1,layer_tool_path.cols());
                        start_pt = layer_tool_path.row(0);
                        start_pt(0,2) = start_pt(0,2)+50;
                        Eigen::MatrixXd end_pt(1,layer_tool_path.cols());
                        end_pt = layer_tool_path.row(layer_tool_path.rows()-1);
                        end_pt(0,2) = end_pt(0,2)+50;
                        tool_path.block(tool_path_strt_idx,0,1,layer_tool_path.cols()) = start_pt;
                        ++tool_path_strt_idx;
                        tool_path.block(tool_path_strt_idx,0,layer_tool_path.rows(),layer_tool_path.cols()) = layer_tool_path;
                        tool_path_strt_idx = tool_path_strt_idx + layer_tool_path.rows();
                        tool_path.block(tool_path_strt_idx,0,1,layer_tool_path.cols()) = end_pt;
                        ++tool_path_strt_idx;
                }
                else
                {
                        ////////////////////////////////
                        // Making toolpath layer by layer
                        ////////////////////////////////
                        tool_path.conservativeResize(layer_tool_path.rows()+tool_path_strt_idx,layer_tool_path.cols());
                        tool_path.block(tool_path_strt_idx,0,(layer_tool_path.rows()),layer_tool_path.cols()) = Eigen::MatrixXd::Constant((layer_tool_path.rows()),layer_tool_path.cols(),0);
                        tool_path.block(tool_path_strt_idx,0,layer_tool_path.rows(),layer_tool_path.cols()) = layer_tool_path;
                        tool_path_strt_idx = tool_path_strt_idx + layer_tool_path.rows();
                }
        }
        std::string out_file=path+".csv";
        file_rw::file_write(out_file,tool_path);

        return out_file;
}
