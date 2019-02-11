#include "transformation_utilities.hpp"

int main(int argc,char* argv[])
{

        if(argc>=2) {
                std::string stl_file_name = argv[1];
                rtf::generate_tool_path(stl_file_name);
        }else{

                std::cout << "Please pass tool path. ";
                return 0;
        }
        return 0;
}
