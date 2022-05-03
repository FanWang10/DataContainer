#include <iostream>
#include "Grid.hpp"
#include "Solution.hpp"
#include "Field.hpp"
#include <netcdf>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;



double calcultae_function(double x, double y, double z, double a, double b)
{

    return (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z -1) * (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z -1) * (x * x + ((1 + b) * y) * ((1 + b) * y) + z * z - 1) - 
                            (x * x) * (z * z * z) - a * (y * y) * (z * z * z);
;
}

int main()
{

    /*
    **  Task 1
    **
    **
    */

    /*Task 1 Grid Memeber Vectors*/
    // Task 1 Grid Min Vector
    vector<double> min_grid_task_1 {-2, -2, -2};
    // Task 1 Grid Max Vector
    vector<double> max_grid_task_1 {2, 2, 2};
    // Task 1 Num of pts Vector
    vector<int> num_of_pts_grid_task_1 {512, 512, 512};
    // Task 1 Delta Vector
    vector<double> delta_grid_task_1 {4.0 / 511.0, 4.0 / 511.0, 4.0 / 511.0}; 

    /*Task 1 Grid*/
    Grid grid_task_1(min_grid_task_1, max_grid_task_1, num_of_pts_grid_task_1, delta_grid_task_1);

    /*Task 1 Solution Member*/
    vector<double> solution_member_task_1(512 * 512 * 512, 0);

    /*Task 1 Solution*/
    Solution solution_task_1(solution_member_task_1, 512, 512, 512);

    /*Pop Solution With Data*/
    double a, b;
    a = 1;
    b = 1;

    double x, y, z;
    x = min_grid_task_1[0];
    y = min_grid_task_1[1];
    z = min_grid_task_1[2];

    for(int i = 0; i <= num_of_pts_grid_task_1[2] - 1; i++)
    {
        for(int j = 0; j <= num_of_pts_grid_task_1[1] - 1; j++)
        {
            for(int k = 0; k <= num_of_pts_grid_task_1[0] - 1; k++)
            {
                solution_task_1.assign(k, j, i, calcultae_function(x, y, z, a, b));
                x += delta_grid_task_1[0];
            }
            x = min_grid_task_1[0];
            y += delta_grid_task_1[1];
        }
        y = min_grid_task_1[1];
        z += delta_grid_task_1[2];
    }
    /*Pop Solution With Data*/

    /*Create Task 1 Field*/
    Field field_task_1(solution_task_1, grid_task_1);

    /*Get Data From The Field*/
    std::vector<double> task_1_a_data_vector = solution_task_1.get_solution_vector();

    /*Task 1 a*/
    netCDF::NcFile dataFile_task1_a("task1_all_volume.nc", netCDF::NcFile::replace);
    auto xDim = dataFile_task1_a.addDim("x", 512);
    auto yDim = dataFile_task1_a.addDim("y", 512);
    auto zDim = dataFile_task1_a.addDim("z", 512);
    auto data_task_1_a = dataFile_task1_a.addVar("data", netCDF::ncDouble, {zDim, yDim, xDim});
    data_task_1_a.putVar(task_1_a_data_vector.data());
    /*Task 1 a*/


    /*Begin Task 1 b*/
    /*Task 1 b Read From task1_random.txt*/
    ifstream task_1_b_file("task1_random.txt", fstream::in);
    task_1_b_file >> a;
    task_1_b_file.ignore(1);
    task_1_b_file >> b;

    /*If a != 1 && b != 1, re-pop container with value*/
    if(a != 1 && b != 1)
    {
            for(int i = 0; i <= num_of_pts_grid_task_1[2] - 1; i++)
            {
                for(int j = 0; j <= num_of_pts_grid_task_1[1] - 1; j++)
                {
                    for(int k = 0; k <= num_of_pts_grid_task_1[0] - 1; k++)
                    {
                        solution_task_1.assign(k, j, i, calcultae_function(x, y, z, a, b));
                        x += delta_grid_task_1[0];
                    }
                    x = min_grid_task_1[0];
                    y += delta_grid_task_1[1];
                }
                y = min_grid_task_1[1];
                z += delta_grid_task_1[2];
            }   
            field_task_1 = Field(solution_task_1, grid_task_1);
    }
    
    /*Task 1 b Write To File*/
    ofstream task_1_b_output_file("task1_random_value.txt", fstream::out);
    if(task_1_b_output_file.is_open())
    {
        while(!task_1_b_file.eof())
        {
            /* Read From task1_random.txt */
            double x_from_file, y_from_file, z_from_file;
            task_1_b_file >> x_from_file;
            task_1_b_file.ignore(1);
            task_1_b_file >> y_from_file;
            task_1_b_file.ignore(1);
            task_1_b_file >> z_from_file;

            /* Output To task1_random_value.txt*/
            task_1_b_output_file << field_task_1.get_value(x_from_file, y_from_file, z_from_file);
            task_1_b_output_file << " \n";
        }
        
    }
    task_1_b_file.close();
    task_1_b_output_file.close();
    /*Task 1b Finish*/

    /*Begin Task 1 c*/
    /*Task 1 c Read From task1_plane.txt*/
    ifstream task_1_c_file("task1_plane.txt", fstream::in);
    task_1_c_file >> a;
    task_1_c_file.ignore(1);
    task_1_c_file >> b;

    /*If a != 1 && b != 1, re-pop container with value*/
    if(a != 1 && b != 1)
    {
            for(int i = 0; i <= num_of_pts_grid_task_1[2] - 1; i++)
            {
                for(int j = 0; j <= num_of_pts_grid_task_1[1] - 1; j++)
                {
                    for(int k = 0; k <= num_of_pts_grid_task_1[0] - 1; k++)
                    {
                        solution_task_1.assign(k, j, i, calcultae_function(x, y, z, a, b));
                        x += delta_grid_task_1[0];
                    }
                    x = min_grid_task_1[0];
                    y += delta_grid_task_1[1];
                }
                y = min_grid_task_1[1];
                z += delta_grid_task_1[2];
            }   
            field_task_1 = Field(solution_task_1, grid_task_1);
    }

    /*Read x, y length of data vector*/
    int x_length_task1c, y_length_task1c;
    task_1_c_file >> x_length_task1c;
    task_1_c_file.ignore(1);
    task_1_c_file >> y_length_task1c;

    // Create Value Array For Task1c Plane
    double values_task_1_c[y_length_task1c][x_length_task1c];

    // Pop array with data
    for(int i = 0; i < x_length_task1c; i++)
    {
        for(int j = 0; j < y_length_task1c; j++)
        {
            double x_plane, y_plane, z_plane;
            task_1_c_file >> x_plane;
            task_1_c_file.ignore(1);
            task_1_c_file >> y_plane;
            task_1_c_file.ignore(1);
            task_1_c_file >> z_plane;
            values_task_1_c[j][i] = field_task_1.get_value(x_plane, y_plane, z_plane);
            
        }
    }

    // Write To task1_plane.nc
    netCDF::NcFile dataFile_task1_c("task1_plane.nc", netCDF::NcFile::replace);
    auto xDim_1c = dataFile_task1_c.addDim("x", x_length_task1c);
    auto yDim_1c = dataFile_task1_c.addDim("y", y_length_task1c);
    auto data_task_1_c = dataFile_task1_c.addVar("data", netCDF::ncDouble, {yDim_1c, xDim_1c});
    data_task_1_c.putVar(values_task_1_c);
    /*Task 1c Finish*/


    /*
    **  Task 2
    **
    **
    */

    /*Task 2 Grid Memeber Vectors*/
    // Task 2 Grid Min Vector
    vector<double> min_grid_task_2 {0, 0, 0};
    // Task 2 Grid Max Vector
    vector<double> max_grid_task_2 {255, 255, 255};
    // Task 2 Num of pts Vector
    vector<int> num_of_pts_grid_task_2 {256, 256, 256};
    // Task 2 Delta Vector
    vector<double> delta_grid_task_2 {1, 1, 1}; 

    /*Task 2 Grid*/
    Grid grid_task_2(min_grid_task_2, max_grid_task_2, num_of_pts_grid_task_2, delta_grid_task_2);

    /*Task 2 Solution Member*/
    vector<double> solution_member_task_2(256 * 256 * 256, 0);

    /*Task 2 Solution*/
    Solution solution_task_2(solution_member_task_2, 256, 256, 256);

    /*Read From raw And Pop Data Into Solution*/
    string task2_raw_name = "task2.raw";
    solution_task_2.read_from_raw(task2_raw_name);

    /*Create Task 2 Field*/
    Field field_task_2(solution_task_2, grid_task_2);

    /*********Task 2a*********/
    std::vector<double> solution_vector_task2a = solution_task_2.get_solution_vector();
    netCDF::NcFile dataFile_task2_a("task2_all_volume.nc", netCDF::NcFile::replace);
    auto xDim_2a = dataFile_task2_a.addDim("x", 256);
    auto yDim_2a = dataFile_task2_a.addDim("y", 256);
    auto zDim_2a = dataFile_task2_a.addDim("z", 256);
    auto data_2a = dataFile_task2_a.addVar("data", netCDF::ncDouble, {xDim_2a, yDim_2a, zDim_2a});
    data_2a.putVar(solution_vector_task2a.data());
    /*********Task 2a*********/


    /*********Task 2b*********/
    ifstream task_2_b_file("task2_random.txt", fstream::in);
    ofstream task_2_b_output_file("task2_random_value.txt", fstream::out);
    if(task_2_b_output_file.is_open())
    {
        while(!task_2_b_file.eof())
        {
            /* Read From task2_random.txt */
            double x_from_file, y_from_file, z_from_file;
            task_2_b_file >> x_from_file;
            task_2_b_file.ignore(1);
            task_2_b_file >> y_from_file;
            task_2_b_file.ignore(1);
            task_2_b_file >> z_from_file;

            /* Output To task1_random_value.txt*/
            task_2_b_output_file << field_task_2.get_value(x_from_file, y_from_file, z_from_file);
            task_2_b_output_file << " \n";
        }
        
    }
    task_2_b_file.close();
    task_2_b_output_file.close();
    /*********Task 2b*********/

    /*********Task 2c*********/
    ifstream task_2_c_file("task2_plane.txt", fstream::in);
    int x_length_task2c, y_length_task2c;
    task_2_c_file >> x_length_task2c;
    task_2_c_file.ignore(1);
    task_2_c_file >> y_length_task2c;

    // Create Value Array For Task2c Plane
    double values_task_2_c[x_length_task2c][y_length_task2c];
    for(int i = 0; i < x_length_task2c; i++)
    {
        for(int j = 0; j < y_length_task2c; j++)
        {
            double x_plane, y_plane, z_plane;
            task_2_c_file >> x_plane;
            task_2_c_file.ignore(1);
            task_2_c_file >> y_plane;
            task_2_c_file.ignore(1);
            task_2_c_file >> z_plane;
            values_task_2_c[i][j] = field_task_2.get_value(x_plane, y_plane, z_plane);
        }
    }

    netCDF::NcFile dataFile_task2_c("task2_plane.nc", netCDF::NcFile::replace);
    auto xDim_2c = dataFile_task2_c.addDim("x", x_length_task2c);
    auto yDim_2c = dataFile_task2_c.addDim("y", y_length_task2c);
    auto data_task_2_c = dataFile_task2_c.addVar("data", netCDF::ncDouble, {xDim_2c, yDim_2c});
    data_task_2_c.putVar(values_task_2_c);
    /*********Task 2c*********/


   


    return 0;
}