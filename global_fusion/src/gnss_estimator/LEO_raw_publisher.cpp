#include <ros/ros.h>
#include <std_msgs/Header.h>
#include <nlosExclusion/GNSS_Raw.h>
#include <nlosExclusion/GNSS_Raw_Array.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include "../../RTKLIB/src/rtklib.h"

std::vector<std::vector<std::string>> readCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<std::string>> data;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> row;
        std::string value;
        while (std::getline(iss, value, ',')) {
            row.push_back(value);
        }
        data.push_back(row);
    }

    return data;
}

void processCSVData(const std::vector<std::vector<std::string>>& data, ros::Publisher& pub) {
    nlosExclusion::GNSS_Raw_Array gnss_raw_array;
    gnss_raw_array.header = std_msgs::Header();
    ros::Rate rate(2.5);  // 2.5 Hz

    size_t i = 0;
    size_t sv_count = 0;
    double total_sv_init = 0;

    while (i < data.size()) {
        if (data[i].empty()) {
            break;
        }

        double gnss_time = std::stod(data[i]);
        total_sv_init = std::stod(data[i+1]);
        double sat_info[6][total_sv_init];
        double sat_clk_info[2][total_sv_init];
        double doppler_shifts[total_sv_init];
        double lambdas[total_sv_init];
        int vsat[total_sv_init];
        while (sv_count < total_sv_init) {
            nlosExclusion::GNSS_Raw gnss_raw;
            gnss_raw.GNSS_time = std::stod(data[i]);
            gnss_raw.total_sv = std::stod(data[i+1]);
            gnss_raw.prn_satellites_index = std::stod(data[i+2]);
            gnss_raw.pseudorange = std::stod(data[i+3]);
            gnss_raw.raw_pseudorange = std::stod(data[i+4]);
            gnss_raw.carrier_phase = std::stod(data[i+5]);
            gnss_raw.lamda = std::stod(data[i+6]);
            gnss_raw.snr = std::stod(data[i+7]);
            gnss_raw.elevation = std::stod(data[i+8]);
            gnss_raw.azimuth = std::stod(data[i+9]);
            gnss_raw.err_tropo = std::stod(data[i+10]);
            gnss_raw.err_iono = std::stod(data[i+11]);
            gnss_raw.sat_clk_err = std::stod(data[i+12]);
            gnss_raw.sat_pos_x = std::stod(data[i+13]);
            gnss_raw.sat_pos_y = std::stod(data[i+14]);
            gnss_raw.sat_pos_z = std::stod(data[i+15]);
            gnss_raw.visable = std::stoi(data[i+16]);
            gnss_raw.sat_system = data[i+17];
            gnss_raw.visable3DMA = std::stoi(data[i+18]);
            gnss_raw.prE3dMA = std::stod(data[i+19]);

            gnss_raw_array.GNSS_Raws.push_back(gnss_raw);
                        
            sat_info[0][sv_count] = std::stod(data[i+13]);
            sat_info[1][sv_count] = std::stod(data[i+14]);
            sat_info[2][sv_count] = std::stod(data[i+15]);
            sat_info[3][sv_count] = std::stod(data[i+20]);
            sat_info[4][sv_count] = std::stod(data[i+21]);
            sat_info[5][sv_count] = std::stod(data[i+22]);
            
            sat_clk_info[0][sv_count] = std::stod(data[i+23]);
            sat_clk_info[1][sv_count] = std::stod(data[i+24]);
            
            doppler_shifts[sv_count]=std::stod(data[i+25]);

            azel[0][sv_count]=std::stod(data[i+9])*D2R;
            azel[1][sv_count]=std::stod(data[i+8])*D2R;

            lambdas[sv_count]=std::stod(data[i+6]);

            vsat[sv_count]=1; // all satellites are visible

            sv_count++;
            i += 26;  // read 26 rows each time

            if (sv_count >= total_sv_init) {
                pub.publish(gnss_raw_array);
                gnss_raw_array = nlosExclusion::GNSS_Raw_Array();
                gnss_raw_array.header = std_msgs::Header();
                psr_spp(0, gnss_time, sat_info, sat_clk_info, doppler_shifts, lambdas, vsat);
                double rr[6] = {0.0, 0.0, 0.0,0.0, 0.0, 0.0};
                // will publish LEO velocity here 
                pntpos_LEO(*doppler_shifts, *lambdas,
                      total_sv_init, *sat_info, *sat_clk_info, *rr,
                      *azel, *vsat)          
                sv_count = 0;
                rate.sleep();
            }
        }
    }
}
int main(int argc, char** argv) {
    ros::init(argc, argv, "gnss_raw_publisher");
    ros::NodeHandle nh;
    ros::Publisher pub = nh.advertise<nlosExclusion::GNSS_Raw_Array>("/gnss_preprocessor_node/LEOPsrCarRov1", 25);

    std::string filename = "/home/gao_yixin/GraphGNSSLib/src/GraphGNSSLib/global_fusion/dataset/2021_0521_0607/StarLink_Whampoa_0521.csv";
    std::vector<std::vector<std::string>> data = readCSV(filename);
    processCSVData(data, pub);

    ros::spin();
    return 0;
}