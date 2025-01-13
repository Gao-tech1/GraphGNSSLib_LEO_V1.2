#!/usr/bin/env python
# Author: GAO Yixin yixin.gao@connect.polyu.hk
# Date: 2025/1/6

import rospy
import csv
import os
from std_msgs.msg import Header
from nlosExclusion.msg import GNSS_Raw, GNSS_Raw_Array

def read_csv_and_publish(csv_file, pub):
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)

        current_time = None
        gnss_raw_array = None
        sv_count = 0
        rate = rospy.Rate(2.5)  # 2.5 Hz

        gnss_raw_array = GNSS_Raw_Array()
        gnss_raw_array.header = Header()
        start_seconds = 0
        nanoseconds = 0
        total_sv_init=0
        current_stamp = rospy.Time(start_seconds, nanoseconds)
        gnss_raw_array.header.stamp = current_stamp
        gnss_raw_array.header.frame_id = ''
        while True:
            
            rows = [next(reader, None) for _ in range(26)]
            print("Read rows: {}".format(rows))  # Print the content of rows for debugging
            if rows[0] is None:
                break
            gnss_time = float(rows[0][0])
            total_sv_init = float(rows[1][0])
            while sv_count<total_sv_init:
            #if sv_count< total_sv:
                # # Publish the previous GNSS_Raw_Array message
                # if gnss_raw_array is not None:
                #     pub.publish(gnss_raw_array)
                # # Update the current time
                # current_time = gnss_time
                # Create a new GNSS_Raw_Array message
            
                    gnss_raw = GNSS_Raw()
                    gnss_raw.GNSS_time = float(rows[0][0])
                    gnss_raw.total_sv = float(rows[1][0])
                    gnss_raw.prn_satellites_index = float(rows[2][0])
                    gnss_raw.pseudorange = float(rows[3][0])
                    gnss_raw.raw_pseudorange = float(rows[4][0])
                    gnss_raw.carrier_phase = float(rows[5][0])
                    gnss_raw.lamda = float(rows[6][0])
                    gnss_raw.snr = float(rows[7][0])
                    gnss_raw.elevation = float(rows[8][0])
                    gnss_raw.azimuth = float(rows[9][0])
                    gnss_raw.err_tropo = float(rows[10][0])
                    gnss_raw.err_iono = float(rows[11][0])
                    gnss_raw.sat_clk_err = float(rows[12][0])
                    gnss_raw.sat_pos_x = float(rows[13][0])
                    gnss_raw.sat_pos_y = float(rows[14][0])
                    gnss_raw.sat_pos_z = float(rows[15][0])
                    gnss_raw.visable = int(rows[16][0])
                    gnss_raw.sat_system = rows[17][0]
                    gnss_raw.visable3DMA = int(rows[18][0])
                    gnss_raw.prE3dMA = float(rows[19][0])

                    gnss_raw_array.GNSS_Raws.append(gnss_raw)
                    sv_count += 1
                    if sv_count< total_sv_init-1:
                        rows = [next(reader, None) for _ in range(26)]
                    else:
                        # Publish the previous GNSS_Raw_Array message
                        pub.publish(gnss_raw_array)
                        gnss_raw_array = GNSS_Raw_Array()
                        gnss_raw_array.header = Header()
                         # current_stamp=current_stamp + rospy.Duration(1)
                        # gnss_raw_array.header.stamp = current_stamp
                        gnss_raw_array.header.frame_id = ''
                        sv_count=0
                        rate.sleep()


        # Publish the last GNSS_Raw_Array message
        #if gnss_raw_array is not None:
        #    pub.publish(gnss_raw_array)

def main():
    rospy.init_node('gnss_raw_publisher', anonymous=True)
    pub = rospy.Publisher('/gnss_preprocessor_node/LEOPsrCarRov1', GNSS_Raw_Array, queue_size=25)
    # add by Yixin, Change to your csv
    csv_file = os.path.expanduser('~/GraphGNSSLib/src/GraphGNSSLib/global_fusion/dataset/2021_0521_0607/StarLink_Whampoa_0521.csv')
    read_csv_and_publish(csv_file, pub)
    # file_processed = False
    # while not rospy.is_shutdown():
    #     read_csv_and_publish(csv_file, pub)
    #     rate.sleep()

if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass