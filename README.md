# GraphGNSSLib
### An Open-source Package for GNSS Positioning and Real-time Kinematic Using Factor Graph Optimization

This repository is the implementation of the open-sourced package, the GraphGNSSLib, which makes use of the factor graph optimization (FGO) to perform the GNSS positioning and real-time kinematic (RTK) positioning. In this package, measurements from the historical and current epochs are structured into a factor graph which is then solved by non-linear optimization. The package is based on C++ which is compatible with the robot operation system (ROS) platform. Meanwhile, this package combines the RTKLIB (**[version: 2.4.3 b33](http://www.rtklib.com/)**) to read/decode the GNSS [RINEX](https://en.wikipedia.org/wiki/RINEX) files. Users from Robotics field can easily have access to GNSS raw data for further study.

**Important Notes**: 
  - Be noted that the **GNSS Positioning** mentioned throughout the package means estimating the positioing of the GNSS receiver based on the combination of pseudorange and Doppler measurements uisng FGO.
  - Be noted that the **GNSS-RTK Positioning** mentioned throughout the package means estimating the positioing (float solution) of the GNSS receiver based on the combination of double-differenced pseudorange, carrier-phase and the Doppler measurements using FGO. Finally, the ambiguity is resolved using LAMBDA algorithm.

**Authors**: [Weisong Wen](https://weisongwen.wixsite.com/weisongwen), [Li-ta Hsu](https://www.polyu-ipn-lab.com/) from the [Intelligent Positioning and Navigation Laboratory](https://www.polyu-ipn-lab.com/), The Hong Kong Polytechnic University. 

**Related Papers:** (paper is not exactly same with code)
  - Wen, W., & Hsu, L. T. (2021, May). [Towards robust GNSS positioning and Real-time kinematic using factor graph optimization](https://ieeexplore.ieee.org/abstract/document/9562037). In 2021 IEEE International Conference on Robotics and Automation (ICRA) (pp. 5884-5890). IEEE. 
  - Wen, W., Zhang, G., & Hsu, L. T. (2021). [GNSS outlier mitigation via graduated non-convexity factor graph optimization](https://ieeexplore.ieee.org/abstract/document/9627801). IEEE Transactions on Vehicular Technology, 71(1), 297-310.
  - Zhong, Y., Wen, W., Ng, H. F., Bai, X., & Hsu, L. T. (2022, September). [Real-time Factor Graph Optimization Aided by Graduated Non-convexity Based Outlier Mitigation for Smartphone Decimeter Challenge](https://www.ion.org/publications/abstract.cfm?articleID=18382). In Proceedings of the 35th International Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS+ 2022) (pp. 2339-2348).

*if you use GraphGNSSLib for your academic research, please cite our related [papers](https://ieeexplore.ieee.org/abstract/document/9562037)*

<p align="center">
  <img width="712pix" src="img/software_flowchart.png">
</p>

<center> Software flowchart of GraphGNSSLib, more information please refer to mannual and paper.</center>

## 0. Docker support
If you are not familiar with ROS, we highly recommend using our docker container to enjoy GraphGNSSLib. 
For the details, please go to the branch [test_docker](https://github.com/weisongwen/GraphGNSSLib/tree/test_docker) step 5.



## 1. Prerequisites
### 1.1 **Ubuntu** and **ROS**
Ubuntu 64-bit 16.04, ROS Kinetic. [ROS Installation](http://wiki.ros.org/ROS/Installation). We only test it on Ubuntu 16.04 with ROS Kinetic. 

### 1.2. **Ceres Solver**
Follow the following instructions to install Ceres-solver instead of using the latest version of Ceres-solver.

**Step 1**: Download the [Ceres-solver](https://github.com/weisongwen/GraphGNSSLib/tree/master/support_files) which is compatible with GraphGNSSLib. 

**Step 2**: make and install
```bash
sudo apt-get install cmake
# google-glog + gflags
sudo apt-get install libgoogle-glog-dev
# BLAS & LAPACK
sudo apt-get install libatlas-base-dev
# Eigen3
sudo apt-get install libeigen3-dev
# make Ceres-solver
mkdir ceres-bin
cd ceres-bin
cmake ../ceres-solver
sudo make -j4
sudo make test
sudo make install
```

### 1.3. **Extra Libraries**
```bash
sudo apt-get install ros-kinetic-novatel-msgs
```
## 2. Build GraphGNSSLib
Clone the repository and catkin_make:
```bash
mkdir GraphGNSSLib/src
cd ~/GraphGNSSLib/src
mkdir result
git clone https://github.com/weisongwen/GraphGNSSLib.git
cd ../
# if you fail in the last catkin_make, please source and catkin_make again
catkin_make
source ~/GraphGNSSLib/devel/setup.bash
catkin_make
```
(**if you fail in this step, try to find another computer with clean system or reinstall Ubuntu and ROS**)

## 3. Run GNSS positioning via FGO using dataset [UrbanNav](https://www.polyu-ipn-lab.com/download)   
The GNSS positioning via FGO is validated using static dataset collected near TST of Hong Kong. Several parameters are as follows:
  - GPS second span: **46701** to **47185**
  - satellite system: **GPS/BeiDou**
  - Window Size: **Batch**
  - measurements considered: double-differenced pseudorange and carrier-phase measurements, Doppler measurements
  - result is saved by default
    ```c++
    ~/GraphGNSSLib/trajectory_psr_dop_fusion.csv
    ```

please enable the following in rtklib.h
```bash
#define RTK_FGO 0
```
- Solution 1 to run the GNSS positioning Demo
  ```bash
  source ~/GraphGNSSLib/devel/setup.bash
  # read GNSS raw data and publish as ROS topic
  # we provide several datasets, enjoy it!
  roslaunch global_fusion dataublox_TST20190428.launch
  # run pseudorange and doppler fusion
  roslaunch global_fusion psr_doppler_fusion.launch
  ```
<p align="center">
  <img width="712pix" src="img/SPP_trajectory1.png">
</p>
<center> Trajectories of three methods (GNSS positioning using WLS with the red curve, GNSS positioning using EKF with the green curve, and GNSS positioning using FGO with blue curve throughout the test. The x-axis and y-axis denote the east and north directions, respectively</center>

<p align="center">
  <img width="712pix" src="img/TSTData.gif">
</p>
<center> TST data collected by the Smartphone: Red dots from WLS, purple curve from FGO</center>

## 4. Run GNSS RTK-FGO using static dataset   
The GNSS RTK-FGO is validated using static dataset collected near TST of Hong Kong. Several parameters are as follows:
  - GPS second span: **270149** to **270306**
  - satellite system: **GPS/BeiDou**
  - Window Size: **Batch**
  - measurements considered: double-differenced pseudorange and carrier-phase measurements, Doppler measurements
  - result is saved by default
    ```c++
    ~/GraphGNSSLib/FGO_trajectoryllh_pdrtk.csv
    ```

please enable the following in rtklib.h
```bash
#define RTK_FGO 1
```
- Solution 1 to run the RTK-FGO Demo
  ```bash
  source ~/GraphGNSSLib/devel/setup.bash
  # read GNSS raw data and publish as ROS topic
  roslaunch global_fusion dataublox_TST20200603.launch
  # run GNSS RTK
  roslaunch global_fusion psr_doppler_car_rtk.launch
  ```
<p align="center">
  <img width="712pix" src="img/RTK_trajectory.png">
</p>
<center> Trajectories of three methods (RTK-EKF with the red dots and RTK-FGO with the blue dots throughout the test. The x-axis and y-axis denote the east and north directions, respectively.</center>


## 5. Docker Support

 To run GraphGNSSLib with docker, first make sure  [docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) are installed on your machine. If you want to use the docker to run the global_fusion:
```bash

cd ~/catkin_ws/src/GraphGNSSLib/docker
make build
sudo -E ./start.bash #Do not delete " -E "
source devel/setup.bash
# run pseudorange and doppler fusion
roslaunch global_fusion psr_doppler_fusion.launch
# you should open another ternimal to enter the docker.
# read GNSS raw data and publish as ROS topic
roslaunch global_fusion dataublox_TST20190428.launch
```

  Also, there is a [video](https://www.youtube.com/watch?v=WMM2de_SxTw) showing the demo after you have built the docker_file in the directory GraphGNSSLib/docker

  If you want to restart the container, please stop it first:
  ```bash
sudo ./stop.bash
#then restart it
sudo -E ./start.bash 
```
 The directory  ~/shared_dir is created to connect the container and the host . In the container, it is located at  ~/graph1/shared_dir, you can also download the code to shared_dir and compile the program in the container (Recommended for those who are interested in making changes to the source code)



## 6. Acknowledgements
We use [Ceres-solver](http://ceres-solver.org/) for non-linear optimization and [RTKLIB](http://www.rtklib.com/) for GNSS data decoding, etc. Some functions are originated from [VINS-mono](https://github.com/HKUST-Aerial-Robotics/VINS-Mono). The [rviz_satellite](https://github.com/nobleo/rviz_satellite) is used for visualization. If there is any thing inappropriate, please contact me through 17902061r@connect.polyu.hk ([Weisong WEN](https://weisongwen.wixsite.com/weisongwen)). Thank you very much for the maintainance by Mr. Zhong Yihan. 

## 7. License
The source code is released under [GPLv3](http://www.gnu.org/licenses/) license. We are still working on improving the code reliability. For any technical issues, please contact Weisong Wen <17902061r@connect.polyu.hk>. For commercial inquiries, please contact Li-ta Hsu <lt.hsu@polyu.edu.hk>.
