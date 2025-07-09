# Safety-Critical Ultra-Wideband 3D Localization with Set-Membership Uncertainty Representation [RA-L'25]

## 1. System Overview

![](./pics/overview.png)

## 2. Prerequisites

### 2.1 Eigen3



### 2.2 G2O
[compatible version here](https://github.com/RainerKuemmerle/g2o/tree/9b41a4ea5ade8e1250b9c1b279f3a9c098811b5a)


## 3. Build
```
cd your_ros_ws
mkdir src
cd src
git clone https://github.com/Zhu-YQ/DB-SMF-UWB.git
cd ..
catkin_make
```


## 4. Run
```
source devel/setup.bash
roslaunch db_smf_uwb run_sim.launch
```


## 5. Citation
