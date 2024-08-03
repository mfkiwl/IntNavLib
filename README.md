# IntNavLib

This is my toy c++ integrated navigation library.
The aim is to have modular pieces of c++ code that can be used in integrated navigation applications.

![image](media/Figure_1.png)

The code is inspired by [Principles of GNSS, Inertial, and Multisensor Integrated Navigation](https://ieeexplore.ieee.org/document/9101092) by Paul Groves.

Reference: 

```
@BOOK{gpintnav,
  author={Groves, Paul},
  booktitle={Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems, Second Edition},
  year={2013}
  }
```

## Build & install instructions

1) Clone the library in your workspace

2) Install Eigen 3. On Ubuntu and other Debian-based Linux distributions, just: `sudo apt install libeigen3-dev`

3) Build & install the Cmake project:
```
cd intnavlib
mkdir build
cd build
cmake ..
make install
```

## Using IntNavLib

You can find apps using IntNavLib in the [apps](/apps/) directory.

For example, to build and run the inertial navigation demo using ECEF navigation equations, run: 

```
cd apps/ins_ecef
mkdir build && cd build
cmake ..
make -j4
./ins_ecef <path_to>/data/Profile_3.csv
```

Find Python scripts to plot results in Groves' style in the [utils](/utils/) directory. You can launch them like this:

```
python3 plot_profile.py <path_to_profile_csv> # to plot a motion profile
python3 plot_errors.py <path_to_errors_csv> # to plot nav errors
python3 plot_errors_sigmas_ecef.py <path_to_ecef_errors_sigmas_csv> # to plot errors + estimated uncertainty 

```


## TODOs

- tic toc apps

- profile code with valgrind

- add gnss sensor model

- GNSS INS (LC)

- GNSS INS (TC)

- Remove dependencies: place eigen in include dir

## Future development

- Fuse map matching

- Fuse Lidar

- Fuse CVN