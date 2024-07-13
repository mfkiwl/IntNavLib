# IntNavLib

This is my toy c++ integrated navigation library.
The aim is to have modular pieces of c++ code that can be used in integrated navigation applications.

The code is inspired by "Principles of GNSS, Inertial, and Multisensor Integrated Navigation" by Paul Groves.

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

## TODOs

- Integrate GNSS (LC)

- Integrate GNSS (TC)

- Remove dependencies: place eigen in include dir, glog?

## Future development

- Fuse map matching

- Fuse Lidar

- Fuse CVN