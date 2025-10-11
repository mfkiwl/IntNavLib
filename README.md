# IntNavLib

IntNavLib is a C++ integrated navigation library with a focus on readability, modularity, and documentation.
Integrated navigation example applications using the IntNavLib API, including a ROS2 app, are in the [apps](/apps/) folder.

The code is inspired by the books [Principles of GNSS, Inertial, and Multisensor Integrated Navigation](https://ieeexplore.ieee.org/document/9101092) and [Estimation with applications to tracking and navigation: theory algorithms and software](https://books.google.it/books/about/Estimation_with_Applications_to_Tracking.html?id=xz9nQ4wdXG4C&redir_esc=y).

References: 

```
@BOOK{gpintnav,
  author={Groves, Paul},
  booktitle={Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems, Second Edition},
  year={2013}
  }

@BOOK{bar2001estimation,
  title={Estimation with applications to tracking and navigation: theory algorithms and software},
  author={Bar-Shalom, Yaakov and Li, X Rong and Kirubarajan, Thiagalingam},
  year={2001},
  publisher={John Wiley \& Sons}
}
```

![image](media/Figure_1.png)

## Documentation

Documentation is enabled by Doxygen-compatible comments. Use the provided Doxyfile to generate project documentation.

```
mkdir -p docs/doxygen
doxygen Doxyfile
firefox docs/doxygen/html/index.html
```

## Build & install instructions

1) Clone the library in your workspace

2) Install Eigen 3. On Ubuntu and other Debian-based Linux distributions, just: `sudo apt install libeigen3-dev`

3) Use the build.sh script to build & install the IntNavLib library locally, and build the apps in the [apps](/apps/) folder.

## Development Guidelines

- Building: see build.sh file

- Linting: see .clang-tidy file

- Testing: write Google tests

- Documentation: document code using Doxygen-compatible comments

- **Hooks**: paste hook.sh to .git/hooks/ so builds, linting, and tests are run before committing

## Using IntNavLib

You can find apps using IntNavLib in the [apps](/apps/) folder.

### ROS 2 App

![image](media/Figure_2.png)

This library also comes with an example of how to use it in a ROS2 project.
To launch a ROS2 node performing integrated loosely-coupled GNSS/INS navigation, first install ROS2 (Humble supported, other distributions might require minor refactoring), then run:

```
cd apps/ros2_lc_ins_gnss_ecef
source install/setup.bash
ros2 launch ros2_lc_ins_gnss_ecef Profile_3_launch.py log_level:=ERROR
```

Find Python and C++ scripts to plot results in the [utils](/utils/) directory. You can launch them like this:

```
python3 plot_profile.py <path_to_profile_csv> # to plot a motion profile
python3 plot_errors.py <path_to_errors_csv> # to plot nav errors
python3 plot_errors_sigmas_ecef.py <path_to_ecef_errors_sigmas_csv> # to plot errors + estimated uncertainty 
```