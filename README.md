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

## TODOs

- Finish unit testing & debugging ECEF inertial navigation

- Integrate GNSS (LC)

- Integrate GNSS (TC)

- Remove dependencies: place eigen in include dir, glog?

## Future development

- Fuse map matching

- Fuse Lidar

- Fuse CVN