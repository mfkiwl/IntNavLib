# IntNavLib

This is my c++ integrated navigation library.
The aim is to have modular pieces of c++ code that can e used in integrated navigation applications.

Code is still WIP and not yet a true library.

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