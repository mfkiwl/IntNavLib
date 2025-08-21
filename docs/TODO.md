## TODOs

### Core Development

- better variable names: no stupid random uppercases

- integrated test after done refactoring app

- bad euler / dcm conversion

- avoid duplication in apps, tests. write interface functions / classes to init sensors, run sim

- no matrix inversions!

- ros2 app: do fast imu prop so state is not delayed.

- Visual Odometry (simplified groves/px4 style)

- Multiple hypotheses 

- Particle filter

- real time whiteness check

- isolate ekf update function?

- random error init in apps, w gaussian dist

- plot innovations

- NEES evaluation

- monte carlo

- make single app, pass desired setup by arg

- in-house trajectory generator

### Support & Infrastructure

- profile code with valgrind

- write tests + test bash script

- Remove dependencies: place eigen in include dir as submodule

- clag tidy: cleanup code

- doxygen: add latex to documentation