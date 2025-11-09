## TODOs

### Core Development

- ros2 app: do fast imu prop so state is not delayed.

- Visual Odometry (simplified groves/px4 style)

- Multiple hypotheses 

- Particle filter

- real time whiteness check

- better eval: display NEES, NIS, MAE, RMSE

- generic core templated kalman update func

- monte carlo bash script + eval script

- in-house trajectory generator

- no dynamic memory allocation

- Joseph update

### Support & Infrastructure

- profile code with valgrind

- write unit tests + test bash script

- Remove dependencies: place eigen in include dir as submodule

- clang tidy: cleanup code

- doxygen: add more latex to documentation, complete docs stubs

- why tests fail with float build? But apps seem to work.

# DONE

- template double float filter funcs and classes

- modular update functions, no duplications

- Git hooks

- no matrix inversions!