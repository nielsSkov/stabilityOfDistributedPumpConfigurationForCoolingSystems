## Stability of Distributed Pump Configuration for Cooling Systems

### run_main
This script sets up and runs the non-linear system simulation with linear control, using the function *simSys_nEq4()*

#### simSys_nEq4()
Within this function a python script, *pyMinimize.py*, is called to solve an optimization problem. Some required python setup is found in the top of *run_main*.

### run_reducedLinSys
This script sets up and runs the linear system with state feedback control, using the function *reducedLinSys()*
