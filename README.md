# FATROP

FATROP stands for Fast Trajectory Optimizer. 
It is an efficient and reliable solver for nonlinear optimal control problems with stagewise constraints, aimed at online applications.

For MUMPS general sparse linear solver
install 
sudo apt-get install intel-mkl-full libmetis-dev libscotch-dev openmpi-bin libopenmpi-dev
use -lmkt_rt for blas and lapack

When benchmarking (or when you want best Fatrop-app performance), you can reserve CPU cores and assign them to the Fatrop-app:

_Dedicate a Whole CPU Core to a Particular Program
While taskset allows a particular program to be assigned to certain CPUs, that does not mean that no other programs or processes will be scheduled on those CPUs. If you want to prevent this and dedicate a whole CPU core to a particular program, you can use isolcpus kernel parameter, which allows you to reserve the CPU core during boot.
Add the kernel parameter isolcpus=<CPU_ID> to the boot loader during boot or GRUB configuration file. Then the Linux scheduler will not schedule any regular process on the reserved CPU core(s), unless specifically requested with taskset. For example, to reserve CPU cores 0 and 1, add isolcpus=0,1 kernel parameter. Upon boot, then use taskset to safely assign the reserved CPU cores to your program._

