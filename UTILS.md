# MUMPS

For MUMPS general sparse linear solver, install:
* `sudo apt-get install intel-mkl-full libmetis-dev libscotch-dev openmpi-bin libopenmpi-dev`
* use `-lmkt_rt` for blas and lapack

# Using Fatrop in another CMake project

* add `find_package(fatrop)` to your project's CMakeLists.txt file (or in one of the subdirectories where it is needed)
* add Fatrop to the target link libraries and add the Fatrop include directories: `target_link_libraries(${example} fatrop)` and `target_include_directories(${example} PRIVATE ${fatrop_INCLUDE_DIR})`
* when running `ccmake`, set the fatrop_DIR to your Fatrop directory, the default is `/usr/local/cmake/fatrop`
* build your project, happy solving :-)

# Reserving CPU cores

When benchmarking (or when you want best Fatrop-app performance), you can reserve CPU cores and assign them to the Fatrop-app. This can be done with isolcpu as explained below. An alternative is to use [cset-shield](http://manpages.ubuntu.com/manpages/trusty/man1/cset-shield.1.html). This latter allows for dynamically changing which CPUs are shielded, and also permits the scheduler to (automatically) balance loads in the shielded CPUs, which is often more convenient in multithreaded applications.

_Dedicate a Whole CPU Core to a Particular Program
While taskset allows a particular program to be assigned to certain CPUs, that does not mean that no other programs or processes will be scheduled on those CPUs. If you want to prevent this and dedicate a whole CPU core to a particular program, you can use isolcpus kernel parameter, which allows you to reserve the CPU core during boot.
Add the kernel parameter isolcpus=<CPU_ID> to the boot loader during boot or GRUB configuration file. Then the Linux scheduler will not schedule any regular process on the reserved CPU core(s), unless specifically requested with taskset. For example, to reserve CPU cores 0 and 1, add isolcpus=0,1 kernel parameter. Upon boot, then use taskset to safely assign the reserved CPU cores to your program._

Getting started:

* check the topology of your cpus with `likwid-topology -g`
* open /etc/default/grub
* modify line setting `GRUB_CMDLINE_LINUX_DEFAULT` to something like `GRUB_CMDLINE_LINUX_DEFAULT="quiet splash isolcpus=5,11"`, here to isolate cores 5 and 11. Change these numbers to the desired ones.
* save the file and update grub, by running `sudo update-grub`
* reboot your computer
* you can check the isolated cpus through `cat /sys/devices/system/cpu/isolated`
* you can check the affinity list of process with id 1 through `taskset -cp 1`. This should now exclude the CPUs that you have isolated.
* you can start a new process as follows, to have it running on CPU 5: `taskset -c 5 <executable>`