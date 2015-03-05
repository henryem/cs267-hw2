# cs267-hw2

=Mac setup instructions=

Visualization: Download SDL (http://www.libsdl.org/) and the visualization tool (http://www.eecs.berkeley.edu/~penpornk/cs267.spr15/hw2/visualize.tar.gz).

GCC: Install Homebrew (http://brew.sh).  Then run:
  brew install gcc
This will install gcc-4.9.  (gcc is still symlinked to the Clang compiler.)

OpenMP: Install GCC as above, then run:
  brew install openmp

CUDA: Install XCode (for Clang):
  xcode-select --install
Then install the NVIDIA CUDA tools, following this document:
  http://developer.download.nvidia.com/compute/cuda/6_5/rel/docs/CUDA_Getting_Started_Mac.pdf
