# isingmodel
Fortran implementation of the 2D Ising model w/ Metropolis' algorithm and periodic boundary conditions (PBC). I have uploaded the naive implementation as well as a faster version.

**Naive implementation.**
This approximation to a better solution checks if each node is a boundary node to ensure PBC. This is a really slow approach. As a rough approximation: a 50x50 lattice takes an average of 10834 milliseconds per Monte Carlo iteration, which roughly implies 4.3 milliseconds per node. Since this problem shows a time complexity of O(n^2), being n the number of nodes in each side of the square lattice, it is easy to see that runtime will escalate rapidly. For a 128x128 matrix and 1000 Monte Carlo iterations, it would take approximately 18-19 hours to complete the task. These estimations are based on my computer's performance, obviously. Therefore, a quicker approach needs to be implemented. I am currently working on it. 

The repository currently includes the program itself (with a suitable random number generator included) and an R script that I wrote for my own use. The R script needs to be tweaked if the size of the lattice changes.


*Update (31/03/2020)*
**Faster implementation.**
After following my professor's advice, the program now uses 1D-indexing (the lattice is stored in a 1D-array and each node has a unique number ranging from 1 to the size of the lattice squared). Also, the bottleneck caused by checking the neighbours of each node has been removed by calculating the neighbours only once. Comparatively, my program now takes about forty seconds to calculate 1000 MC steps for a 128x128 lattice instead of the estimated 18 hours I was talking about earlier.

I have maintained the first version just like it was for everyone to be able to check both. The native implementation should be good enough for pedagogical purposes, in case anyone ends up here after a Google search.
