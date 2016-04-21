# MultistateDesignOptimization

#### 4/13:

* Use Armadillo unless someone says otherwise (could try Eigen but my personal experience with it was convoluted)
  * get libraries working locally and on clusters
* Implement low-hanging fruit first
* Divide and conquor the C++ implementation
  * take care of the bits and pieces and things first
* MapReduce? Sorting of objects
* Other libraries? Mantella or other optimized optimization libraries
* Long-term: load balancing and where bottlenecks happens
  * scheduling
* Serial optimization of calculations]
* Draw up program structure description
* More ways to optimize the program

#### 4/20:

* I started the project in Visual Studio, hence the nested folders, project files, and stdafx precomiled header
* However, I see no reason why it wouldn't work straight outta the box on other machines since we're not using any OS-specific libraries
