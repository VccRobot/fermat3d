# 3D Fermat Spiral for [Husky](https://github.com/VccRobot/husky)


## Install:
On Windows, Using Git command line (https://git-scm.com)

    git clone --recursive https://github.com/VccRobot/fermat3d.git
   

## Compile:

Using CMake GUI on Windows, the project can be compiled following the standard cmake routine. First, create new directory inside where the repo is downloaded and name it `build`. Then in CMake GUI, set the `source code` to where this repo is downloaded and the `binaries` to the `build` directoy. Hit `Configure` and make sure to specify `Viusal Studio 15 2017` as the generator for this project, then hit `Generate`. This will create a `.sln` file under the `build` directory that you can open using Visual Studio 2017. 

## Run:

 In Visual Studio, the project is named `fermat3d_bin`. Make sure to select this project and hit F5 on Release mode. After a while, you should see a window pop up with seperate fermat spirals constructed on a terrain mesh.
 
 ### Fixing Runtime Error:
 During the construction of isolines, libigl uses a very small number to detect the duplicated vertices. If you receive a run-time error, that means you need to increase this number. You can do this by going to 'fermat3d\libigl\include\igl\isolines.cpp' and change the number in line 106. On my machine, I changed it to 2.2204e-9.  


## Understanding the Code:

Constructing of the Connected Fermat Spiral follows the steps in (Connected Fermat Spirals for Layered Fabrication)[https://www.cs.sfu.ca/~haoz/pubs/zhao_sig16_fermat.pdf]. Instead of using the offsets contours from the boundary of the input shape, we use the isolines from the height function (compute_surface_func.h). The isolines are computed using Libigl. There a set of post-processing operations that needs to be done in order to prepare the isolines to construct the Fermat Spirals. The 
post-processing of the isolines includes 

- Removing zero length isolines 
- Creating separate spirallable regions 

All the operations related to isolines are implemented in isolines_ds.h and isolines_util.h. Since we are using the height function, this creates isolines with non-uniform spacing. The proper way is to compute the offset from the boundary of the input mesh based on the geodesic function as done in (DSCarver: decompose-and-spiral-carve for subtractive manufacturing)[https://www.cs.sfu.ca/~haoz/pubs/zhao_sig18_cnc.pdf] and then use this to compute the isolines. 

After computing the isolines and recognizing the spirallable regions, we are ready to construct the Fermat spiral. For each spirallable regions, we create a single spiral and then convert it to Fermat spiral (fermat3d_imp.h).
