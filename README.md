# 3D Fermat Spiral for [Husky](https://github.com/VccRobot/husky)


## Install
On Windows, Using Git command line (https://git-scm.com)

    git clone --recursive https://github.com/VccRobot/fermat3d.git
   

## Compile

Using CMake GUI on Windows, the project can be compiled following the standard cmake routine. First, create new directory inside where the repo is downloaded and name it `build`. Then in CMake GUI, set the `source code` to where this repo is downloaded and the `binaries` to the `build` directoy. Hit `Configure` and make sure to specify `Viusal Studio 15 2017` as the generator for this project, then hit `Generate`. This will create a `.sln` file under the `build` directory that you can open using Visual Studio 2017. 

## Run

 In Visual Studio, the project is named `fermat3d_bin`. Make sure to select this project and hit F5 on Release mode. After a while, you should see a window pop up with seperate fermat spirals constructed on a terrain mesh.
