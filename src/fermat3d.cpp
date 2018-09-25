#include <igl/readOBJ.h>

#include "isolines_ds.h"
#include "fermat3d_imp.h"
#include "compute_surface_func.h"

#ifdef VIZ

#include <igl/opengl/glfw/Viewer.h>
#include "Viz.h"

igl::opengl::glfw::Viewer viewer;
#endif


int main(int argc, char *argv[]){

    //Load mesh    
    //../mesh/terrain_yen
    //../mesh/hemispherical.obj
    //../mesh/terrain_rough.obj

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;    
    igl::readOBJ("../mesh/terrain_yen.obj", V, F);  


    Eigen::VectorXd vertex_func;
    compute_surface_func(V, F, vertex_func,HEIGHT);
    

#ifdef VIZ
    reset_viewer(viewer, V, F, vertex_func);
#endif


    int num_contours = 100;
    IsolinesDS iso_lines(V, F, vertex_func, num_contours);
        
    double tol = 0.0003*sqrt(dist(V.col(0).maxCoeff(), V.col(1).maxCoeff(),
        V.col(2).maxCoeff(), V.col(0).minCoeff(), V.col(1).minCoeff(),
        V.col(2).minCoeff()));

    std::vector<std::vector<Eigen::MatrixXd>> fermat = 
        fermat3d(&iso_lines, tol);

#ifdef VIZ
    viewer.launch();
#endif



    return 0;
}

