#ifndef _COMPUTE_SURFACE_FUNC_
#define _COMPUTE_SURFACE_FUNC_

#include <igl/exact_geodesic.h>
enum FUNCTION {
    HEIGHT,
    GEODESIC
};

void compute_surface_func(const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F, Eigen::VectorXd & vertex_func,
    FUNCTION my_func) {

    //compute per vertex function and store it in vertex_func

    
    if (my_func == HEIGHT) {
        //using the height value as the per-vertex function
        vertex_func = V.col(1);
    }
    else if (my_func == GEODESIC) {
        int l = 0;
    }


    



}


#endif /*_COMPUTE_SURFACE_FUNC_*/
