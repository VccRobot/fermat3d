#ifndef _VIZ_
#define _VIZ_

#include <igl/opengl/glfw/Viewer.h>
#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/jet.h>
#include <igl/colormap.h>

#include <fstream>



//********************** reset_viewer()
void reset_viewer(igl::opengl::glfw::Viewer&viewer,
    const Eigen::MatrixXd&V, const Eigen::MatrixXi&F, 
    const Eigen::VectorXd&vertex_func ,
    const bool is_colorize_vertex = true,
    const double line_width = 5.0) {

    //set the mesh by clearing any data that is currently present 
    viewer.data().clear();

    viewer.data().show_lines = false;//don't wanna show mesh edges    s
    viewer.data().set_face_based(true);
    viewer.data().line_width = line_width;
    viewer.data().set_mesh(V, F);
    viewer.core.align_camera_center(V, F);

    if (is_colorize_vertex) {
        //colorize the mesh vertices with the function 
        Eigen::MatrixXd func_colors;
        igl::jet(vertex_func, true, func_colors);
        viewer.data().set_colors(func_colors);
    }
   
}
//******************************************************************************


//********************** add_edge_to_viewer()
void add_edge_to_viewer(igl::opengl::glfw::Viewer&viewer,
	const Eigen::MatrixXd&iso_value,
	const std::vector<Eigen::MatrixXi>&iso_edges,
    const bool is_random_colors = true,
    const float r = 0.0f, const float g = 0.0f, const float b = 0.0f){
    
    
	Eigen::VectorXd val(iso_edges.size());
    Eigen::MatrixXd color;

    if (is_random_colors) {
        for (int i = 0; i < iso_edges.size(); i++) {
            Eigen::MatrixXd isoVV1(1, 3);
            igl::slice(iso_value, iso_edges.at(i).row(0).col(1),
                Eigen::RowVector3d(0, 1, 2), isoVV1);
            val(i) = isoVV1(1);
        }
        igl::colormap(igl::COLOR_MAP_TYPE_JET, val, true, color);
    }


	//pass the isolines edges 
	for (int i = 0; i < iso_edges.size(); i++) {

		const size_t num_edges = iso_edges.at(i).rows();

		Eigen::MatrixXd isoVV1(num_edges, 3), isoVV2(num_edges, 3);

		//source 
		igl::slice(iso_value, iso_edges.at(i).col(0),
			Eigen::RowVector3d(0, 1, 2), isoVV1);
		//target
		igl::slice(iso_value, iso_edges.at(i).col(1),
			Eigen::RowVector3d(0, 1, 2), isoVV2);

		viewer.data().add_edges(
			isoVV1,
			isoVV2,
			(is_random_colors) ? 
            Eigen::RowVector3d(color.row(i)) : Eigen::RowVector3d(r, g, b));
	}
}
//******************************************************************************

//********************** add_edge_to_viewer()
void add_colorful_edge_to_viewer(igl::opengl::glfw::Viewer&viewer,
    const Eigen::MatrixXd&iso_value,
    const std::vector<Eigen::MatrixXi>&iso_edges) {

    //colorful edges with the edge id in each isoline 
    
    for (int i = 0; i < iso_edges.size(); i++) {
        const size_t num_edges = iso_edges.at(i).rows();

        Eigen::VectorXd val(num_edges);
        
       // val(0) = 1.0f;//big number (1/0)

        for (int e =0 ; e < num_edges; e++) {
            val(e) = double(e);
        }

        Eigen::MatrixXd color;
        igl::colormap(igl::COLOR_MAP_TYPE_JET, val, true, color);


        for (int e = 0; e < num_edges; e++) {
            Eigen::MatrixXd isoVV1(1, 3), isoVV2(1, 3);
            int v0(iso_edges.at(i)(e, 0)), v1(iso_edges.at(i)(e, 1));

            isoVV1 = iso_value.row(v0);
            isoVV2 = iso_value.row(v1);

            
            viewer.data().add_edges(isoVV1, isoVV2,
                Eigen::RowVector3d(color.row(e)));
        }
    }
}
//******************************************************************************


//********************** add_edge_to_viewer()
void add_edge_to_viewer(igl::opengl::glfw::Viewer&viewer, 
	const std::vector<Eigen::MatrixXd>&points,
	bool is_connected = true) {
	//use this to draw the spiral 
	//each row in the points contains points that compose a contour 
	//listed in order 


	Eigen::MatrixXd contour_color = (Eigen::MatrixXd::Random(1, 3) +
		Eigen::MatrixXd::Constant(1, 3, 1.f)) / 2.0f;

	for (int i = 0; i < points.size(); i++) {

		for (int p = 0; p < points.at(i).rows() - 1; p++) {

			viewer.data().add_edges(points.at(i).row(p),
				points.at(i).row(p + 1), contour_color);
		}
		contour_color = (Eigen::MatrixXd::Random(1, 3) +
			Eigen::MatrixXd::Constant(1, 3, 1.f)) / 2.0f;
		if (is_connected) {
			//connect the last point in this conour witht the next one 
			if (i < points.size() - 1) {
				viewer.data().add_edges(points.at(i).row(points.at(i).rows() - 1),
					points.at(i + 1).row(0), contour_color);
			}
		}
	}

}
//******************************************************************************


//********************** add_vertices_id_text_to_viewer()
void add_vertices_id_text_to_viewer(igl::opengl::glfw::Viewer&viewer,
	const Eigen::MatrixXd&iso_value){

	for (int i = 0; i < iso_value.rows(); i++) {
		std::stringstream st;
		st << i;
		viewer.data().add_label(iso_value.row(i), st.str());
	}

}
//******************************************************************************


//**********************isolines_obj()
void isolines_obj(Eigen::MatrixXd&iso_value,
	std::vector<Eigen::MatrixXi>&iso_edges) {

	std::fstream file("isolines.obj", std::ios::out);
	file.precision(20);
	for (int i = 0; i < iso_value.rows(); i++) {
		file << "v " << iso_value(i, 0) << "  " << iso_value(i, 1) << "  "
			<< iso_value(i, 2) << std::endl;
	}

	for (int i = 0; i < iso_edges.size(); i++) {

		const size_t num_edges = iso_edges.at(i).rows();
		for (int j = 0; j < num_edges; j++) {
			int p0 = iso_edges.at(i)(j, 0);
			int p1 = iso_edges.at(i)(j, 1);
			file << "l " << p0 + 1 << " " << p1 + 1 << std::endl;

		}
	}
}
//******************************************************************************


//**********************spiral_obj()
void spiral_obj(std::vector<Eigen::MatrixXd>&points, bool is_connected = true)
{
	std::fstream file("spirals.obj", std::ios::out);
	file.precision(20);

	int total_edges(0);
	for (int i = 0; i < points.size(); i++) {
		for (int j = 0; j < points.at(i).rows() - 1; j++) {
			file << "v " << points.at(i)(j, 0) << " " << points.at(i)(j, 1)
				<< " " << points.at(i)(j, 2) << std::endl;
			file << "v " << points.at(i)(j + 1, 0) << " " 
				<< points.at(i)(j + 1, 1)
				<< " " << points.at(i)(j + 1, 2) << std::endl;
			total_edges++;
			file << "l " << total_edges * 2 - 1
				<< " " << total_edges * 2 << std::endl;
		}
		if (is_connected) {
			if (i < points.size() - 1) {
				file << "v " << points.at(i)(points.at(i).rows() - 1, 0)
					<< " " << points.at(i)(points.at(i).rows() - 1, 1)
					<< " " << points.at(i)(points.at(i).rows() - 1, 2)
					<< std::endl;

				file << "v " << points.at(i + 1)(0, 0)
					<< " " << points.at(i + 1)(0, 1)
					<< " " << points.at(i + 1)(0, 2)
					<< std::endl;

				total_edges++;
				file << "l " << total_edges * 2 - 1
					<< " " << total_edges * 2 << std::endl;

			}
		}
	}

}
//******************************************************************************

#endif /*_VIZ_*/
