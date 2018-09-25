#ifndef _ISOLINES_HELPER_
#define _ISOLINES_HELPER_

#include <random>

#include <igl/isolines.h>
#include <igl/project_to_line_segment.h>

#include "common.h"
#include "Viz.h"


extern igl::opengl::glfw::Viewer viewer;
extern Eigen::MatrixXd V;
extern Eigen::MatrixXi F;
extern Eigen::VectorXd vertex_func;


//********************** project_point_to_isoline()✔
void project_point_to_isoline(const Eigen::VectorXd & point,
    const Eigen::MatrixXd & isoV,
    const Eigen::MatrixXi & isoE,
    Eigen::VectorXd & t_min,
    Eigen::VectorXd & sqrD_min,
    int &edge_p0, int&edge_p1, int&edge_id) {

    //a variant of point projection to isolines where the isoline is passed by 
    //as a mesh (with line segments as mesh elements). the mesh elements are 
    //stored in isoE and the vertices can be indexed from isoV

    sqrD_min(0) = std::numeric_limits<double>::max();

    double min_dist = std::numeric_limits<double>::max();
    int pp;

    for (int i = 0; i < isoE.rows(); i++) {
        //for each edge in the isoline

        //the end points of the edge 
        int ed_p0(isoE(i, 0)), ed_p1(isoE(i, 1));

        //For some reason the 3rd parameter in project_to_line_segment
        //can be isoV.row(ed_p1) so we need another matrix for it

        Eigen::MatrixXd dest(1, 3), source(1,3), pt(1,3);
        dest(0) = isoV(ed_p1, 0);
        dest(1) = isoV(ed_p1, 1);
        dest(2) = isoV(ed_p1, 2);

        source(0) = isoV(ed_p0, 0);
        source(1) = isoV(ed_p0, 1);
        source(2) = isoV(ed_p0, 2);

        pt(0) = point(0);
        pt(1) = point(1);
        pt(2) = point(2);
                
        Eigen::VectorXd t(1), sqrD(1);

        igl::project_to_line_segment(pt, source, dest, t, sqrD);


        /*double dddd1 = dist(point(0), point(1), point(2), isoV(ed_p1, 0),
            isoV(ed_p1, 1), isoV(ed_p1, 2));
        double dddd0 = dist(point(0), point(1), point(2), isoV(ed_p0, 0),
            isoV(ed_p0, 1), isoV(ed_p0, 2));
        if (dddd1 < min_dist) { min_dist = dddd1; pp = ed_p1; }
        if (dddd0 < min_dist) { min_dist = dddd0; pp = ed_p0; }*/

        

        if (t(0) > -0.000001 && t(0) < 1.0 - 0.000001 &&
            sqrD_min(0) > sqrD(0)) {
            sqrD_min(0) = sqrD(0);
            t_min(0) = t(0);
            edge_p0 = ed_p0;
            edge_p1 = ed_p1;
            edge_id = i;
        }
    }

    /*std::cout << " pp=" << pp << "  edge_p0= " << edge_p0 << "  edge_p1= "
        << edge_p1 << std::endl;

    viewer.data().add_points(isoV.row(pp),
        Eigen::RowVector3d(1.0, 0.0, 1.0));
    viewer.data().add_points(isoV.row(edge_p0),
        Eigen::RowVector3d(0.0, 0.0, 0.0));
    viewer.data().add_points(isoV.row(edge_p1),
        Eigen::RowVector3d(0.0, 0.0, 0.0));*/
}
//******************************************************************************


//********************** project_point_to_isoline()✔
void project_point_to_isoline(const Eigen::VectorXd & point,
    const Eigen::MatrixXd & path,
    Eigen::VectorXd & t_min,
    Eigen::VectorXd & sqrD_min,
    int &edge_p0, int&edge_p1) {

    //variant of project_point_to_isoline that takes a list 
    //of points connected like a chain 
    //note that path does not form a closed chain 

    sqrD_min(0) = std::numeric_limits<double>::max();

    for (int i = 0; i < path.rows() - 1; i++) {
        //for each edge in the isoline

        //For some reason the 3rd parameter in project_to_line_segment
        //can be isoV.row(ed_p1) so we need another matrix for it

        Eigen::MatrixXd dest(1, 3), source(1, 3), pt(1, 3);
        dest(0) = path(i + 1, 0);
        dest(1) = path(i + 1, 1);
        dest(2) = path(i + 1, 2);

        source(0) = path(i, 0);
        source(1) = path(i, 1);
        source(2) = path(i, 2);

        pt(0) = point(0);
        pt(1) = point(1);
        pt(2) = point(2);

        Eigen::VectorXd t(1), sqrD(1);
        
        igl::project_to_line_segment(pt, source, dest, t, sqrD);

        if (t(0) > -0.000001 && t(0) < 1.0 - 0.000001 &&
            sqrD_min(0) > sqrD(0)) {
            sqrD_min(0) = sqrD(0);
            t_min(0) = t(0);
            edge_p0 = i;
            edge_p1 = i + 1;
        }
    }
}
//******************************************************************************


//********************** isoline_isoline_min_distance()✔
double isoline_isoline_min_distance(const Eigen::MatrixXi & source,
    const Eigen::MatrixXi & target,
    const Eigen::MatrixXd & isoV) {
    //compute the min distance from the source isoline to the target isoline
    //by picking NUM_SAMPLES samples from the source and find the min distance
    //to the target isoline (by projection)

    const int NUM_SAMPLES = 100;

    std::vector<int> indices(source.rows());
    std::iota(indices.begin(), indices.end(), 0);

    double min_dist = std::numeric_limits<double>::max();

    if (source.rows() > NUM_SAMPLES) {
        unsigned seed =
            std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(indices.begin(), indices.end(),
            std::default_random_engine(seed));
    }
    int sam_count = std::min(NUM_SAMPLES, int(source.rows()));

    for (int i = 0; i < sam_count; i++) {
        int pt_id = indices.at(i);

        Eigen::VectorXd  t(1), sqrD(1);
        int edge_p0, edge_p1, edge_id;
        //project this point onto the target
        Eigen::VectorXd pt(3);
        pt = isoV.row(source(pt_id));

        project_point_to_isoline(pt.transpose(),
            isoV, target, t, sqrD, edge_p0, edge_p1, edge_id);

        min_dist = std::min(min_dist, double(sqrD(0)));
    }

    return min_dist;

}
//******************************************************************************


//********************** sort_indexes()✔
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    //https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}
//******************************************************************************


//********************** remove_zero_length_isoE_points()✔
inline Eigen::MatrixXi remove_zero_length_isoE_points(
	const Eigen::MatrixXi & isoE) {

    //some edges has the same start and end points
    //we remove such edges here

	Eigen::MatrixXi isoE_filtered;
	isoE_filtered.resize(isoE.rows(), isoE.cols());
	int num_edges = 0;
	//remove edges with the same starting and end point
	for (int irow = 0; irow < isoE.rows(); irow++) {
		int ed_p0 = isoE(irow, 0);
		int ed_p1 = isoE(irow, 1);
		if (ed_p0 != ed_p1) {
			isoE_filtered.row(num_edges) = isoE.row(irow);
			num_edges++;
		}
	}
	//resize to fit 
	isoE_filtered.conservativeResize(num_edges, isoE_filtered.cols());
	return isoE_filtered;
}
//******************************************************************************


//********************** group_isolines()✔
std::vector<Eigen::MatrixXi> group_isolines(const Eigen::MatrixXd& isoV,
    const Eigen::MatrixXi& isoE,
    const Eigen::VectorXd& isoVal,
    const bool is_sort_ascend = false) {
    //reorder the isolines and return a vector such that each row in this 
    //vector is a single isoline with constant value 
    //set is_sort_ascend to true if you want the smallest value isoline to be 
    //the first in the returned vector 
    
    assert(isoV.rows() == isoVal.rows());

    //record  isoline values 
    //TODO use hash function instead of vector
    std::vector<double >isoVal_seen;
    for (int i = 0; i < isoV.rows(); i++) {
        bool seen = false;
        for (int j = 0; j < isoVal_seen.size(); j++) {
            if (abs(isoVal(i) - isoVal_seen[j]) < 0.0001) {
                seen = true;
                break;                
            }
        }
        if (!seen) {
            isoVal_seen.push_back(isoVal(i));
        }
    }
  

    if (is_sort_ascend) {
        std::sort(isoVal_seen.begin(), isoVal_seen.end());
    }
    else {
        std::sort(isoVal_seen.begin(), isoVal_seen.end(), std::greater<>());
    }

    //the return vector 
    std::vector < Eigen::MatrixXi> isoE_vec(isoVal_seen.size());


    for (int e = 0; e < isoE.rows(); e++) {
        //start and end point of the edge 
        int p1(isoE(e, 0)), p2(isoE(e, 1));

        //the function value on the edge start and end point 
        double v1(isoVal(p1)), v2(isoVal(p2));

        assert(abs(v1 - v2) < 2.2204e-10
            && "the two end points of an edge should have the same value");

        int iso_num = -1;
        //what the iso index 
        //TODO use hash function
        for (int s = 0; s < isoVal_seen.size(); s++) {
            if (abs(v1 - isoVal_seen.at(s)) < 2.2204e-10) {
                iso_num = s;
                break;
            }
        }

        //if we have not seen it before
        if (iso_num < 0) {
            /*iso_num = isoVal_seen.size();
            isoVal_seen.push_back(v1);

            Eigen::MatrixXi new_iso;
            isoE_vec.push_back(new_iso);*/
            std::cout<<"Error::group_isolines() can not find the "
                << "iso value index"<<std::endl;
            system("pause");
        }

        isoE_vec.at(iso_num).conservativeResize(
            isoE_vec.at(iso_num).rows() + 1, 2);
        isoE_vec.at(iso_num).row(isoE_vec.at(iso_num).rows() - 1) = 
            isoE.row(e);
        
    }

    //remove the empty isoE_vec
    for (int i = isoE_vec.size() - 1; i >= 0; i--) {
        if (isoE_vec.at(i).rows() == 0) {
            isoE_vec.erase(isoE_vec.begin() + i);                   
        }
        
    }

    return isoE_vec;
    
}
//******************************************************************************


//********************** create_chain()✔
void create_chain(std::vector<Eigen::MatrixXi>&seperate_isoE,
    const Eigen::MatrixXd& isoV) {
    //contour reorder by making sure that each contour is a single chain of
    //connected vertices (not just punch of seperate edges)
    //expected the contour to not have edges with same start and end point 
    //i.e., execute remove_zero_length_isoE_points()

    for (int cont = 0; cont < seperate_isoE.size(); cont++) {


        Eigen::MatrixXi new_cont;

       // std::cout << "\n\n \n " << seperate_isoE.at(cont) << std::endl;

        for (int ed_cur = 0; ed_cur < seperate_isoE.at(cont).rows() - 1;
            ed_cur++) {

            int ed_cur_p0 = seperate_isoE.at(cont)(ed_cur, 0);
            int ed_cur_p1 = seperate_isoE.at(cont)(ed_cur, 1);

            int ed_next_p0, ed_next_p1, ed_next;
            ed_next_p0 = ed_next_p1 = -1;

            for (int ed = ed_cur + 1; ed < seperate_isoE.at(cont).rows();
                ed++) {

                if (ed == ed_cur) { continue; }

                int ed_p0 = seperate_isoE.at(cont)(ed, 0);
                int ed_p1 = seperate_isoE.at(cont)(ed, 1);

                if ((ed_p1 == ed_cur_p1 || ed_p0 == ed_cur_p1)
                    || (ed_cur == 0 &&
                    (ed_p1 == ed_cur_p0 || ed_p0 == ed_cur_p0))) {
                    if ((ed_next_p0 >= 0 || ed_next_p1 >= 0) && ed_cur != 0) {
                        std::cout << "Error::create_chain() more than one" <<
                            "one edge to connect to" << std::endl;
                        system("pause");
                    }
                    ed_next_p0 = ed_p0;
                    ed_next_p1 = ed_p1;
                    ed_next = ed;
                }
            }

            if (ed_next_p0 < 0 || ed_next_p1 < 0) {
               // std::cout << "\n\n \n " << seperate_isoE.at(cont) << std::endl;
                std::cout << "Error::create_chain() can not find next edge"
                    << std::endl;
                system("pause");
            }


            seperate_isoE.at(cont).row(ed_next).swap(
                seperate_isoE.at(cont).row(ed_cur + 1));


            if (seperate_isoE.at(cont)(ed_cur, 1) !=
                seperate_isoE.at(cont)(ed_cur + 1, 0)) {
                //if the last point in the current edge does not match 
                //the first point in the first point in the next edge,
                //then we need to find the common point (same_p) and put it 
                //as first point in the next edge. 
                //if this is does not fix the point, then we can flip the 
                //current edge (mave the first point as last point) only if
                //this is the first edge 

                int same_p, diff_p;
                if (seperate_isoE.at(cont)(ed_cur, 0) ==
                    seperate_isoE.at(cont)(ed_cur + 1, 0) ||
                    seperate_isoE.at(cont)(ed_cur, 1) ==
                    seperate_isoE.at(cont)(ed_cur + 1, 0)) {
                    same_p = seperate_isoE.at(cont)(ed_cur + 1, 0);
                    diff_p = seperate_isoE.at(cont)(ed_cur + 1, 1);
                }
                else {
                    same_p = seperate_isoE.at(cont)(ed_cur + 1, 1);
                    diff_p = seperate_isoE.at(cont)(ed_cur + 1, 0);
                }


                seperate_isoE.at(cont)(ed_cur + 1, 0) = same_p;
                seperate_isoE.at(cont)(ed_cur + 1, 1) = diff_p;

                if (seperate_isoE.at(cont)(ed_cur, 1) !=
                    seperate_isoE.at(cont)(ed_cur + 1, 0)) {
                    if (ed_cur == 0) {
                        //it is okay to flip the very first edge 
                        std::swap(seperate_isoE.at(cont)(ed_cur, 0),
                            seperate_isoE.at(cont)(ed_cur, 1));
                    }
                    else {
                        std::cout << "Error::create_chain() can't find"
                            " good order !!!"
                            << std::endl;
                        system("pause");
                    }
                }
            }
        }
    }
}
//******************************************************************************


//********************** create_vertex_edge()✔
void create_vertex_edge(const Eigen::MatrixXi isoE,
    Eigen::MatrixXi & vertex_edge) {

    //Construct vertex-edge matrix 
    vertex_edge = Eigen::MatrixXi::Constant(vertex_edge.rows(), 2, -1);
    
    for (int e = 0; e < isoE.rows(); e++) {
        int v0 = isoE(e, 0);
        int v1 = isoE(e, 1);

        if (vertex_edge(v0, 0) == -1) { vertex_edge(v0, 0) = e; }
        else if (vertex_edge(v0, 1) == -1) { vertex_edge(v0, 1) = e; }
        else {
            std::cout << "Error::multi_spirallable_region() " <<
                "vertex [" << v0 << "] connected to more than two edges"
                << std::endl;
            system("pause");
        }

        if (vertex_edge(v1, 0) == -1) { vertex_edge(v1, 0) = e; }
        else if (vertex_edge(v1, 1) == -1) { vertex_edge(v1, 1) = e; }
        else {
            std::cout << "Error::multi_spirallable_region() " <<
                "vertex [" << v1 << "] connected to more than two edges"
                << std::endl;
            system("pause");
        }
    }   

}
//******************************************************************************


//********************** create_chain_use_vertex_edge()✔
void create_chain_use_vertex_edge(const Eigen::MatrixXi & iso_edges,
    const  Eigen::MatrixXi & vertex_edge,
    const Eigen::MatrixXd & isoV,
    std::vector<Eigen::MatrixXi> & chain_iso_edges) {
    //A varaint of create_chain where we use the vertex_edge data structure 
    //iso_edges is a collection of edges whose vertices have the same iso value
    //these edges may belong to one or more spirallable regions 
    //vertex_edge gives the two edges a vertex belongs to 
    //chain_iso_edges is the returned vector where each row is a seperate 
    //chain (belong to different spirallable region) and the index represent
    //edge in iso_edges 


    int starting_edge(0), starting_v(iso_edges(0, 0));
    //Edge id that is being processed
    int my_edge = starting_edge;

    //we stop when all edges are taken
    int num_edge_taken = 0;

    std::vector<bool> taken(iso_edges.rows(), false);
    std::vector<bool> taken_v(vertex_edge.rows(), false);

    while (num_edge_taken < iso_edges.rows()) {

        //The current chain we are building 
        Eigen::MatrixXi my_chain(iso_edges.rows()- num_edge_taken +1, 2);
        int num_edges_chain = 0;
        while (true) {         
          
            //the start and end of the edge 
            int v0 = iso_edges(my_edge, 0);
            int v1 = iso_edges(my_edge, 1);

            //add the edge (its vertices) to the chain
            my_chain.row(num_edges_chain) = iso_edges.row(my_edge);
            if (num_edges_chain > 1) {
                //make sure that the start point of this edge is the end 
                //point of the previous edge 
                if (my_chain(num_edges_chain, 0) !=
                    my_chain(num_edges_chain - 1, 1)) {
                    std::swap(my_chain(num_edges_chain, 0),
                        my_chain(num_edges_chain, 1));
                }
            }
            num_edges_chain++;
            num_edge_taken++;
                        
            //mark the edge as being taken
            taken.at(my_edge) = true;

#if 1
            if (vertex_edge(v1, 0) != my_edge &&
                vertex_edge(v1, 1) != my_edge) {
                std::cout << "Error::create_chain_use_vertex_edge() "
                    << "the end vertex of an edge does not contain the edge\n"
                    << " v0= " << v0 << " v1= " << v1 << " my_edge= " 
                    << my_edge << std::endl;
                system("pause");
            }
#endif 
            //move to the next edge 
            //my_edge = (vertex_edge(v1, 0) != my_edge) ? vertex_edge(v1, 0) :
            //    vertex_edge(v1, 1);

            if ((v1 == starting_v || v0 == starting_v) &&
                my_edge != starting_edge) {
                taken_v[v0] = taken_v[v1] = true;
                break;
            }
           

            int my_v = v1;
            if (taken_v[my_v]) { my_v = v0; }
            if (taken_v[my_v] && v1!= starting_v && v0!= starting_v) {                
                std::cout<<"Error::create_chain_use_vertex_edge()"
                    <<" both vertices of an edge are added!!!"<<std::endl;
                system("pause");
            }
                       
            my_edge = (vertex_edge(my_v, 0) != my_edge) ?
                vertex_edge(my_v, 0) : vertex_edge(my_v, 1);

            //mark the taken vertices 
            taken_v[v0] = taken_v[v1] = true;

            if (my_edge == starting_edge) { break; }
        }

        /*
        //Make sure that the edges are ordered in the same direction as other 
        //chains. First compute the center mass of the chain and then check
        //the angle of the end points of the first edge 
        Eigen::RowVector3d chain_center(3,0);
        for (int v = 0; v < my_chain.rows(); v++) {
        //    chain_center += my_chain.row(v);
        }
        chain_center /= my_chain.rows();

        //get angle between st0 - center - end0
        //where st0 and end0 are the start and end points of the first edge in
        //the chain
        double s1 = dist(isoV(my_chain(0, 1), 0),
            isoV(my_chain(0, 1), 1), isoV(my_chain(0, 1), 2),
            isoV(my_chain(0, 0), 0),
            isoV(my_chain(0, 0), 1), isoV(my_chain(0, 0), 2));

        double s2 = dist(chain_center(0), chain_center(1),
            chain_center(2), isoV(my_chain(0, 1), 0),
            isoV(my_chain(0, 1), 1), isoV(my_chain(0, 1), 2));

        double s3 = dist(chain_center(0), chain_center(1),
            chain_center(2), isoV(my_chain(0, 0), 0),
            isoV(my_chain(0, 0), 1), isoV(my_chain(0, 0), 2));

        
        double angle = acos((s3 + s2 - s1) / (2.*sqrt(s3*s2)));
        */

        //add the newly-created chain only if it is big enough 
        if (num_edges_chain > MIN_NUM_POINT_PER_CONTOUR) {
            chain_iso_edges.push_back(my_chain.topRows(num_edges_chain));

        }

        //if we have processed all the edges 
        if (num_edge_taken >= iso_edges.rows()) { break; }

        //look for the next not-taken edge to be our starting edge 
        starting_edge = -1;
        for (int i = 0; i < taken.size(); i++) {
            if (!taken.at(i)) {
                starting_edge = i;
                break;
            }
        }

        if (starting_edge < 0) { break; }
        my_edge = starting_edge;
        starting_v = iso_edges(starting_edge, 0);    
     
    }
}
//******************************************************************************


//********************** isolate_isolines()✔
void isolate_isolines(const std::vector<Eigen::MatrixXi>&all_seperate_isoE,
    const Eigen::MatrixXd & isoV,
    std::vector <std::vector<Eigen::MatrixXi>>&isolate_isoE) {

    //Isolate the isolines into spirallable regions. Each spirallable region 
    //is set of isolines such that each isolines (its vertices) have a unique 
    //values (different than another isoline in the same region)
    //Input: in each row in all_seperate_isoE, we store the set of edges 
    //that has the same iso value (they are not ordered or sorted in any way)

    //Output: isolate_isoE[X][Y] where X is index for the spirallable region 
    //Y is the contour index (betwee 0 and num_contours -1). isolate_isoE[X][Y]
    //gives the edges as Eigen::MatrixXi of contour Y in the spiralle region X

    //TODO filter out isolines with few points (less than 
    //MIN_NUM_POINT_PER_CONTOUR)

    //We first sort the isolines such that each isolines form a connected 
    //closed chain. Now we have one or more chain for the current processed 
    //iso value. We then add the chain to the correct spirallable region 


    //because the first one is the outer boundary edges 
   // assert(!isolate_isoE.at(0).empty());
    int region_id = 1;
    int num_regions = 1;
    //vertex_edge data structure is build for each isoline (isoline whose 
    //vertices has the same iso value) such that each vertex points to the 
    //two edges it belongs to
    Eigen::MatrixXi vertex_edge(isoV.rows(), 2);

    //We keep track of the regions that we have last added new isolines to
    //because these are the regions that will grow further.
    std::vector<int> last_region(1, 0);


    for (int i = 0; i < all_seperate_isoE.size(); i++) {

        std::vector<Eigen::MatrixXi> chains;

        create_vertex_edge(all_seperate_isoE.at(i), vertex_edge);

        create_chain_use_vertex_edge(all_seperate_isoE.at(i), vertex_edge,
            isoV, chains);

        if (false) {
            reset_viewer(viewer, V, F, vertex_func);
            add_edge_to_viewer(viewer, isoV, chains);
            viewer.launch();
        }

        if (chains.size() == 1 &&
            (isolate_isoE.size() == 0 || isolate_isoE.size() == 1)) {
            //Special case
           //The first couple of isolines may not need checks. This is the case
            //of haveing one region and the chains are just one so it is always
            //added to the exisiting region (creating one spirallalbe region)
            if (isolate_isoE.size() == 0) {
                std::vector<Eigen::MatrixXi> ch;
                ch.push_back(chains.at(0));
                isolate_isoE.push_back(ch);
            }
            else {
                isolate_isoE.at(0).push_back(chains.at(0));
            }
            
        }
        else {
            //General case  

            //For each chain, we store its distance to all the region that we 
            //have updated last time along with
            //the min distance as pair (which is used later for sorting
            //to figure out the closest region / containing region of the chain)
            std::vector<std::vector<std::pair<double, int>>>
                chain_dist_to_region;
            chain_dist_to_region.reserve(chains.size());

            for (int c = 0; c < chains.size(); c++) {
                std::vector < std::pair<double, int>> my_chain_distances;

                //find the min distance to all the regions that has been 
                //updated or created last time 
                for (int r = 0; r < last_region.size(); r++) {
                    int reg = last_region.at(r);
                    double min_dist = isoline_isoline_min_distance(
                        chains.at(c), isolate_isoE.at(reg).back(), isoV);

                    my_chain_distances.push_back(
                        std::make_pair(min_dist, reg));

                }
                chain_dist_to_region.push_back(my_chain_distances);
            }


            std::vector<std::vector<int>> add_to_region(isolate_isoE.size());
            //Sort all the collected distances for each chain 
            //from there we know which region the chain should be added to 
            for (int c = 0; c < chains.size(); c++) {
                //sort key-value pairs 
                std::sort(chain_dist_to_region.at(c).begin(),
                    chain_dist_to_region.at(c).end());

                int closest_region = chain_dist_to_region.at(c).front().second;

                //add that we are going to add this chain to that region
                add_to_region.at(closest_region).push_back(c);
            }

            std::vector<int> current_region;

            //actually add the chain to the region only if we are going to
            //add only one chain to this region
            for (int r = 0; r < isolate_isoE.size(); r++) {
                if (add_to_region.at(r).size() == 1) {
                    int ch = add_to_region.at(r).front();
                    isolate_isoE.at(r).push_back(chains.at(ch));

                    current_region.push_back(r);
                }
            }

            //if we found that we need to add more than one chain to the region
            //then we don't add them to these regions but instead open new 
            //region using them 
            for (int r = isolate_isoE.size() - 1; r >= 0; r--) {
                if (add_to_region.at(r).size() > 1) {
                    for (int c = 0; c < add_to_region.at(r).size(); c++) {

                        int ch = add_to_region.at(r).at(c);

                        std::vector<Eigen::MatrixXi> mat_V;
                        mat_V.push_back(chains.at(ch));

                        isolate_isoE.push_back(mat_V);

                        current_region.push_back(isolate_isoE.size() - 1);
                    }
                }
            }

            last_region.swap(current_region);

        }

        region_id++;
    }

    for (int i = 0; i < isolate_isoE.size(); i++) {
        std::cout << "i= " << i << " --> " << isolate_isoE.at(i).size() << std::endl;
    }


    if (false) {
        reset_viewer(viewer, V, F, vertex_func);
        for (int c = 0; c < isolate_isoE.size(); c++) {
            if (false) {
                //draw each region in a different color
                float r, g, b(0.0f);
                if (isolate_isoE.at(c).size() > 1) {
                    r = float(rand()) / float(RAND_MAX);
                    g = float(rand()) / float(RAND_MAX);
                    b = float(rand()) / float(RAND_MAX);
                }
                else {
                    r = 0.0f;
                    g = 0.0f;
                    b = 0.0f;
                }
                add_edge_to_viewer(viewer, isoV, isolate_isoE.at(c), false,
                    r, g, b);
            }
            else {
                //draw each isoline in raincolor relevent to the increasing 
                //order of its edges 
                add_colorful_edge_to_viewer(viewer, isoV, isolate_isoE.at(c));
            }
        }
        viewer.launch();
    }  
  
}
//******************************************************************************
#endif /*_ISOLINES_HELPER_*/