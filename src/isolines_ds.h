#ifndef _ISOLINES_DS_
#define _ISOLINES_DS_

#include <chrono>
#include <random>
#include <numeric>

#include <igl/isolines.h>
#include <igl/project_to_line_segment.h>

#include "common.h"


#ifdef VIZ
#include "Viz.h"
extern igl::opengl::glfw::Viewer viewer;
#endif 



class IsolinesDS
{
public:
    IsolinesDS() {};

    //Populate the basic data structure that is used for further processing 
    //of the isolines (isoV and isoE)
    IsolinesDS(Eigen::MatrixXd&V, Eigen::MatrixXi&F,
        Eigen::VectorXd& vertex_func, const int num_contours);

    //Find the projection of a point onto a single isoline (my_isoline)
    //whose the indeices of its edges can index isoV
    void project_point_to_isoline(const Eigen::VectorXd & point,
        const Eigen::MatrixXi & my_isoline,
        Eigen::VectorXd & t_min,
        Eigen::VectorXd & sqrD_min,
        int &edge_p0, int&edge_p1, int&edge_id);

    //Find the projection of a point onto a single path defined by a set of 
    //points (given their xy coordinates)
    void project_point_to_isoline(const Eigen::VectorXd & point,
        const Eigen::MatrixXd & path,
        Eigen::VectorXd & t_min,
        Eigen::VectorXd & sqrD_min,
        int &edge_p0, int&edge_p1);

    //*** getter
    inline int get_num_spiral_region() {
        //return the number of spirallal regions 
        return spiral_regions.size();
    }

    inline int get_num_contours(int r) {
        //return the number of contours in region number r 
        return spiral_regions.at(r).size();
    }

    inline std::vector<Eigen::MatrixXi> get_region(int id) {
        //return the region by its index
        return spiral_regions.at(id);
    }

    inline Eigen::MatrixXd get_coordinates(int v_id) {
        //return the vertex cooridnates. vertex here is an isoline vertex 
        //not a mesh vertex
        return isoV.row(v_id);
    }


    ~IsolinesDS() {};

private:

    //******** VARIABLES *********//

    Eigen::MatrixXd V;//Mesh input vertices 
    Eigen::MatrixXi F;//Mesh input faces
    Eigen::VectorXd vertex_func; //Function defined over the mesh 

    Eigen::MatrixXd isoV;//isolines vertices 
    Eigen::MatrixXi isoE;//isolines edges 

    int num_contours; //number of contours used to get the isolines 

    //groups the isolines by their isovalue
    std::vector<Eigen::MatrixXi> isolines_groups; 
    
    //groups of spirallable region
    //spiral_region[X][Y].row(Z) gives the end points of the Zth edges 
    //of the Yth isolines in the Xth spirallable region
    std::vector<std::vector<Eigen::MatrixXi>> spiral_regions;



    //******** FUNCTIONS *********//

    //Get the isovalue of the vertex
    //vertex here is a vertex on the isolines not a mesh vertex 
    auto isovalue(int vertex);

    //Populates isolines_groups
    void grouping_isolines(const bool is_sort_ascend = false);

    //Clean up the iso edge from edges that has the same start and end point 
    void remove_zero_length_isoE_points();
    

    //Find the spirallable region (grouping_isolines() should be called first)
    void spirallable_region();

    //Find the one-sided Hausdorff distance between two isolines (source and 
    //target) by picking NUM_SAMPLES of the source vertices and find their 
    //projection on the target 
    double isoline_isoline_min_distance(const Eigen::MatrixXi & source,
        const Eigen::MatrixXi & target, const int NUM_SAMPLES = 100);


    //Creating vertex edge data structure which gives for each vertex 
    //what edges it is connected to. Vertex here is a vertex in the isolines
    //not a mesh vertex
    void create_vertex_edge(const Eigen::MatrixXi isoE,
        Eigen::MatrixXi & vertex_edge);

    //Building chains from edges (their vertices) have the same isovalue
    void create_chain_use_vertex_edge(const Eigen::MatrixXi & iso_edges,
        const  Eigen::MatrixXi & vertex_edge,
        std::vector<Eigen::MatrixXi> & chain_iso_edges);

};

//********************** IsolinesDS::IsolinesDS()
IsolinesDS::IsolinesDS(Eigen::MatrixXd&V, Eigen::MatrixXi&F,
    Eigen::VectorXd& vertex_func,  const int num_contours) :
    V(V),F(F),vertex_func(vertex_func), num_contours(num_contours) {

    //Construct the isolines 
    igl::isolines(V, F, vertex_func, num_contours, isoV, isoE);

    //clean up by removing the duplicated edges (edges with same start and end
    //point)
    remove_zero_length_isoE_points();

    //group the isolines by their isovalue 
    grouping_isolines(true);

    //Creating spirallable regions 
    spirallable_region();

    //for (int i = 0; i < spiral_regions.size(); i++) {
    //    add_edge_to_viewer(viewer, isoV, spiral_regions.at(i), false,
    //        double(rand()) / double(RAND_MAX),
    //        double(rand()) / double(RAND_MAX), 
    //        double(rand()) / double(RAND_MAX));
    //}
    //viewer.launch();

}
//*****************************************************************************


//********************** IsolinesDS::isovalue()
auto IsolinesDS::isovalue(int vertex) {
    return isoV(vertex, 1);
}
//*****************************************************************************

//********************** IsolinesDS::remove_zero_length_isoE_points()
void IsolinesDS::remove_zero_length_isoE_points() {

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
    isoE.swap(isoE_filtered);
}
//******************************************************************************


//********************** IsolinesDS::seperate_isolines()
void IsolinesDS::grouping_isolines(const bool is_sort_ascend) {
    //populates seperate_isolines such that each row in this 
    //vector is a single isoline with the same isovalue 
    //set is_sort_ascend to true if you want the smallest value isoline to be 
    //the first in the returned vector 

    
    //record  isoline values 
    //TODO use hash function instead of vector
    std::vector<double >isoVal_seen;
    for (int i = 0; i < isoV.rows(); i++) {
        bool seen = false;
        for (int j = 0; j < isoVal_seen.size(); j++) {
            if (abs(isovalue(i) - isoVal_seen[j]) < 0.0001) {
                seen = true;
                break;
            }
        }
        if (!seen) {
            isoVal_seen.push_back(isovalue(i));
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
        double v1(isovalue(p1)), v2(isovalue(p2));

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
            std::cout << "Error::group_isolines() can not find the "
                << "iso value index" << std::endl;
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

    isolines_groups.clear();
    isoE_vec.swap(isolines_groups);

}
//*****************************************************************************


//********************** IsolinesDS::isoline_isoline_min_distance()
double IsolinesDS::isoline_isoline_min_distance(const Eigen::MatrixXi & source,
    const Eigen::MatrixXi & target, const int NUM_SAMPLES) {
    //compute the min distance from the source isoline to the target isoline
    //by picking NUM_SAMPLES samples from the source and find the min distance
    //to the target isoline (by projection)

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
            target, t, sqrD, edge_p0, edge_p1, edge_id);

        min_dist = std::min(min_dist, double(sqrD(0)));
    }

    return min_dist;
}
//*****************************************************************************


//********************** IsolinesDS::project_point_to_isoline()
void IsolinesDS::project_point_to_isoline(const Eigen::VectorXd & point, 
    const Eigen::MatrixXi & my_isoline,
    Eigen::VectorXd & t_min,
    Eigen::VectorXd & sqrD_min,
    int &edge_p0, int&edge_p1, int&edge_id) {

    //a variant of point projection to isolines where the isoline is passed by 
    //as a mesh (with line segments as mesh elements). the mesh elements are 
    //stored in isoE and the vertices can be indexed from isoV

    sqrD_min(0) = std::numeric_limits<double>::max();

    for (int i = 0; i < my_isoline.rows(); i++) {
        //for each edge in the isoline

        //the end points of the edge 
        int ed_p0(my_isoline(i, 0)), ed_p1(my_isoline(i, 1));

        //For some reason the 3rd parameter in project_to_line_segment
        //can be isoV.row(ed_p1) so we need another matrix for it

        Eigen::MatrixXd dest(1, 3), source(1, 3), pt(1, 3);
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

        if (t(0) > -0.000001 && t(0) < 1.0 - 0.000001 &&
            sqrD_min(0) > sqrD(0)) {
            sqrD_min(0) = sqrD(0);
            t_min(0) = t(0);
            edge_p0 = ed_p0;
            edge_p1 = ed_p1;
            edge_id = i;
        }
    }
}
//*****************************************************************************


//********************** IsolinesDS::project_point_to_isoline()
void IsolinesDS::project_point_to_isoline(const Eigen::VectorXd & point,
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
//*****************************************************************************


//********************** create_vertex_edge()
void IsolinesDS::create_vertex_edge(const Eigen::MatrixXi isoE,
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
//*****************************************************************************


//********************** IsolinesDS::create_chain_use_vertex_edge()
void IsolinesDS::create_chain_use_vertex_edge(
    const Eigen::MatrixXi & iso_edges,
    const  Eigen::MatrixXi & vertex_edge,
    std::vector<Eigen::MatrixXi> & chain_iso_edges) {

    //contour reorder by making sure that each contour is a single chain of
    //connected vertices (not just punch of seperate edges)
    //expected the contour to not have edges with same start and end point 
    //i.e., execute remove_zero_length_isoE_points()

    //we use the vertex_edge data structure 
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
        Eigen::MatrixXi my_chain(iso_edges.rows() - num_edge_taken + 1, 2);
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
            if (taken_v[my_v] && v1 != starting_v && v0 != starting_v) {
                std::cout << "Error::create_chain_use_vertex_edge()"
                    << " both vertices of an edge are added!!!" << std::endl;
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
//*****************************************************************************

//********************** spirallable_region()
void IsolinesDS::spirallable_region() {

    //Isolate the isolines into spirallable regions. Each spirallable region 
    //is set of isolines such that each isolines (its vertices) have a unique 
    //values (different than another isoline in the same region)
    //Input: in each row in isolines_groups, we store the set of edges 
    //that has the same iso value (they are not ordered or sorted in any way)

    //Output: spiral_regions[X][Y] where X is index for the spirallable region 
    //Y is the contour index (betwee 0 and num_contours -1).
    //spiral_regions[X][Y] gives the edges as Eigen::MatrixXi of 
    //contour Y in the spiralle region X

    //TODO filter out isolines with few points (less than 
    //MIN_NUM_POINT_PER_CONTOUR)

    //We first sort the isolines such that each isolines form a connected 
    //closed chain. Now we have one or more chain for the current processed 
    //iso value. We then add the chain to the correct spirallable region 


    //because the first one is the outer boundary edges 
   // assert(!spiral_regions.at(0).empty());
    int region_id = 1;
    int num_regions = 1;
    //vertex_edge data structure is build for each isoline (isoline whose 
    //vertices has the same iso value) such that each vertex points to the 
    //two edges it belongs to
    Eigen::MatrixXi vertex_edge(isoV.rows(), 2);

    //We keep track of the regions that we have last added new isolines to
    //because these are the regions that will grow further.
    std::vector<int> last_region(1, 0);


    for (int i = 0; i < isolines_groups.size(); i++) {

        std::vector<Eigen::MatrixXi> chains;

        create_vertex_edge(isolines_groups.at(i), vertex_edge);

        create_chain_use_vertex_edge(isolines_groups.at(i), vertex_edge, 
            chains);

        if (chains.size() == 1 &&
            (spiral_regions.size() == 0 || spiral_regions.size() == 1)) {
            //Special case
           //The first couple of isolines may not need checks. This is the case
            //of haveing one region and the chains are just one so it is always
            //added to the exisiting region (creating one spirallalbe region)
            if (spiral_regions.size() == 0) {
                std::vector<Eigen::MatrixXi> ch;
                ch.push_back(chains.at(0));
                spiral_regions.push_back(ch);
            }
            else {
                spiral_regions.at(0).push_back(chains.at(0));
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
                        chains.at(c), spiral_regions.at(reg).back());

                    my_chain_distances.push_back(
                        std::make_pair(min_dist, reg));

                }
                chain_dist_to_region.push_back(my_chain_distances);
            }


            std::vector<std::vector<int>> add_to_region(spiral_regions.size());
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
            for (int r = 0; r < spiral_regions.size(); r++) {
                if (add_to_region.at(r).size() == 1) {
                    int ch = add_to_region.at(r).front();
                    spiral_regions.at(r).push_back(chains.at(ch));

                    current_region.push_back(r);
                }
            }

            //if we found that we need to add more than one chain to the region
            //then we don't add them to these regions but instead open new 
            //region using them 
            for (int r = spiral_regions.size() - 1; r >= 0; r--) {
                if (add_to_region.at(r).size() > 1) {
                    for (int c = 0; c < add_to_region.at(r).size(); c++) {

                        int ch = add_to_region.at(r).at(c);

                        std::vector<Eigen::MatrixXi> mat_V;
                        mat_V.push_back(chains.at(ch));

                        spiral_regions.push_back(mat_V);

                        current_region.push_back(spiral_regions.size() - 1);
                    }
                }
            }

            last_region.swap(current_region);

        }

        region_id++;
    }

    for (int i = 0; i < spiral_regions.size(); i++) {
        std::cout << "i= " << i << " --> " <<
            spiral_regions.at(i).size() << std::endl;
    }

#ifdef VIZ
    if (false) {
        reset_viewer(viewer, V, F, vertex_func);
        for (int c = 0; c < spiral_regions.size(); c++) {
            if (false) {
                //draw each region in a different color
                float r, g, b(0.0f);
                if (spiral_regions.at(c).size() > 1) {
                    r = float(rand()) / float(RAND_MAX);
                    g = float(rand()) / float(RAND_MAX);
                    b = float(rand()) / float(RAND_MAX);
                }
                else {
                    r = 0.0f;
                    g = 0.0f;
                    b = 0.0f;
                }
                add_edge_to_viewer(viewer, isoV, spiral_regions.at(c), false,
                    r, g, b);
            }
            else {
                //draw each isoline in raincolor relevent to the increasing 
                //order of its edges 
                add_colorful_edge_to_viewer(viewer, isoV, 
                    spiral_regions.at(c));
            }
        }
        viewer.launch();
    }
#endif

}
//*****************************************************************************

#endif /*_ISOLINES_DS_*/

