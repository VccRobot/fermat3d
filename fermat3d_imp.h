#ifndef _FERMAT3D_IMP_
#define _FERMAT3D_IMP_

#include "common.h"
#include "isolines_ds.h"


#ifdef VIZ
#include "Viz.h"
extern igl::opengl::glfw::Viewer viewer;
#endif 


//********************** get_spiral3d_one_region()
std::vector<Eigen::MatrixXd> get_spiral3d_one_region(IsolinesDS * iso_lines,
    const int reg_id,   
    const double tol) {

    const std::vector<Eigen::MatrixXi> region_isoE = iso_lines->get_region(reg_id);

    //Generate a single spiral from a (spirallable) isolines in region_isoE
    //each row in region_isoE stores the points of single contours 

    //Starting with the outer contour (indesx 0) we keep rerouting to form 
    //the spiral

    std::vector<Eigen::MatrixXd> spiral;

    int str_edge_id, end_edge_id;

    Eigen::VectorXd p_end;

    for (int i = 0; i < region_isoE.size(); i++) {
        const int num_edges = region_isoE.at(i).rows();

        Eigen::MatrixXd contour(num_edges, 3);
        
        int direction;
        int p_str_id;
       

        if (i == 0) {
            //for the first contour, we can pick any vertex to be our starting
            //point 
            str_edge_id = 0;
            p_str_id = region_isoE.at(i)(str_edge_id, 0);
            
            direction = 1;
            end_edge_id = num_edges - 1;
        }
        else {
            //For subsequent contours, we have to choose the right direction 
            //so that the spiral can continue growing 
            //Also we pick the next edge (not the edge that we projected 
            //on it from previous isoline)    
           
            //1) using distance 
            int str_edge_id_p0, str_edge_id_p1;
            str_edge_id_p0 = region_isoE.at(i)(str_edge_id, 0);
            str_edge_id_p1 = region_isoE.at(i)(str_edge_id, 1);

            /*double dd0 = dist(isoV(p_str_id, 0),
                isoV(p_str_id, 1),
                isoV(p_str_id, 2),
                isoV(str_edge_id_p0, 0),
                isoV(str_edge_id_p0, 1),
                isoV(str_edge_id_p0, 2));

            double dd1 = dist(isoV(p_str_id, 0),
                isoV(p_str_id, 1),
                isoV(p_str_id, 2),
                isoV(str_edge_id_p1, 0),
                isoV(str_edge_id_p1, 1),
                isoV(str_edge_id_p1, 2));*/

            double dd0 = dist(
                iso_lines->get_coordinates(p_str_id)(0),
                iso_lines->get_coordinates(p_str_id)(1),
                iso_lines->get_coordinates(p_str_id)(2),
                iso_lines->get_coordinates(str_edge_id_p0)(0),
                iso_lines->get_coordinates(str_edge_id_p0)(1),
                iso_lines->get_coordinates(str_edge_id_p0)(2));

            double dd1 = dist(
                iso_lines->get_coordinates(p_str_id)(0),
                iso_lines->get_coordinates(p_str_id)(1),
                iso_lines->get_coordinates(p_str_id)(2),
                iso_lines->get_coordinates(str_edge_id_p1)(0),
                iso_lines->get_coordinates(str_edge_id_p1)(1),
                iso_lines->get_coordinates(str_edge_id_p1)(2));

                  end_edge_id = str_edge_id;


            if (dd1 < dd0) {
                //set direction 
                p_str_id = str_edge_id_p1;
                direction = 1;

                //get next edge 
                str_edge_id = (str_edge_id == region_isoE.at(i).rows() - 1)
                    ? 0 : str_edge_id + 1;
            }
            else {
                //set direction
                p_str_id = str_edge_id_p0;
                direction = -1;

                //get next edge 
                str_edge_id = (str_edge_id == 0) ?
                    region_isoE.at(i).rows() - 1 : str_edge_id - 1;
            }

           /*
           //2) using vectors 
           Eigen::VectorXd p_proj = spiral.back().row(
                spiral.back().rows() - 1);

            Eigen::VectorXd va0 = isoV.row(str_edge_id_p0) - p_proj.transpose();
            Eigen::VectorXd va1 = isoV.row(str_edge_id_p1) - p_proj.transpose();
            Eigen::VectorXd vb;
            if (spiral.size() == 1) {
                vb = isoV.row(p_str_id) - p_end.transpose();
            }
            else {
                Eigen::MatrixXd last_pt_prv_contour =
                    spiral.at(spiral.size() - 2).row(
                        spiral.at(spiral.size() - 2).rows() - 1);
                vb = isoV.row(p_str_id) - last_pt_prv_contour;
            }

            va0.normalize();
            va1.normalize();
            vb.normalize();
                        
            double dir0, dir1;
            dir0 = va0.dot(vb);
            dir1 = va1.dot(vb);

            //***
            if (false) {
                add_edge_to_viewer(viewer, spiral,false);

                viewer.data().add_points(isoV.row(p_str_id),
                    Eigen::RowVector3d(0.0, 0, 0.0));
               
                viewer.data().add_points(p_proj.transpose(),
                    Eigen::RowVector3d(1.0, 0, 1.0));
                viewer.data().add_points(isoV.row(str_edge_id_p0),
                    Eigen::RowVector3d(1.0, 0, 0.0));
                viewer.data().add_points(isoV.row(str_edge_id_p1),
                    Eigen::RowVector3d(0.0, 0, 1.0));
                
                std::cout << isoV.row(str_edge_id_p0) << std::endl;
                std::cout << isoV.row(str_edge_id_p1) << std::endl;
                std::cout << p_proj.transpose() << std::endl;

                if (i == 1) {
                    viewer.data().add_points(
                        p_end.transpose(),
                        Eigen::RowVector3d(0.0, 1.0, 0.0));
                }
                else {
                    viewer.data().add_points(
                        spiral.at(spiral.size() - 2).row(
                            spiral.at(spiral.size() - 2).rows() - 1),
                        Eigen::RowVector3d(0.0, 1.0, 0.0));
                }


                viewer.launch();
            }


            if (dir0 *dir1 > 0.0000001) {
                std::cout << "Error::get_spiral3d() can not decide which"
                    << "direction to take!!!" << std::endl;
                system("pause");
            }
            
            end_edge_id = str_edge_id;

            if (dir1 > 0.0) {
                //set direction 
                p_str_id = str_edge_id_p1;
                direction = 1;

                //get next edge 
                str_edge_id = (str_edge_id == region_isoE.at(i).rows() - 1)
                    ? 0 : str_edge_id + 1;
            }
            else {
                //set direction
                p_str_id = str_edge_id_p0;
                direction = -1;

                //get next edge 
                str_edge_id = (str_edge_id == 0) ?
                    region_isoE.at(i).rows() - 1 : str_edge_id - 1;
            }  */     

            
        }

        //this is a guide point that we use such that if we get too close
        //to it, we re-route to the next isoline 
        //p_end = 0.5*(isoV.row(region_isoE.at(i)(end_edge_id, 0)) +
        //        isoV.row(region_isoE.at(i)(end_edge_id, 1)));

        p_end = 0.5*(iso_lines->get_coordinates(region_isoE.at(i)(end_edge_id, 0)) +
                     iso_lines->get_coordinates(region_isoE.at(i)(end_edge_id, 1)));

        //**
        //viewer.data().add_edges(isoV.row(region_isoE.at(i)(end_edge_id, 0)),
        //    isoV.row(region_isoE.at(i)(end_edge_id, 1)),
        //    Eigen::RowVector3d(0.5, 1.0, 0.3));


        int num_point_take_from_contour = 0;
        int ed = str_edge_id;

        while (true) {

            //The edge end points id
            int ed_p0(region_isoE.at(i)(ed, 0)),
                ed_p1(region_isoE.at(i)(ed, 1));

            //if (direction < 0) { std::swap(ed_p0, ed_p1); }

            //Edge mid point 
            //Eigen::VectorXd p_mid = 0.5*(isoV.row(ed_p0) + isoV.row(ed_p1));

            Eigen::VectorXd p_mid = 0.5*(
                iso_lines->get_coordinates(ed_p0) +
                iso_lines->get_coordinates(ed_p1));

            //The distance to last point
            double dd = dist(p_mid(0), p_mid(1), p_mid(2),
                p_end(0), p_end(1), p_end(2));


            //If we are too close to the end point and have taken enough points
            //from this contour 
            if (dd < tol &&
                num_point_take_from_contour>MIN_NUM_POINT_PER_CONTOUR - 3) {

                if (i == region_isoE.size() - 1) {
                    //if it is the las contour, we add the mid point and quit
                    //since no projection needed 
                    contour.row(num_point_take_from_contour) = p_mid;
                    num_point_take_from_contour++;
                    break;
                }

                //We should reroute here by projecting to the nexxt isoline  
                Eigen::VectorXd t(1), sqrD(1);                
                int edge_p0(0), edge_p1(0);

                //project_point_to_isoline(p_mid, isoV, region_isoE.at(i + 1),
                //    t, sqrD, edge_p0, edge_p1, str_edge_id);

                iso_lines->project_point_to_isoline(p_mid,
                    region_isoE.at(i + 1),
                    t, sqrD, edge_p0, edge_p1, str_edge_id);


                //str_edge_id for the next contour
                if (abs(sqrD(0, 0) - std::numeric_limits<double>::max()) < 1) {
                    std::cout << "Error:: get_spiral3d" <<
                        " Can not find intersection" << std::endl;
                    system("pause");
                }

                Eigen::VectorXd p_proj;
                p_proj.resize(3);

                //p_proj = isoV.row(edge_p0).transpose()
                //    + t(0)*(isoV.row(edge_p1).transpose() -
                //        isoV.row(edge_p0).transpose());

                p_proj = iso_lines->get_coordinates(edge_p0).transpose()
                    + t(0)*(iso_lines->get_coordinates(edge_p1).transpose() -
                        iso_lines->get_coordinates(edge_p0).transpose());


                contour.row(num_point_take_from_contour) = p_mid;
                contour.row(num_point_take_from_contour) = p_proj;
                num_point_take_from_contour += 2;

                //**
                /*if (false) {
                    //viewer.data().add_points(p_mid.transpose(),
                    //    Eigen::RowVector3d(0.0, 0.0, 1.0));
                    //viewer.data().add_points(p_proj.transpose(),
                    //    Eigen::RowVector3d(1.0, 0.0, 1.0));
                    if (false) {
                        add_edge_to_viewer(viewer, spiral, true);
                        viewer.launch();
                    }
                }
                if (false) {
                    viewer.launch();
                }*/
                
                break;
            }
            else {
                //Add the point to the contour
                contour.row(num_point_take_from_contour) = p_mid;
                num_point_take_from_contour++;
            }
            if (direction == 1) {
                ed = (ed == region_isoE.at(i).rows() - 1) ? 0 : ed + 1;
            }
            else {
                ed = (ed == 0) ? region_isoE.at(i).rows() - 1 : ed - 1;
            }
        }

        //only push to the spiral the true points 
        spiral.push_back(contour.topRows(num_point_take_from_contour - 1));
        
    }

    //add_edge_to_viewer(viewer, spiral, true);
    //viewer.launch();

    return spiral;


}
//*****************************************************************************


//********************** get_fermat3d_one_region()
std::vector<Eigen::MatrixXd> get_fermat3d_one_region(IsolinesDS * iso_lines,
    std::vector<Eigen::MatrixXd> spiral,
    const double tol) {
    //construct fermat spiral from contour spiral 
    //each row spiral is a seperate contour whose last point is connected 
    //to the first point of the next contour 
    int MAX_SIZE = 0;
    for (int i = 0; i < spiral.size(); i++) {
        MAX_SIZE += spiral.at(i).rows();
    }

    Eigen::MatrixXd inward(MAX_SIZE, 3), outward(MAX_SIZE, 3);
    int num_inward_pts(0), num_outward_pts(0);

    int current_path_id = 0;
    int nxt_pt_id = 1;
    int END = 3;
    Eigen::MatrixXd p_in = spiral.at(0).row(0);
    Eigen::MatrixXd p_out = spiral.at(0).row(spiral.at(0).rows() - END);

    Eigen::MatrixXd nxt_pt;
    nxt_pt = p_in;

#ifdef VIZ
    //**
    viewer.data().add_points(p_in, Eigen::RowVector3d(0.0, 0.0, 1.0));
    viewer.data().add_points(p_out, Eigen::RowVector3d(1.0, 0.0, 0.0));
#endif 

    Eigen::MatrixXd contour_end_pt = spiral.at(0).row(
        spiral.at(0).rows() - (END + 6));

    //**
    //viewer.data().add_points(contour_end_pt, Eigen::RowVector3d(1.0, 0, 1.0));

    //viewer.launch();

    std::vector<int> reroute_pt;//index of the rerouting points in the inward 
                                //link. used to construct the outward line

    inward.row(num_inward_pts) = nxt_pt;
    num_inward_pts++;

    int num_points_taken_from_path = 1;

    //did the inward link walked on the last contour
    bool walked_on_last_contour = false;
    while (true) {
        nxt_pt = spiral.at(current_path_id).row(nxt_pt_id);


        //The distance to last point
        double dd = dist(contour_end_pt(0),
            contour_end_pt(1),
            contour_end_pt(2),
            nxt_pt(0), nxt_pt(1), nxt_pt(2));

        if (dd < tol &&
            num_points_taken_from_path>MIN_NUM_POINT_PER_CONTOUR - 6) {

            //reroute 
            inward.row(num_inward_pts) = nxt_pt;
            num_inward_pts++;
            num_points_taken_from_path++;
            inward.row(num_inward_pts) = contour_end_pt;
            num_inward_pts++;
            num_points_taken_from_path++;

            //**
            //viewer.data().add_points(nxt_pt, Eigen::RowVector3d(0, 0, 0));
            //viewer.data().add_points(contour_end_pt,
            //    Eigen::RowVector3d(1.0, 0.0, 1.0));

            //index of reroute point
            if (current_path_id > 0) {
                reroute_pt.push_back(num_inward_pts - 1);
                //**
                //viewer.data().add_points(contour_end_pt,
                //       Eigen::RowVector3d(0, 1.0, 1.0));
            }

            if (current_path_id == spiral.size() - 1) {
                walked_on_last_contour = true;
            }

            //move to the next path
            if (current_path_id != spiral.size() - 1) {
                current_path_id++;
                num_points_taken_from_path = 0;
            }
            

            //reroute and prepare the next contour by projecting to the next
            //path 
            Eigen::VectorXd t(1, 1), sqrD(1, 1);
            int seg_st(0), seg_end(0);

            //project_point_to_isoline(contour_end_pt.transpose(),
            //    spiral.at(current_path_id), t, sqrD, seg_st, seg_end);

            iso_lines->project_point_to_isoline(contour_end_pt.transpose(),
                spiral.at(current_path_id), t, sqrD, seg_st, seg_end);


            if (abs(sqrD(0, 0) - std::numeric_limits<double>::max()) < 1) {
                std::cout << "Error:: get_fermat3d" <<
                    " Can not find intersection (0)" << std::endl;
               
                //**
                //viewer.data().add_points(contour_end_pt,
                //    Eigen::RowVector3d(0, 0, 1.0));

                //std::vector<Eigen::MatrixXd> fermat;
                //fermat.push_back(inward.topRows(num_inward_pts));
                //add_edge_to_viewer(viewer, fermat, true);
                //viewer.launch();
                system("pause");
            }
            Eigen::VectorXd p_proj(3);

            p_proj = spiral.at(current_path_id).row(seg_st).transpose()
                + t(0)*(spiral.at(current_path_id).row(seg_end).transpose() -
                    spiral.at(current_path_id).row(seg_st).transpose());

            //add the projected point to the link
            inward.row(num_inward_pts) = p_proj;
            num_inward_pts++;
            num_points_taken_from_path++;


            //**
            //viewer.data().add_points(p_proj.transpose(),
            //    Eigen::RowVector3d(0, 1.0, 0));

            //the next id is the end of the segment we projected onto since 
            // we are moving ahead 
            nxt_pt_id = seg_end;


            if (current_path_id < spiral.size() - 1) {
                //get the contour_end_pt by moving a little backward from 
                 //the projected point and then project

                int end_id = seg_st;
                for (int hh = 0; hh < END + 10; hh++) {
                    //we should not step into another path 
                    if (end_id == 0) {
                        std::cout << "Should not get here. Fix by moving to" <<
                            " another path or just terminate" << std::endl;
                        //std::vector<Eigen::MatrixXd> fermat;
                        //fermat.push_back(inward.topRows(num_inward_pts));
                        //add_edge_to_viewer(viewer, fermat, true);
                        //viewer.launch();

                        system("pause");
                    }
                    end_id--;
                    //end_id = (end_id == 0) ?
                    //    spiral.at(current_path_id).rows() - 1 : end_id - 1;
                }

                //The point after moving a little backward from the project 
                //point 
                contour_end_pt = spiral.at(current_path_id).row(end_id);

                //**           
               //viewer.data().add_points(contour_end_pt,
               //    Eigen::RowVector3d(1.0, 0.0, 1.0));

                //Now project to get the actual contour_end_pt
                Eigen::VectorXd t1(1, 1), sqrD1(1, 1);
                int st1(0), end1(0);
                
                //project_point_to_isoline(contour_end_pt.transpose(),
                //    spiral.at(current_path_id + 1), t1, sqrD1, st1, end1);

                iso_lines->project_point_to_isoline(contour_end_pt.transpose(),
                    spiral.at(current_path_id + 1), t1, sqrD1, st1, end1);


                if (abs(sqrD1(0, 0) - std::numeric_limits<double>::max()) < 1){
                    std::cout << "Error:: get_fermat3d" <<
                        " Can not find porjection" << std::endl;
                    //std::vector<Eigen::MatrixXd> fermat;
                    //fermat.push_back(inward.topRows(num_inward_pts));
                    //add_edge_to_viewer(viewer, fermat, true);
                    //viewer.launch();
                    system("pause");
                }
                contour_end_pt = (spiral.at(
                    current_path_id + 1).row(st1).transpose()
                    + t1(0)*(spiral.at(
                        current_path_id + 1).row(end1).transpose() -
                        spiral.at(
                            current_path_id + 1).row(st1).transpose())
                    ).transpose();

                //**           
                //viewer.data().add_points(contour_end_pt,
                //    Eigen::RowVector3d(1.0, 0.0, 1.0));

            }
            else {
                //If we are near the end, we don't care much about 
                //contour_end_pt but we add it to compute the distance to know
                //when should we stop

                //if walked_on_last_contour == true
                //in this case, we have reached the end, and thus the                 
                //outward link will start by projecting and walking                 
                //along path_id = spiral.size()-2
                //this happens when the number of contours is odd (but not
                //necessarily so)
                
                //if walked_on_last_contour == false 
                //in this case, we still have one final path below the last 
                //one that the inward link has reached. for that, the 
                //outward link starts by walking along the path with 
                //id = current_path_id (which the inward has not walked on
                //to yet) and we start with the seg_st because we are 
                //moving in the opposite direction 
                

                if (walked_on_last_contour) {
                    //std::cout<<"Error:: get_fermat3d() Make the number of "
                    //    << "contours even to scape this."<<std::endl;
                    //system("pause");
                    
                    //here we project seg_end                    
                    //viewer.data().add_points(
                    //    spiral.at(current_path_id).row(seg_end),
                    //    Eigen::RowVector3d(1.0, 0, 0));

                   

                    Eigen::VectorXd t2(1, 1), sqrD2(1, 1);                    
                    //project_point_to_isoline(
                    //    spiral.at(current_path_id).row(seg_end).transpose(),
                    //    spiral.at(current_path_id - 1),
                    //    t2, sqrD2, seg_st, seg_end);

                    iso_lines->project_point_to_isoline(
                        spiral.at(current_path_id).row(seg_end).transpose(),
                        spiral.at(current_path_id - 1),
                        t2, sqrD2, seg_st, seg_end);


                    if (abs(sqrD2(0, 0) 
                        - std::numeric_limits<double>::max()) < 1) {
                        std::cout << "Error:: get_fermat3d" <<
                            " Can not find intersection (1)" << std::endl;
                        system("pause");
                    }         
                    current_path_id--;
                    reroute_pt.pop_back();
                }
                
                nxt_pt_id = seg_st;

                //**
                //viewer.data().add_points(spiral.at(
                //    current_path_id).row(seg_st),
                //    Eigen::RowVector3d(0.5, 0.25, 0.5));

                /*viewer.data().add_points(spiral.at(
                    current_path_id).row(seg_end),
                    Eigen::RowVector3d(0, 0, 0.5));*/

                //viewer.launch();
                break;
               
            }


            //**           
            //viewer.data().add_points(contour_end_pt,
            //    Eigen::RowVector3d(0.5, 0.5, 0.5));

        }
        else {
            inward.row(num_inward_pts) = nxt_pt;
            num_inward_pts++;
            num_points_taken_from_path++;


            //**
            //viewer.data().add_points(nxt_pt, Eigen::RowVector3d(0, 0, 0));


            if (nxt_pt_id == spiral.at(current_path_id).rows() - 1) {
                //move to next conour 
                if (current_path_id == spiral.size() - 1) {
                    //this was the last one alreadt
                    break;
                }

                current_path_id++;
                nxt_pt_id = 0;
            }
            else {
                nxt_pt_id++;
            }
        }
    }

    

    //outward link
    //follow the same logic but with going in the opposite direction 
    if (reroute_pt.size() > 0) {
        contour_end_pt = inward.row(reroute_pt.back());
        reroute_pt.pop_back();
    }
    else {
        //if there is no enough rerouting points i.e, the inward linke is just
        //one path, the the outward is one path also ends at p_out
        contour_end_pt = p_out;
    }
    num_points_taken_from_path = 0;

    while (true) {
        nxt_pt = spiral.at(current_path_id).row(nxt_pt_id);

        //The distance to last point
        double dd = dist(contour_end_pt(0),
            contour_end_pt(1),
            contour_end_pt(2),
            nxt_pt(0), nxt_pt(1), nxt_pt(2));
        if (dd < tol &&
            num_points_taken_from_path > MIN_NUM_POINT_PER_CONTOUR - 3) {

            outward.row(num_outward_pts) = nxt_pt;
            num_outward_pts++;

            //**           
            //viewer.data().add_points(outward.row(num_outward_pts - 1),
            //    Eigen::RowVector3d(0.5, 0.5, 0.5));

            //decrement the path id
            if (current_path_id > 1) {
                //we don't want reach path with id =0 because we know that this 
                //path is already occupoed by the inward link
                current_path_id--;
                num_points_taken_from_path = 0;
            }
            else {
                //add the exit point and exit               
                outward.row(num_outward_pts) = contour_end_pt;
                num_outward_pts++;
                break;
            }

            

            //reroute by projection on the current path (we already updated it)
            Eigen::VectorXd t(1, 1), sqrD(1, 1);
            int seg_st(0), seg_end(0);           
            
            //project_point_to_isoline(nxt_pt.transpose(),
            //    spiral.at(current_path_id), t, sqrD, seg_st, seg_end);

            iso_lines->project_point_to_isoline(nxt_pt.transpose(),
                spiral.at(current_path_id), t, sqrD, seg_st, seg_end);


            if (abs(sqrD(0, 0) - std::numeric_limits<double>::max()) < 1) {
                std::cout << "Error:: get_fermat3d" <<
                    " Can not find intersection" << std::endl;
                system("pause");
            }
            Eigen::VectorXd p_proj(3);

            p_proj = spiral.at(current_path_id).row(seg_st).transpose()
                + t(0)*(spiral.at(current_path_id).row(seg_end).transpose() -
                    spiral.at(current_path_id).row(seg_st).transpose());

            //add the projected point to the link
            outward.row(num_outward_pts) = p_proj;
            num_outward_pts++;
            num_points_taken_from_path++;

            //**           
            //viewer.data().add_points(outward.row(num_outward_pts - 1),
            //    Eigen::RowVector3d(0.1, 1.0, 1.0));


            //the next id is the start of the segment since we are moving 
            //in the opposite direction
            nxt_pt_id = seg_st;

            //update the contour_end_pt so we know where to stop next time 
            if (reroute_pt.size() == 0) {
                contour_end_pt = p_out;
            }
            else {
                contour_end_pt = inward.row(reroute_pt.back());
                reroute_pt.pop_back();
            }
        }
        else {
            //if we are still far from the rerouting point, we add the next_pt 
            //and decrement the point id (and the path if necessary)
            outward.row(num_outward_pts) = nxt_pt;
            num_outward_pts++;
            num_points_taken_from_path++;

            //**       
            //viewer.data().add_points(outward.row(num_outward_pts - 1),
            //    Eigen::RowVector3d(0.5, 0.5, 0.5));
            
            if (nxt_pt_id == 0) {
                if (current_path_id == 0) { break; }
                current_path_id--;
                nxt_pt_id = spiral.at(current_path_id).rows() - 1;
            }
            else {
                nxt_pt_id--;
            }

        }
    }

    std::vector<Eigen::MatrixXd> fermat;    
    fermat.push_back(inward.topRows(num_inward_pts));
    fermat.push_back(outward.topRows(num_outward_pts));
    
#ifdef VIZ
    add_edge_to_viewer(viewer, fermat, true);
#endif 

    //viewer.launch();

    return fermat;


}
//*****************************************************************************


//********************** fermat3d()
std::vector<std::vector<Eigen::MatrixXd>> fermat3d(
    IsolinesDS * iso_lines, const double tol) {


    //Create the fermat spiral for a bunch of spirallable regions 
    int num_regions = iso_lines->get_num_spiral_region();

    //std::vector<std::vector<Eigen::MatrixXd>> spiral(num_regions);
    //for (int r = 0; r < num_regions; r++) {
    //    spiral.at(r) = get_spiral3d_one_region(iso_lines, r, tol);
    //}
    std::vector<std::vector<Eigen::MatrixXd>> spiral;
    for (int r = 0; r < num_regions; r++) { 
        if (iso_lines->get_num_contours(r) < 6) { continue; }
        spiral.push_back(get_spiral3d_one_region(iso_lines, r, tol));
    }

    //std::vector<std::vector<Eigen::MatrixXd>> fermat(spiral.size());
    //for (int s = 0; s < spiral.size(); s++) {
    //    fermat.at(s) = get_fermat3d_one_region(iso_lines, spiral.at(s), tol);
    //}
    std::vector<std::vector<Eigen::MatrixXd>> fermat;
    for (int s = 0; s < spiral.size(); s++) {
        fermat.push_back(get_fermat3d_one_region(iso_lines,spiral.at(s), tol));
    }

    return fermat;

}
//*****************************************************************************
#endif /*_FERMAT3D_IMP_*/
