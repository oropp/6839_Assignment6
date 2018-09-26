#include "GCodeConverter.hpp"
#include "FabSlicer.hpp"
#include "generate_commands.hpp"
#include "server.hpp"

// sample main code for FabSlicer, you should create your own main code
 int main() {
     // load a stl file and convert to a connected mesh (obj)
     mesh::TriMesh<double> tri_mesh(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.stl", 0.01);

     // you can visualize your mesh as an obj file to check if it's loaded correctly
     tri_mesh.WriteToObj(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.obj");

     // create a FabSlicer instance
     //fab_translation::FabSlicer<double> fab(tri_mesh, 0.0, 2.0, 0.03, 0.05);
    fab_translation::FabSlicer<double> fab(tri_mesh, 0.0, 2.0, 0.03, 0.05);

     std::vector<std::vector<std::vector<Eigen::Vector3d>>> contour;
     std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> infill_edges;

     fab.RunTranslation(contour, infill_edges);

     // visualize your results with ply format
     //std::string contour_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-contour.ply";
     //fab.VisualizeContour(contour_file, 0.001, contour);

     //std::string infill_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-infill.ply";
     //fab.VisualizeInfill(infill_file, 0.001, infill_edges);

     std::vector<std::vector<Eigen::Vector3d>> paths;
     paths.clear();
     network_communication::GenerateCommands::ResetCommand();
     for (int i = 0; i < contour.size(); ++i) {
         for (int j = 0; j < contour[i].size(); ++j) {
             std::vector<Eigen::Vector3d> path;
             path.clear();
             for (int k = 0; k < contour[i][j].size(); ++k) {
                 float x = contour[i][j][k][0];
                 float y = contour[i][j][k][2];
                 float z = contour[i][j][k][1];
                 path.push_back(Eigen::Vector3d(x, y, z));
             }
             paths.push_back(path);
         }
         for (int j = 0; j < infill_edges[i].size(); ++j) {
             std::vector<Eigen::Vector3d> path;
             path.clear();

             float x_first = infill_edges[i][j].first[0];
             float y_first = infill_edges[i][j].first[2];
             float z_first = infill_edges[i][j].first[1];

             float x_sec = infill_edges[i][j].second[0];
             float y_sec = infill_edges[i][j].second[2];
             float z_sec = infill_edges[i][j].second[1];
             path.push_back(Eigen::Vector3d(x_first, y_first, z_first));
             path.push_back(Eigen::Vector3d(x_sec, y_sec, z_sec));
             paths.push_back(path);
         }
     }
    fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);

     return 0;
 }

// sample code for GCode generation
//int main() {
//    std::vector<std::vector<Eigen::Vector3d>> paths;
//    paths.clear();
//    //// reset UI
//    network_communication::GenerateCommands::ResetCommand();
//    //// translate to gcode
//    std::vector<Eigen::Vector3d> path;
//    path.clear();
//    for (int i = 0; i < 10; ++i) {
//	    for (int j = 0; j < 36; ++j) {
//	        path.push_back(Eigen::Vector3d(sin(j*(6.28/35.0))*i*0.01, 0.01, i*0.01*cos(j*(6.28/35.0))));
//	    }
//	    paths.push_back(path);
//	    path.clear();
//    }
//    for (int i = 0; i < 30; ++i) {
//	    for (int j = 0; j < 36; ++j) {
//	        path.push_back(Eigen::Vector3d(sin(j*(6.28/35.0))*0.1, i*0.01, cos(j*(6.28/35.0))*0.1));
//	    }
//	    paths.push_back(path);
//	    path.clear();
//    }
//
//    fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);
//
//    return 0;
//}

/* Implement your code here */
// int main() {
//    
// }
