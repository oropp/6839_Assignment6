#include "GCodeConverter.hpp"
#include "FabSlicer.hpp"
#include "generate_commands.hpp"
#include "server.hpp"

int main() {
    int t0 = std::clock();
    mesh::TriMesh<double> tri_mesh(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.stl", 0.01);
    int t1 = std::clock();
    printf("load mesh success... %.6lf seconds\n", (double)(t1 - t0) / 1000000.0);

    tri_mesh.WriteToObj(std::string(PROJECT_SOURCE_DIR) + "/data/bunny_watertight.obj");
    fab_translation::FabSlicer<double> fab(tri_mesh, 0.0, 2.0, 0.03, 0.05);

    std::vector<std::vector<std::vector<Eigen::Vector3d>>> contour;
    std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> infill_edges;

    fab.RunTranslation(contour, infill_edges);

    std::string contour_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-contour.ply";
    fab.VisualizeContour(contour_file, 0.001, contour);

    std::string infill_file = std::string(PROJECT_SOURCE_DIR) + "/data/bunny-infill.ply";
    fab.VisualizeInfill(infill_file, 0.001, infill_edges);

    std::vector<std::vector<Eigen::Vector3d>> paths;
    paths.clear();
    network_communication::GenerateCommands::ResetCommand();
    // translate to gcode
    for (int i = 0;i < contour.size();++i) { // for each layer
        // first print contour
        paths.clear();
        for (int j = 0;j < contour[i].size();++j) {
            // first, move to the first contour point without extrude
            std::vector<Eigen::Vector3d> path;
            path.clear();
            path.push_back(contour[i][j][0]);
            // then print contour
            path.clear();
            for (int k = 0;k < contour[i][j].size();++k)
                path.push_back(contour[i][j][k]);
            path.push_back(contour[i][j][0]);
            for (int k = 0;k < path.size();++k)
                std::swap(path[k](1), path[k](2));
            paths.push_back(path);
        }
        fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);
        paths.clear();
        // then print infill
        for (int j = 0;j < infill_edges[i].size();++j) {
            // first, move to the end of segment
            std::vector<Eigen::Vector3d> path;
            path.clear();
            path.push_back(infill_edges[i][j].first);
            path.push_back(infill_edges[i][j].second);
            for (int k = 0;k < path.size();++k)
                std::swap(path[k](1), path[k](2));
            paths.push_back(path);
        }
        fab_translation::GCodeConverter::ConvertToGCode(paths, nullptr);
    }


    return 0;
}
