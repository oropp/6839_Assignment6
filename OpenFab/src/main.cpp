#include "voxelizer.hpp"
#include "marching_cube.hpp"
#include <iostream>
#include <ctime>
#include <string>

int main(int argc, char* argv[]) {
    // Usage:
    // Slow version:
    // ./OpenFab bunny.stl 2.0
    // ./OpenFab fandisk.stl 0.05
    // ./OpenFab spot.stl 0.125
    // ./OpenFab dragon.stl 0.05
    // Fast version:
    // ./OpenFab bunny.stl 2.0 fast
    // ./OpenFab fandisk.stl 0.05 fast
    // ./OpenFab spot.stl 0.125 fast
    // ./OpenFab dragon.stl 0.05 fast
    // Approximation:
    // ./OpenFab bunny_with_hole.stl 2.0 approx
    // ./OpenFab spot_with_whole.stl 0.125 approx
    // Marching cube version:
    // ./OpenFab bunny_voxel_info.txt
    // ./OpenFab dragon_voxel_info.txt
    // ./OpenFab fandisk_voxel_info.txt
    // ./OpenFab spot_voxel_info.txt
    if (argc == 2) {
        // Marching cube version.
        const std::string info_file(argv[1]);
        mesh::MarchingCube<double> mc(std::string(PROJECT_SOURCE_DIR) + "/data/assignment2/" + info_file);
        mc.BuildMesh();
        const std::string name = info_file.substr(0, info_file.size() - std::string("_voxel_info.txt").size());
        mc.ExportMeshToFile(name + "_mc.stl");
        return 0;
    }

    int t0 = std::clock();
    const std::string stl_name(argv[1]);
    const double dx = std::stod(argv[2]);
    mesh::Voxelizer<double> voxelizer(std::string(PROJECT_SOURCE_DIR) + "/data/assignment2/" + stl_name, dx);
    int t1 = std::clock();
    std::cout << "load mesh success... " << (double)(t1 - t0) / 1000000.0 << " seconds." << std::endl;
    std::cout << "Bounding box: " << voxelizer.pmin().transpose() << ", " << voxelizer.pmax().transpose() << std::endl;
    std::cout << "Number of voxels: " << voxelizer.voxel_num().transpose() << std::endl;
    if (argc == 3) {
        voxelizer.BasicVoxelization();
    } else {
        const std::string flag(argv[3]);
        if (flag == "fast")
            voxelizer.AdvancedVoxelization();
        else if (flag == "approx")
            voxelizer.AdvancedVoxelizationWithApproximation();
        else {
            std::cout << "ERROR: unexpected flag" << std::endl;
            exit(0);
        }
    }
    std::cout << "Voxelization done..." << std::endl;
    // Export results to mesh.
    const std::string stl_prefix = stl_name.substr(0, stl_name.size() - 4);
    const std::string voxel_file_name = stl_prefix + "_voxel.stl";
    std::cout << "Saving results to " << voxel_file_name << std::endl;
    voxelizer.WriteVoxelToMesh(voxel_file_name);
    std::cout << "Results saved..." << std::endl;
    return 0;
}
