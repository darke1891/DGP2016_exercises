//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

using std::min;
using std::max;
using namespace surface_mesh;

typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Viewer::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature = mesh.vertex_property<Scalar>("v:unicurvature", 0);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate mean curvature using
    // the length of the uniform Laplacian approximation
    // Save your approximation in unicurvature vertex property of the mesh.
    // ------------- IMPLEMENT HERE ---------
    for (auto vertex: mesh.vertices()) {
        if (mesh.is_boundary(vertex))
            continue;
        int n = mesh.valence(vertex);
        if (n == 0)
            continue;
        Vec3 v_sum(0,0,0);
        for (auto neighbor_vertex: mesh.vertices(vertex))
            v_sum += mesh.position(neighbor_vertex) - mesh.position(vertex);
        v_sum /= n;
        v_unicurvature[vertex] = norm(v_sum) / 2;
    }
}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Viewer::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
    Mesh::Vertex_property<Scalar>  v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // ------------- IMPLEMENT HERE ---------
    calc_weights();
    for (auto vertex: mesh.vertices()) {
        if (mesh.is_boundary(vertex))
            continue;
        int n = mesh.valence(vertex);
        if (n == 0)
            continue;
        Vec3 v_sum(0,0,0);
        for (auto halfedge: mesh.halfedges(vertex)) {
            auto edge = mesh.edge(halfedge);
            auto neighbor = mesh.to_vertex(halfedge);
            v_sum += e_weight[edge] * (mesh.position(neighbor) - mesh.position(vertex));
        }
        v_sum *= v_weight[vertex];
        v_curvature[vertex] = norm(v_sum);
    }
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Viewer::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // ------------- IMPLEMENT HERE ---------
    calc_weights();
    for (auto vertex: mesh.vertices()) {
        if (mesh.is_boundary(vertex))
            continue;
        Scalar sum = 2 * M_PI;
        for (auto face: mesh.faces(vertex)) {
            Scalar l[3];
            for (auto face_edge: mesh.halfedges(face)){
                if (mesh.to_vertex(face_edge) == vertex)
                    l[0] = mesh.edge_length(mesh.edge(face_edge));
                else if (mesh.from_vertex(face_edge) == vertex)
                    l[1] = mesh.edge_length(mesh.edge(face_edge));
                else
                    l[2] = mesh.edge_length(mesh.edge(face_edge));
            }
            sum -= acos((l[0] * l[0] + l[1] * l[1] - l[2] * l[2]) / 2 / l[0] / l[1]);
        }
        v_gauss_curvature[vertex] = sum * 2 * v_weight[vertex];
    }
}

// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Viewer::uniform_smooth(unsigned int n_iters) {

    Vec3 *Lu = new Vec3[mesh.vertices_size()];
    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        // ------------- IMPLEMENT HERE ---------
        for (auto vertex: mesh.vertices()) {
            Lu[vertex.idx()] = Vec3(0, 0, 0);
            if (mesh.is_boundary(vertex))
                continue;
            int n = mesh.valence(vertex);
            if (n == 0)
                continue;
            for (auto neighbor_vertex: mesh.vertices(vertex))
                Lu[vertex.idx()] += mesh.position(neighbor_vertex) - mesh.position(vertex);
            Lu[vertex.idx()] /= n * 2;
        }
        for (auto vertex: mesh.vertices()) {
            mesh.position(vertex) += Lu[vertex.idx()];
        }

    }

    // update face and vertex normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Viewer::smooth(unsigned int n_iters) {

    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // Perform Laplace-Beltrami smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        // 2) for each non-boundary vertex, update its position using the normalized Laplace-Beltrami operator
        //    (Hint: use the precomputed edge weights in the edge property "e:weight")
        // ------------- IMPLEMENT HERE ---------
    }


    // update face and vertex normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}

// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::uniform_laplacian_enhance_feature(int enhancement_smoothing_iterations,
                                               float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the uniform Laplacian operator:
    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------

    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::laplace_beltrami_enhance_feature(int enhancement_smoothing_iterations,
                                              float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the Laplace-Beltrami operator:
    // 1) perform Laplace-Beltrami smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
