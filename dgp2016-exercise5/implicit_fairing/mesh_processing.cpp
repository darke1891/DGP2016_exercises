//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#include "mesh_processing.h"
#include <set>

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using std::min;
using std::max;
using std::cout;
using std::endl;

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

// ============================================================================
// EXERCISE 5.1
// ============================================================================
void MeshProcessing::implicit_smoothing(const double timestep) {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd B(n,3);

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.1 HERE
    // ========================================================================
    for (auto vertex: mesh_.vertices()) {
        auto areaInvInv = 1 / area_inv[vertex];
        auto vertex_position = mesh_.position(vertex);
        if (mesh_.is_boundary(vertex)) {
            B.row(vertex.idx()) =  Eigen::RowVector3d(vertex_position[0], vertex_position[1], vertex_position[2]);
            triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(),1));
        }
        else {
            triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(),areaInvInv));
            B.row(vertex.idx()) =  Eigen::RowVector3d(vertex_position[0], vertex_position[1], vertex_position[2]) * areaInvInv;
            for (auto h: mesh_.halfedges(vertex)){
                auto d = timestep * cotan[mesh_.edge(h)];
                auto vertex2 = mesh_.to_vertex(h);
                if (mesh_.is_boundary(vertex2)){
                    auto vertex2_position = mesh_.position(vertex2);
                    B.row(vertex.idx()) += Eigen::RowVector3d(vertex2_position[0], vertex2_position[1], vertex2_position[2]) * d;
                }
                else
                    triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex2.idx(), -d));
                triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(), d));
            }
        }
    }


    // build sparse matrix from triplets
    A.setFromTriplets(triplets.begin(), triplets.end());

    // solve A*X = B
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(A);
    if (solver.info () != Eigen::Success) {
        printf("linear solver init failed.\n");
    }

    Eigen::MatrixXd X = solver.solve(B);

    if (solver.info () != Eigen::Success) {
        printf("linear solver failed.\n");
    }

    // copy solution
    for (int i = 0; i < n; ++i)
    {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = X(i, dim);
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

// ============================================================================
// EXERCISE 5.2
// ============================================================================
void MeshProcessing::minimal_surface() {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");
    auto points_init = mesh_init_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> L (n, n);
    Eigen::MatrixXd rhs (Eigen::MatrixXd::Zero (n, 3));

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets_L;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.2 HERE
    // ========================================================================

    std::vector< Eigen::Triplet<double> > triplets;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.1 HERE
    // ========================================================================
    for (auto vertex: mesh_.vertices()) {
        auto areaInv = 1 / area_inv[vertex];
        if (mesh_.is_boundary(vertex)) {
            triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(),1));
        }
        else {
            triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(),areaInv));
            for (auto h: mesh_.halfedges(vertex)){
                auto d =  cotan[mesh_.edge(h)];
                auto vertex2 = mesh_.to_vertex(h);
                if (!mesh_.is_boundary(vertex2))
                    triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex2.idx(), -d * areaInv));
                triplets.push_back(Eigen::Triplet<double>(vertex.idx(),vertex.idx(), d * areaInv));
            }
        }
    }
    L.setFromTriplets (triplets_L.begin (), triplets_L.end ());

    // solve A*X = B
    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver(L);
    if (solver.info () != Eigen::Success) {
        printf("linear solver init failed.\n");
    }

    Eigen::MatrixXd X = solver.solve(rhs);
    if (solver.info () != Eigen::Success) {
        printf("linear solver failed.\n");
    }

    // copy solution
    for (int i = 0; i < n; ++i) {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim) {
            points[v][dim] += 1. * (X(i, dim) - points[v][dim]);
        }
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}
void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
}

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
}

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
}

void MeshProcessing::uniform_smooth(const unsigned int iterations) {

    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    }
}

void MeshProcessing::smooth(const unsigned int iterations) {

    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    }
}

void MeshProcessing::uniform_laplacian_enhance_feature(const unsigned int iterations,
                                                       const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
}

void MeshProcessing::laplace_beltrami_enhance_feature(const unsigned int iterations,
                                                      const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
}

void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges())
    {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0))
        {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::compute_mesh_properties() {
    Mesh::Vertex_property<Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for(auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
                mesh_.position(v).y,
                mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
                vertex_normal[v].y,
                vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
                v_color_valence[v].y,
                v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                v_color_unicurvature[v].y,
                v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
                v_color_curvature[v].y,
                v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                v_color_gaussian_curv[v].y,
                v_color_gaussian_curv[v].z;
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
                                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}

MeshProcessing::~MeshProcessing() {}
}