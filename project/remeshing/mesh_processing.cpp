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
#include "mesh_processing.h"
# include "HoleDigger.h"
#include <set>
#include <vector>
#include <map>

using namespace std;

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using surface_mesh::Vec3;
using surface_mesh::Normal;
using surface_mesh::cross;
using surface_mesh::norm;
using std::min;
using std::max;
using std::cout;
using std::endl;

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

MeshProcessing::~MeshProcessing() {
    // TODO
}

void MeshProcessing::hullify(double magnitude, double reduce)
{
    compute_mesh_properties();
    std::map<Mesh::Vertex,Mesh::Vertex> innerVertexes;
    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    std::vector<Mesh::Vertex> inner;


    if(reduce > 0){
        for(auto v : mesh_.vertices()){
            auto innerVPos = mesh_.position(v) - normals[v]*reduce;
            //mesh_.position(v) = innerVPos;
            mesh_.position(v) = Point(innerVPos.x, innerVPos.z, innerVPos.y);
        }
        uniform_smooth(1);
    }
    compute_mesh_properties();
    for(auto f : mesh_.faces()) {

        inner.clear();
        auto vs = mesh_.vertices(f);
        for(auto v : vs) {
            if(innerVertexes.count(v) == 0){
                auto innerVPos = mesh_.position(v) - normals[v]*magnitude;

                auto innerV = mesh_.add_vertex(innerVPos);
                innerVertexes[v] = innerV;
                inner.insert(inner.begin(),innerV);

                mesh_.position(innerV)  = innerVPos;


            }else{
                inner.insert(inner.begin(),innerVertexes[v]);
            }
        }
        mesh_.add_face(inner);

    }
    compute_mesh_properties();
    for(auto e : mesh_.halfedges()){
        if(!mesh_.face(e).is_valid()){

            std::vector<Mesh::Vertex> newFace;

            auto f =mesh_.from_vertex(e);
            auto t = mesh_.to_vertex(e);
            auto fi = innerVertexes[mesh_.from_vertex(e)];
            auto ti = innerVertexes[mesh_.to_vertex(e)];

            newFace.push_back(f);
            newFace.push_back(t);
            newFace.push_back(fi);
            mesh_.add_face(newFace);
            newFace.clear();

            newFace.push_back(ti);
            newFace.push_back(fi);
            newFace.push_back(t);

            mesh_.add_face(newFace);





        }
    }
    compute_mesh_properties();
}
void MeshProcessing::uniform_smooth(unsigned int n_iters) {

    Vec3 *Lu = new Vec3[mesh_.vertices_size()];
    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        // ------------- IMPLEMENT HERE ---------
        for (auto vertex: mesh_.vertices()) {
            Lu[vertex.idx()] = Vec3(0, 0, 0);
            if (mesh_.is_boundary(vertex))
                continue;
            int n = mesh_.valence(vertex);
            if (n == 0)
                continue;
            for (auto neighbor_vertex: mesh_.vertices(vertex))
                Lu[vertex.idx()] += mesh_.position(neighbor_vertex) - mesh_.position(vertex);
            Lu[vertex.idx()] /= n * 2;
        }
        for (auto vertex: mesh_.vertices()) {
            mesh_.position(vertex) += Lu[vertex.idx()];
        }

    }
    delete []Lu;

    // update face and vertex normals
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}





// ========================================================================
// EXERCISE 1.2
// ========================================================================
void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // ------------- IMPLEMENT HERE ---------

    calc_weights();
    for (auto vertex: mesh_.vertices()) {
        if (mesh_.is_boundary(vertex))
            continue;
        int n = mesh_.valence(vertex);
        if (n == 0)
            continue;
        Vec3 v_sum(0,0,0);
        for (auto halfedge: mesh_.halfedges(vertex)) {
            auto edge = mesh_.edge(halfedge);
            auto neighbor = mesh_.to_vertex(halfedge);
            v_sum += e_weight[edge] * (mesh_.position(neighbor) - mesh_.position(vertex));
        }
        v_sum *= v_weight[vertex];
        v_curvature[vertex] = norm(v_sum);
    }
}

// ========================================================================
// EXERCISE 1.3
// ========================================================================
void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // ------------- IMPLEMENT HERE ---------

    calc_vertices_weights();
    for (auto vertex: mesh_.vertices()) {
        if (mesh_.is_boundary(vertex))
            continue;
        Scalar sum = 2 * M_PI;
        for (auto face: mesh_.faces(vertex)) {
            Scalar l[3];
            for (auto face_edge: mesh_.halfedges(face)){
                if (mesh_.to_vertex(face_edge) == vertex)
                    l[0] = mesh_.edge_length(mesh_.edge(face_edge));
                else if (mesh_.from_vertex(face_edge) == vertex)
                    l[1] = mesh_.edge_length(mesh_.edge(face_edge));
                else
                    l[2] = mesh_.edge_length(mesh_.edge(face_edge));
            }
            std::complex<double> z((l[0] * l[0] + l[1] * l[1] - l[2] * l[2]) / 2 / l[0] / l[1], 0);
            sum -= acos(z).real();
        }
        v_gauss_curvature[vertex] = sum * 2 * v_weight[vertex];
    }
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
surface_mesh::Color const MeshProcessing::non_mark_ = surface_mesh::Color(1.0, 1.0, 1.0);

void MeshProcessing::delete_marked_faces()
{
    std::vector<Mesh::Face> marked_faces;
    marked_faces.reserve(mesh_.n_faces());

    for (auto face : mesh_.faces())
    {
        if (is_marked(face)) marked_faces.push_back(face);
    }

    marked_faces.shrink_to_fit();

    for (auto face : marked_faces)
    {
        if (!mesh_.is_deleted(face))
        {
            mesh_.delete_face(face);
        }
    }

    mesh_.garbage_collection();
}

bool MeshProcessing::is_marked(Mesh::Vertex vertex)
{
    Mesh::Vertex_property<Color> colors = mesh_.vertex_property<Color>("v:color");
    return colors[vertex] != non_mark_;
}

bool MeshProcessing::is_marked(Mesh::Face face)
{
    if (mesh_.is_deleted(face))
    {
        return false;
    }

    unsigned int marked_vertices_count = 0;

    for (auto vertex : mesh_.vertices(face))
    {
        if (is_marked(vertex)) ++marked_vertices_count;
    }

    // We assume that we deal with a triangular mesh.
    // Therefore, if more than one vertex out of three
    // is marked, then we mark the face as deleted.
    //
    return marked_vertices_count > 1;
}
void MeshProcessing::digHole(double distance, double diameter){
    auto digger = HoleDigger(mesh_);
    digger.sample(distance);
    digger.digHole(diameter);
}
/*
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
*/
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
void MeshProcessing::save_mesh(string file)
{
    bool const written_successfully = mesh_.write(file);

    if (written_successfully) std::cout << "Saved successfully" << std::endl;
    else                      std::cout << "Saving failed" << std::endl;
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
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
                 8 /* max */);
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
                                  Mesh::Vertex_property<Color> color_prop, Scalar min_value,
                                  Scalar max_value, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    if (min_value == 0.0 && max_value == 0.0) {
        // discard upper and lower bound
        unsigned int n = values.size()-1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        min_value = values[i];
        max_value = values[n-1-i];
    }

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

void MeshProcessing::octreeSample(double distance, vector<Mesh::Vertex> &input,
                                  vector<Mesh::Vertex> &to_pick) {
    if (input.size() == 0)
        return;
    if (input.size() == 1) {
        to_pick.push_back(input[0]);
        return;
    }
    float ma[3];
    float mi[3];
    float length;
    bool first_vertex = true;
    for (auto vertex : input) {
        auto pos =mesh_.position(vertex);
        if (first_vertex) {
            first_vertex = false;
            for (int i = 0; i < 3; i++) {
                ma[i] = pos.data()[i];
                mi[i] = pos.data()[i];
            }
        }
        else {
            for (int i = 0; i < 3; i++) {
                ma[i] = max(ma[i], pos.data()[i]);
                mi[i] = min(mi[i], pos.data()[i]);
            }
        }
    }
    length = 0;
    for (int i = 0; i < 3; i++)
        length += pow(ma[i] - mi[i], 2);
    length = sqrt(length);
    if (length > distance) {
        vector<Mesh::Vertex> subRange[8];
        for (int i = 0; i < 8; i++)
            subRange[i].clear();
        for (auto vertex : input) {
            auto pos =mesh_.position(vertex);
            int index = 0;
            for (int i = 0; i < 3; i++)
                if (pos.data()[i] < (mi[i] + ma[i]) / 2)
                    index += (1 << i);
            subRange[index].push_back(vertex);
        }
        for (int i = 0; i < 8; i++)
            octreeSample(distance, subRange[i], to_pick);
    }
    else {
        Point center;
        for (int i = 0; i < 3; i++)
            center[i] = (mi[i] + ma[i]) / 2;
        auto to_add = findNearest(center, input);
        to_pick.push_back(to_add);
    }
}

Mesh::Vertex MeshProcessing::findNearest(Point center, vector<Mesh::Vertex> &range) {
    float les = -1;
    Mesh::Vertex result;
    for (auto vertex : range) {
        float dis = norm(mesh_.position(vertex) - center);
        if (les < 0) {
            les = dis;
            result = vertex;
        }
        else if (dis < les) {
            les = dis;
            result = vertex;
        }
    }
    return result;
}


}


