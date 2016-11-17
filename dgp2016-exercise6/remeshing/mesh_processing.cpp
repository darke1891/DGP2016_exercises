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
#include <set>

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

void MeshProcessing::remesh (const REMESHING_TYPE &remeshing_type,
                             const int &num_iterations) {
    calc_weights ();
    calc_mean_curvature ();
    calc_uniform_mean_curvature ();
    calc_gauss_curvature ();
    calc_target_length (remeshing_type);

    // main remeshing loop
    for (int i = 0; i < num_iterations; ++i)
    {
        split_long_edges ();
        collapse_short_edges ();
        equalize_valences ();
        tangential_relaxation ();
    }
}

void MeshProcessing::calc_target_length (const REMESHING_TYPE &remeshing_type)
{
    Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator  vv_c, vv_end;
    Scalar                   length;
    Scalar                   mean_length;
    Scalar                   H;
    Scalar                   K;

    int curvatureSmoothing = 1; //Number of time the curvature is smoothed before being used for remeshing.

    Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Vertex_property<Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    Mesh::Vertex_property<Scalar> target_new_length  = mesh_.vertex_property<Scalar>("v:newlength", 0);


    if (remeshing_type == AVERAGE)
    {
        length = 0;
        int n = 0;
        for (auto edge: mesh_.edges()) {
            length += mesh_.edge_length(edge);
            n++;
        }
        if (n == 0)
            length = 1;
        else
            length /= n;
        length *= 0.8;
        for (auto vertex: mesh_.vertices()) {
            target_length[vertex] = length;
        }
    }
    else if (remeshing_type == CURV)
    {

        length = 0;
        int n = 0;
        for (auto edge: mesh_.edges()) {
            length += mesh_.edge_length(edge);
            n++;
        }
        if (n == 0)
            length = 1;
        else
            length /= n;
        length *= 0.8;
        for (auto vertex: mesh_.vertices()) {
            target_length[vertex] = length;
        }

        //Only one of these is useful
        calc_mean_curvature();
        cout <<"a"<<endl;
        //Smooth curvature
        for(int i = 0; i < curvatureSmoothing;i++){
            for(auto vertex : mesh_.vertices()) {
                auto valence = mesh_.valence(vertex);

                for(auto neighboor : mesh_.vertices(vertex)){
                    curvature[vertex] += curvature[neighboor];
                }
                // Average between the curvature of the vertex and all its neighbors. So there are valence + 1 points.

                curvature[vertex]/=valence+1.0;

            }
        }
        auto average_smoothed_curv = 0.0;
        for(auto vertex : mesh_.vertices()) {
            average_smoothed_curv += curvature[vertex];
        }
        average_smoothed_curv /= mesh_.n_vertices();

        for (auto vertex : mesh_.vertices())
        {
            //The user have no way to chose a mean length so I use this.
            auto c = average_smoothed_curv/curvature[vertex];
            target_length[vertex] *=  c;
            //If he would, I would use this :
            /*
             * target_length[vertex]  = target_mean * average_smoothed_curv/curvature[vertex];
             */
        }
    }
    else if (remeshing_type == HEIGHT)
    {
        length = 0;
        int n = 0;
        for (auto edge: mesh_.edges()) {
            length += mesh_.edge_length(edge);
            n++;
        }
        if (n == 0)
            length = 1;
        else
            length /= n;
        Scalar min, max;
        auto vertex0 = *(mesh_.vertices().begin());
        min = mesh_.position(vertex0)[1];
        max = mesh_.position(vertex0)[1];
        for (auto vertex: mesh_.vertices()) {
            Scalar h = mesh_.position(vertex)[1];
            if (min > h)
                min = h;
            if (max < h)
                max = h;
        }
        for (auto vertex: mesh_.vertices()) {
            Scalar h = mesh_.position(vertex)[1];
            Scalar scale = (h - min) / (max - min) * 4 + 1;
            target_length[vertex] = length * scale;
        }
    }

}

void MeshProcessing::split_long_edges ()
{
    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    bool finished;
    int i;

    for (finished = false, i = 0; !finished && i < 100; ++i)
    {
        // We presume that there is no more splitting left to do.
        //
        finished = true;


        //
        // Lookup and split the edges.
        //

        for (auto edge : mesh_.edges())
        {
            Scalar target_edge_length = (target_length[mesh_.vertex(edge, 0)] + target_length[mesh_.vertex(edge, 1)]) / 2.0f;

            if (mesh_.edge_length(edge) > (target_edge_length * 4.0f / 3.0f))
            {
                // Calculate the position of the vertex that will be positioned
                // at the point where the edge is split.
                //
                // We position the new vertex at the center of the split edge, so
                // calculate the center of the edge to get the vertex position.
                //
                Point new_vertex_position = (mesh_.position(mesh_.vertex(edge, 0)) + mesh_.position(mesh_.vertex(edge, 1))) / 2.0f;

                // Add the new vertex to the mesh at the calculated position.
                //
                auto new_vertex = mesh_.add_vertex(new_vertex_position);

                // Configure the properties of the new vertex.
                //
                target_length[new_vertex] = target_edge_length;
                normals[new_vertex] = mesh_.compute_vertex_normal(new_vertex);

                // Finally, split the edge at the new vertex.
                //
                mesh_.split(edge, new_vertex);

                // Since we have just splitted an edge, we start presuming that there may
                // be more edges to split.
                //
                finished = false;
            }
        }
    }

    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();

    if (i == 100) std::cerr << "split break\n";
}

bool MeshProcessing::check_collapse_ok(Mesh::Halfedge v0v1) {
    if (!mesh_.is_collapse_ok(v0v1))
        return false;
    Mesh::Vertex    v0(mesh_.from_vertex(v0v1));
    Mesh::Vertex    v1(mesh_.to_vertex(v0v1));
    Mesh::Vertex    vv, vv2, v_first;
    Mesh::Vertex_around_vertex_circulator vv_it, vv_end;
    vv_it = vv_end = mesh_.vertices(v0);
    v_first = *vv_it;
    do
    {
        vv = *vv_it;
        ++vv_it;
        if (vv_it != vv_end)
            vv2 = *vv_it;
        else
            vv2 = v_first;
        if ((vv == v1) || (vv2 == v1))
            continue;
        auto n0 = calc_triangle_norm(mesh_.position(v0), mesh_.position(vv), mesh_.position(vv2));
        auto n1 = calc_triangle_norm(mesh_.position(v1), mesh_.position(vv), mesh_.position(vv2));
        if ((norm(n1) < norm(n0)) && (norm(n1) < small_triangle))
            return false;
        n0 = normalize(n0);
        n1 = normalize(n1);
        auto dot01 = dot(n0, n1);
        if (dot01 < small_flip)
            return false;
    }
    while (vv_it != vv_end);

    return true;
}

void MeshProcessing::collapse_short_edges ()
{
    Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    bool            finished;
    int             i;

    for (finished = false, i = 0; !finished && i < 100; ++i)
    {
        // We assume that there is no more splitting left to do.
        //
        finished = true;

        // Lookup and collapse the edges.
        //
        for (auto halfedge : mesh_.halfedges())
        {
            auto edge = mesh_.edge(halfedge);
            
            // Do not collapse halfedges if their correspoding edges are already deleted.
            //
            if (mesh_.is_deleted(edge))
            {
                continue;
            }

            // Do not collapse edges that lie on the boundary.
            //
            if (mesh_.is_boundary(edge))
            {
                continue;
            }

            auto from_vertex = mesh_.from_vertex(halfedge);
            auto to_vertex = mesh_.to_vertex(halfedge);

            if (mesh_.is_boundary(from_vertex))
                continue;
            if (mesh_.is_boundary(to_vertex))
                continue;

            // Get the target edge length as the average of the corresponding vertices.
            //
            Scalar target_edge_length = (target_length[from_vertex] + target_length[to_vertex]) / 2.0f;

            // Do not collapse edges that are greater or equal than the 4/5 of the target edge length.
            //
            if (mesh_.edge_length(edge) >= (target_edge_length * 0.8))
            {
                continue;
            }

            auto opposite_halfedge = mesh_.opposite_halfedge(halfedge);


            // Start collapsing.
            //
            if (check_collapse_ok(halfedge) && check_collapse_ok(opposite_halfedge))
            {
                // If both halfedges are collapsible, then collapse the lower valence vertex into the higher one.
                //
                auto collapsed_halfedge = mesh_.valence(from_vertex) < mesh_.valence(to_vertex) ? halfedge : opposite_halfedge;
                mesh_.collapse(collapsed_halfedge);
            }
            else if (check_collapse_ok(halfedge))
            {
                mesh_.collapse(halfedge);
            }
            else if (check_collapse_ok(opposite_halfedge))
            {
                mesh_.collapse(opposite_halfedge);
            }
            else
            {
                // Do not collapse any of the two halfedges if neither of them is collapsible.
                //
                continue;
            }

            // Since we have just collapsed an edge, we start presuming that there may
            // be more edges to collapse.
            //
            finished = false;
        }
        mesh_.garbage_collection();
        mesh_.update_face_normals();
        mesh_.update_vertex_normals();
    }

    if (i == 100) std::cerr << "collapse break\n";
}

unsigned int MeshProcessing::get_ideal_valence(Mesh::Vertex vertex) const
{
    return mesh_.is_boundary(vertex) ? 4 : 6;
}

unsigned int MeshProcessing::calc_valence_deviation_squared(Mesh::Vertex vertex, bool flipped, int correction) const
{
    int const real_valence = !flipped ? mesh_.valence(vertex) : static_cast<int>(mesh_.valence(vertex)) + correction;
    int const ideal_valence = get_ideal_valence(vertex);

    int const valence_deviation = real_valence - ideal_valence;

    return abs(valence_deviation);
}

Vec3 MeshProcessing::calc_triangle_norm(Point a0, Point a1, Point a2) {
    auto edge01 = a1 - a0;
    auto edge02 = a2 - a0;
    auto edge12 = a2 - a1;
    edge01 = normalize(edge01);
    edge02 = normalize(edge02);
    edge12 = normalize(edge12);
    auto normal0 = cross(edge01, edge02);
    auto normal1 = cross(edge12, -edge01);
    auto normal2 = cross(-edge02, -edge12);
    if (norm(normal0) < norm(normal1)) {
        if (norm(normal0) < norm(normal2))
            return normal0;
        else
            return normal2;
    }
    else if (norm(normal1) < norm(normal2))
        return normal1;
    else
        return normal2;
}

bool MeshProcessing::check_flip_edge_ok(Mesh::Vertex up, Mesh::Vertex down, Mesh::Vertex left, Mesh::Vertex right) {
    auto norm_old0 = calc_triangle_norm(mesh_.position(up), mesh_.position(left), mesh_.position(down));
    auto norm_old1 = calc_triangle_norm(mesh_.position(up), mesh_.position(down), mesh_.position(right));
    auto norm_new0 = calc_triangle_norm(mesh_.position(up), mesh_.position(left), mesh_.position(right));
    auto norm_new1 = calc_triangle_norm(mesh_.position(down), mesh_.position(right), mesh_.position(left));
    if (norm(norm_new0) < small_triangle)
        return false;
    if (norm(norm_new1) < small_triangle)
        return false;
    norm_old0 = normalize(norm_old0);
    norm_old1 = normalize(norm_old1);
    norm_new0 = normalize(norm_new0);
    norm_new1 = normalize(norm_new1);
    if (dot(norm_new0, norm_new1) < 0)
        return false;
    if (dot(norm_old0, norm_new0) < small_flip)
        return false;
    if (dot(norm_old0, norm_new1) < small_flip)
        return false;
    if (dot(norm_old1, norm_new0) < small_flip)
        return false;
    if (dot(norm_old1, norm_new1) < small_flip)
        return false;
    return true;
}

void MeshProcessing::equalize_valences ()
{
    bool            finished;
    int             i;

    for (finished = false, i = 0; !finished && i < 100; ++i)
    {
        // We assume that there is no more flipping left to do.
        //
        finished = true;

        for (auto halfedge : mesh_.halfedges())
        {
            auto edge = mesh_.edge(halfedge);

            // Do not flip an edge if it is not suitable for it.
            //
            if (!mesh_.is_flip_ok(edge))
            {
                continue;
            }

            // Consider "Connectivity query" slide from SurfaceMeshTutorial.pdf
            // In this case down_vertex is v0 from the slide, up_vertex is v1 and so on.
            //
            auto down_vertex = mesh_.from_vertex(halfedge);
            auto up_vertex = mesh_.to_vertex(halfedge);
            auto left_vertex = mesh_.to_vertex(mesh_.next_halfedge(halfedge));
            auto right_vertex = mesh_.to_vertex((mesh_.next_halfedge(mesh_.opposite_halfedge(halfedge))));

            // Compute the initial local valence deviation.
            //
            auto const initial_local_valence_deviation
                    = calc_valence_deviation_squared(down_vertex)
                    + calc_valence_deviation_squared(up_vertex)
                    + calc_valence_deviation_squared(left_vertex)
                    + calc_valence_deviation_squared(right_vertex);

            // Compute the local valence after the potential flip.
            // Here, we sort of simulate the flip by accrodingly adjusting the valences.
            //
            auto const flipped_local_valence_deviation
                    = calc_valence_deviation_squared(down_vertex, true, -1)
                    + calc_valence_deviation_squared(up_vertex, true, -1)
                    + calc_valence_deviation_squared(left_vertex, true, +1)
                    + calc_valence_deviation_squared(right_vertex, true, +1);

            // Do not flip the edges if it will not bring any improvement in the local valence.
            //
            if (flipped_local_valence_deviation >= initial_local_valence_deviation)
            {
                continue;
            }

            // make sure that fliping edge doesn't flip faces
            if (!check_flip_edge_ok(up_vertex, down_vertex, left_vertex, right_vertex)) {
                continue;
            }

            mesh_.flip(edge);

            // Since we have just flipped an edge, we start presuming that there may
            // be more edges to flip.
            //
            finished = false;
        }
    }
    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();

    if (i == 100) std::cerr << "flip break\n";
}

bool MeshProcessing::check_relaxation_ok(Mesh::Vertex vertex, Vec3 updated) {
    auto v0 = mesh_.position(vertex);
    auto v1 = v0 + updated;
    Mesh::Vertex    vv, vv2, v_first;
    Mesh::Vertex_around_vertex_circulator vv_it, vv_end;
    vv_it = vv_end = mesh_.vertices(vertex);
    v_first = *vv_it;
    do
    {
        vv = *vv_it;
        ++vv_it;
        if (vv_it != vv_end)
            vv2 = *vv_it;
        else
            vv2 = v_first;
        auto n0 = calc_triangle_norm(v0, mesh_.position(vv), mesh_.position(vv2));
        auto n1 = calc_triangle_norm(v1, mesh_.position(vv), mesh_.position(vv2));
        if ((norm(n1) < norm(n0)) && (norm(n1) < small_triangle))
            return false;
        n0 = normalize(n0);
        n1 = normalize(n1);
        auto dot01 = dot(n0, n1);
        if (dot01 < small_flip)
            return false;
    }
    while (vv_it != vv_end);
    return true;
}

void MeshProcessing::tangential_relaxation ()
{
    Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");

    for (auto vertex : mesh_.vertices())
    {
        update[vertex] = Vec3(0,0,0);
        if (mesh_.is_boundary(vertex))
        {
            continue;
        }

        auto const vertex_valence = mesh_.valence(vertex);

        if (vertex_valence == 0)
        {
            continue;
        }

        // Calculate the uniform Laplacian smoothing vector.
        //
        Vec3 lu_smoothing_direction(0.0f, 0.0f, 0.0f);

        for (auto neighbor_vertex : mesh_.vertices(vertex))
        {
            lu_smoothing_direction += mesh_.position(neighbor_vertex) - mesh_.position(vertex);
        }
        lu_smoothing_direction /= vertex_valence;

        // Get the unit normal vector for the given surface at the current vertex point.
        // The normals in the corresponding vertex property should be already normalized.
        //
        Normal const & unit_normal = normals[vertex];

        // Compute the tangential component of the uniform laplacian smoothing vector
        // using the formula:
        // v_t = -n x (n x lu), where n is unit normal, x is the cross product and lu is our lu vector.
        //
        Vec3 const tangential_component = cross(-unit_normal, cross(unit_normal, lu_smoothing_direction));

        // Store the tangential component as a value in the update property of the current vertex.
        //
        update[vertex] = tangential_component;
    }

    // Smooth the vertices using their update properties.
    //
    int num_update = 0;
    for (auto vertex : mesh_.vertices())
    {
        if (mesh_.is_boundary(vertex))
        {
            continue;
        }

        if (mesh_.valence(vertex) == 0)
        {
            continue;
        }

        if (!check_relaxation_ok(vertex, update[vertex])) {
            continue;
        }

        mesh_.position(vertex) += update[vertex];
        num_update++;
    }
    cout << "update " << num_update << " vertices" << endl;
    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

// ========================================================================
// EXERCISE 1.1
// ========================================================================
void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate mean curvature using
    // the length of the uniform Laplacian approximation
    // Save your approximation in unicurvature vertex property of the mesh.
    // ------------- IMPLEMENT HERE ---------

    for (auto vertex: mesh_.vertices()) {
        if (mesh_.is_boundary(vertex))
            continue;
        int n = mesh_.valence(vertex);
        if (n == 0)
            continue;
        Vec3 v_sum(0,0,0);
        for (auto neighbor_vertex: mesh_.vertices(vertex))
            v_sum += mesh_.position(neighbor_vertex) - mesh_.position(vertex);
        v_sum /= n;
        v_unicurvature[vertex] = norm(v_sum) / 2;
    }
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


}


