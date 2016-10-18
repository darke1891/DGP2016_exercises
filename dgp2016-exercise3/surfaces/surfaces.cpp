//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#include "viewer.h"
#include <iostream>

using namespace surface_mesh;

void computeValence(Surface_mesh * mesh)
{
    Surface_mesh::Vertex_property<unsigned int> vertex_valence = mesh->vertex_property<unsigned int>("v:valence", 0);

    for (auto vertex : mesh->vertices())
    {
        vertex_valence[vertex] = mesh->valence(vertex);
    }
}

void computeNormalsWithConstantWeights(Surface_mesh * mesh)
{
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_cste_weights_n = mesh->vertex_property<Point>("v:cste_weights_n", default_normal);

    for (auto vertex : mesh->vertices())
    {
        Normal vertex_normal(0.0, 0.0, 0.0);

        for (auto face : mesh->faces(vertex))
        {
            vertex_normal += mesh->compute_face_normal(face);
        }
        v_cste_weights_n[vertex] = vertex_normal.normalize();
    }
}

void computeNormalsByAreaWeights(Surface_mesh * mesh)
{
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_area_weights_n = mesh->vertex_property<Point>("v:area_weight_n", default_normal);

    for (auto vertex : mesh->vertices())
    {
        Normal vertex_normal(0.0, 0.0, 0.0);

        auto initial_incident_vertex = mesh->vertices(vertex);
        auto current_incident_vertex = initial_incident_vertex;

        do
        {
            // Since we need two incident vertices to operate with, save
            // the current vertex and then iterate it to the next one.
            //
            auto previous_incident_vertex = current_incident_vertex;
            ++current_incident_vertex;

            // Get the positions of the current face vertices,
            // i.e. the position of the current vertex and the positions
            // of the two incident vertices that together form a triangular face.
            //
            auto const p0 = mesh->position(vertex);
            auto const p1 = mesh->position(*previous_incident_vertex);
            auto const p2 = mesh->position(*current_incident_vertex);

            // Get the vectors pointing from the current vertex to the incident vertices.
            //
            auto const a = p1 - p0;
            auto const b = p2 - p0;

            auto ab_cross_product = cross(a, b);


            // Compute the area of the incident face.
            //
            // Since the norm of a cross product of two vectors is equal to
            // the area of a parallelogram formed by these two vectors,
            // we divide this norm by 2 to get the area of a triangle,
            // the area of the incident triangular face.
            //
            auto const incident_face_area = norm(ab_cross_product) / 2.0;
            auto const incident_face_normal = ab_cross_product.normalize();

            vertex_normal += incident_face_area * incident_face_normal;
        }
        while (current_incident_vertex != initial_incident_vertex);

        v_area_weights_n[vertex] = vertex_normal.normalize();
    }
}

void computeNormalsWithAngleWeights(Surface_mesh * mesh)
{
    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_angle_weights_n = mesh->vertex_property<Point>("v:angle_weight_n", default_normal);

    for (auto vertex : mesh->vertices())
    {
        Normal vertex_normal(0.0, 0.0, 0.0);

        auto initial_incident_vertex = mesh->vertices(vertex);
        auto current_incident_vertex = initial_incident_vertex;

        do
        {
            // Since we need two incident vertices to operate with, save
            // the current vertex and then iterate it to the next one.
            //
            auto previous_incident_vertex = current_incident_vertex;
            ++current_incident_vertex;

            // Get the positions of the current face vertices,
            // i.e. the position of the current vertex and the positions
            // of the two incident vertices that together form a triangular face.
            //
            auto const p0 = mesh->position(vertex);
            auto const p1 = mesh->position(*previous_incident_vertex);
            auto const p2 = mesh->position(*current_incident_vertex);

            // Get the vectors pointing from the current vertex to the incident vertices.
            //
            auto const a = p1 - p0;
            auto const b = p2 - p0;

            auto const incident_face_angle = acos(dot(a, b) / norm(a) / norm(b));
            auto const incident_face_normal = cross(a, b).normalize();

            vertex_normal += incident_face_angle * incident_face_normal;
        }
        while (current_incident_vertex != initial_incident_vertex);

        v_angle_weights_n[vertex] = vertex_normal.normalize();
    }
}

// #############################################################################
int main(int /* argc */, char ** /* argv */) {
    try {
        nanogui::init();
        {
            // Load the Mesh
            Surface_mesh mesh;
            if (!mesh.read("../data/bunny.off")) {
                std::cerr << "Mesh not found, exiting." << std::endl;
                return -1;
            }

            computeValence(&mesh);
            computeNormalsWithConstantWeights(&mesh);
            computeNormalsByAreaWeights(&mesh);
            computeNormalsWithAngleWeights(&mesh);

            nanogui::ref<Viewer> app = new Viewer(&mesh);
            app->drawAll();
            app->setVisible(true);
            nanogui::mainloop();
        }

        nanogui::shutdown();
    } catch (const std::runtime_error &e) {
        std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(_WIN32)
        MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
        std::cerr << error_msg << endl;
#endif
        return -1;
    }

    return 0;
}
