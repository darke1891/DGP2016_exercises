#include "HoleDigger.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <queue>

using namespace std;
using namespace surface_mesh;

typedef pair<double, Surface_mesh::Vertex> pdv;
typedef pair<int, double> pid;

HoleDigger::HoleDigger()
{
}


HoleDigger::~HoleDigger()
{
}

void HoleDigger::read(string name) {
  mesh.read(name);
}

void HoleDigger::save(std::string name) {
  mesh.write(name);
}

void HoleDigger::octreeSample(double distance, vector<Surface_mesh::Vertex> &input,
  vector<Surface_mesh::Vertex> &to_pick) {
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
    auto pos = mesh.position(vertex);
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
    vector<Surface_mesh::Vertex> subRange[8];
    for (int i = 0; i < 8; i++)
      subRange[i].clear();
    for (auto vertex : input) {
      auto pos = mesh.position(vertex);
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

Surface_mesh::Vertex HoleDigger::findNearest(Point center, vector<Surface_mesh::Vertex> &range) {
  float les = -1;
  Surface_mesh::Vertex result;
  for (auto vertex : range) {
    float dis = norm(mesh.position(vertex) - center);
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


void HoleDigger::sample(double distance) {
  srand(unsigned(time(0)));
  vector<int> picked;
  vector<Surface_mesh::Vertex> to_pick;
  sample_points.clear();
  picked.clear();
  to_pick.clear();
  dis_to_sample.clear();
  vector<Surface_mesh::Vertex> all_vertices;
  all_vertices.clear();
  for (auto vertex : mesh.vertices())
    all_vertices.push_back(vertex);
  octreeSample(distance, all_vertices, to_pick);
  cout << to_pick.size() << endl;

  for (auto vertex: mesh.vertices()) {
    to_pick.push_back(vertex);
    picked.push_back(0);
    dis_to_sample.push_back(make_pair(distance * 2, -1));
  }
  //random_shuffle(to_pick.begin(), to_pick.end());

  for (auto vertex: to_pick) {
    int pick = vertex.idx();
    if (picked[pick] == 1)
      continue;
    sample_points.push_back(pick);
    priority_queue<pdv> to_travel;
    while (!to_travel.empty())
      to_travel.pop();

    to_travel.push(pdv(0, vertex));
    while (!to_travel.empty()) {
      pdv this_v = to_travel.top();
      to_travel.pop();
      auto this_vertex = this_v.second;
      double dist = -this_v.first;
      if (dist < dis_to_sample[this_vertex.idx()].first) {
        picked[this_vertex.idx()] = 1;
        dis_to_sample[this_vertex.idx()] = make_pair(dist, vertex.idx());
        for (auto halfedge: mesh.halfedges(this_vertex)) {
            auto edge = mesh.edge(halfedge);
          auto next_vertex = mesh.to_vertex(halfedge);
          double next_dist = dist + mesh.edge_length(edge);
          if (next_dist < distance)
            if (next_dist < dis_to_sample[next_vertex.idx()].first)
              to_travel.push(pdv(- next_dist, next_vertex));
        }
      }
    }
  }
  cout << sample_points.size() << endl;
}

void HoleDigger::digHole(double distance) {
  for (auto face : mesh.faces()) {
    double dis = 0;
    int range_index = -1;
    bool boundary = false;
    for (auto vertex : mesh.vertices(face)) {
      dis += distance - dis_to_sample[vertex.idx()].first;
      if (range_index == -1)
        range_index = dis_to_sample[vertex.idx()].second;
      else if (dis_to_sample[vertex.idx()].second != -1)
        if (dis_to_sample[vertex.idx()].second != range_index)
          boundary = true;
    }
    if ((dis > 0) && (!boundary)) {
      mesh.delete_face(face);
    }
  }
  mesh.garbage_collection();
  for (auto vertex : mesh.vertices())
    if (mesh.valence(vertex) == 0)
      mesh.delete_vertex(vertex);
  mesh.garbage_collection();
  mesh.update_face_normals();
  mesh.update_vertex_normals();
}