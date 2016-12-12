#pragma once
#include <vector>
#include <string>
#include <set>
#include <surface_mesh\Surface_mesh.h>
#include "Type.h"

class HoleDigger
{
public:
  HoleDigger();
  ~HoleDigger();

  void read(std::string name);
  void save(std::string name);
  void sample(double distance);
  void digHole(double distance);
  void octreeSample(double distance, std::vector<surface_mesh::Surface_mesh::Vertex> &input,
                    std::vector<surface_mesh::Surface_mesh::Vertex> &to_pick);
  surface_mesh::Surface_mesh::Vertex findNearest(surface_mesh::Point center,
                                                 std::vector<surface_mesh::Surface_mesh::Vertex> &range);

  std::vector<int> sample_points;
  std::vector<std::pair<double, int> > dis_to_sample;

  surface_mesh::Surface_mesh mesh;
};

