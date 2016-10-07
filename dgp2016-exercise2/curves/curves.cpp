#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
  PointsRenderer render_points = PointsRenderer ();
  SegmentsRenderer render_segments = SegmentsRenderer ();

  MatMxN points;
  MatMxN points_3d_render;
  int num_points;
  SegmentsRenderer::Segments segments;

// ============================================================================
// Exercise 2 : fill the 2 functions below (see PDF for instructions)
// To test your implementation, use the S key for laplacian smoothing and the
// C key for the osculating circle.
// Hint : try to play with epsilon
// ============================================================================
  // time step for smoothing
  double epsilon = 0.1;

  double getCurveLength(OpenGP::MatMxN const & curve)
  {
	  double length = 0.0;
	  
	  for (auto i = 0; i < num_points; ++i)
	  {
		  // Compute the indices of the points which take part in the current length calculation.
		  //
		  // Due to the fact that the given curve is closed,
		  // the indices can wrap at the edges of the points buffer.
		  //
		  auto const currentpoint_index = i;
		  auto const nextpoint_index = (i + 1) % num_points;
		  
		  // Calculate the length of the current curve segment using the norm of the vector which
		  // goes from the current point to the next one.
		  //
		  length += (curve.col(nextpoint_index) - curve.col(currentpoint_index)).norm();
	  }
	  
	  return length;
  }

  void laplacianSmoothing()
  {
	  // Initialize a matrix for the processed points with zeros.
	  // Every column in the resulting matrix represents the coordinates of a certain point.
	  //
	  Eigen::Index const coords_count = 2;
	  OpenGP::MatMxN smoothed_points = OpenGP::MatMxN::Zero(coords_count, num_points);
	  
	  
	  //
	  // Process the initial points into their corresponding "smoothed" versions.
	  //
	  // This will effectively fill the container for the smoothed points.
	  //
	  
	  for (auto i = 0; i < num_points; ++i)
	  {
		  // Calculate the indices of both the previous and the next points.
		  //
		  // Due to the fact that the given curve is closed,
		  // the indices can wrap at the edges of the points buffer.
		  //
		  auto const prevpoint_index = (i + num_points - 1) % num_points;
		  auto const nextpoint_index = (i + 1) % num_points;
		  
		  // Calculate the center of the line that connects the previous point with the next one.
		  //
		  auto const line_center = (points.col(prevpoint_index) + points.col(nextpoint_index)) / 2.0;
		  
		  // Create a vector that points from the current point to the line center calculated above.
		  //
		  // The vector is scaled by the epsilon value.
		  //
		  auto const currentpoint_index = i;
		  auto const currentpoint_direction = epsilon * (line_center - points.col(currentpoint_index));
		  
		  // Create a smoothed version of the current point.
		  //
		  smoothed_points.col(currentpoint_index) = points.col(currentpoint_index) + currentpoint_direction;
	  }
	  
	  
	  //
	  // Update the initial points with their corresponding "smoothed" versions while preserving the initial curve length.
	  //
	  
	  // Calculate the lengths of both the initial curve and the smoothed one.
	  //
	  auto const initialcurve_length = getCurveLength(points);
	  auto const smoothedcurve_length = getCurveLength(smoothed_points);
		  
	  for (auto i = 0; i < num_points; ++i)
	  {  
		  points.col(i) = smoothed_points.col(i) * initialcurve_length / smoothedcurve_length;
	  }
  }

  void osculatingCircle() {
    // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
  }

// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code)
// ============================================================================

  void generateRandomizedClosedPolyline() {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5*3e-2);

    Vec2 center(3e-2, 2e-3);
    const double radius = 0.3;

    points = MatMxN::Zero(2, num_points);
    for (int i = 0; i < num_points; ++i)
    {
      double frac = static_cast<double>(i) / static_cast<double>(num_points);
      points(0, i) = center(0) + radius * cos (2. * M_PI * frac) + distribution(generator);
      points(1, i) = center(1) + radius * sin (2. * M_PI * frac) + distribution(generator);
    }
  }

  void render () {

    // Prepare the render points
    points_3d_render = MatMxN::Zero(3, points.cols());
    points_3d_render.block(0, 0, 2, points.cols()) = points;

    // Rebuild the segments
    segments.clear();
    for (int i = 0; i < points_3d_render.cols(); ++i) {
      segments.push_back({ points_3d_render.col(i), points_3d_render.col((i+1) % points_3d_render.cols()) });
    }
    render_points.init_data(points_3d_render);
    render_segments.init_data(segments);
  }

  MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
    num_points = 50;
    generateRandomizedClosedPolyline();

    this->scene.add(render_points);
    this->scene.add(render_segments);

    render();
  }

  bool key_callback(int key, int scancode, int action, int mods) override {
    TrackballWindow::key_callback(key, scancode, action, mods);
    if (key == GLFW_KEY_S && action == GLFW_RELEASE)
    {
      laplacianSmoothing();     
    }
    else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
    {
      osculatingCircle();
    }

    render();
    return true;
  }
};


int main(int argc, char** argv)
{
  MainWindow window(argc, argv);
  return window.run();
}
