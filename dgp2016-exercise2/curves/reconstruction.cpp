#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>
#include <cmath>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
    PointsRenderer render_points = PointsRenderer();
    SegmentsRenderer render_segments = SegmentsRenderer();

    double epsilon;
    MatMxN points;
    MatMxN points_3d_render;
    Vec2 tangent;
    SegmentsRenderer::Segments segments;


    // ============================================================================
    // Exercise 3 : fill the function below (see PDF for instructions)
    // To test your implementation, use the SPACE key
    // ============================================================================
    /*
    //If we can calculate Integral(k(s))ds Annalitically :
    double primitiveOfk(double s){

        //k(s) = 1
        //return s;

        //k(s) = s
        return s*s/2;

        //k(s) = s² - 2.19
        //return (s*s/3 - 2.19) * s;

    }
    Vec2 tangentAt(double s){
        double pks = primitiveOfk(s);
        double sini = sin(pks);
        double cosi = cos(pks);
        return Vec2(cosi * tangent(0) - sini * tangent(1),sini * tangent(0) + cosi * tangent(1));
    }
    Vec2 deltaTangentAt(double s){
        return (tangentAt(s+epsilon)+tangentAt(s))/2 * epsilon;
    }
    */

    //If we can't have Integral(k(s))ds Annalitically :
    double k(double s) {

        //k(s) = 1
        //return 1;

        //k(s) = s
        return s*s;

        //k(s) = s² - 2.19
        //return s*s - 2.19;
    }

    Vec2 deltaTangentAt(double s){
        double deltaks = (k(s)+k(s+epsilon))/2*epsilon;
        double sind = sin(deltaks);
        double cosd = cos(deltaks);
        Vec2 prevT = tangent;
        tangent = Vec2(
            tangent(0) * cosd - tangent(1) * sind,
            tangent(1) * cosd + tangent(0) * sind
        );
        return (prevT + tangent)/2*epsilon;
    }

    void reconstructCurveStep() {
        //every time you press space, a new point of the curve is added, so the number of points increases
        points.conservativeResize(2, points.cols() + 1);
        int n = points.cols() - 1;
        double s = n * epsilon;
        Vec2 deltaTangente = deltaTangentAt(s);
        points.col(n) = deltaTangente+points.col(n-1);

    }
    // ============================================================================
    // END OF Exercise 2 (do not thouch the rest of the code except to uncomment
    // the part where you can change the function)
    // ============================================================================

    void render() {
        // Prepare the render points
        points_3d_render = MatMxN::Zero(3, points.cols());
        points_3d_render.block(0, 0, 2, points.cols()) = points * 1e-1;

        // Rebuild the segments
        segments.clear();

        for (int i = 1; i < points_3d_render.cols (); ++i) {
            segments.push_back({ points_3d_render.col(i), points_3d_render.col(i-1)});
        }

        render_points.init_data(points_3d_render);
        render_segments.init_data(segments);
    }

    MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
        points = MatMxN(2, 1);

        // First, second, third curve
        //points.col(0) = Vec2 (0, -1); // first
        //points.col (0) = Vec2 (0, 0.0); // second
        points.col (0) = Vec2 (0, 0.0); // third

        tangent = Vec2 (1.0, 0.0);

        epsilon = 0.1;

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render ();
    }

    bool key_callback(int key, int scancode, int action, int mods) override {
        TrackballWindow::key_callback(key, scancode, action, mods);
        if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
        {
            reconstructCurveStep();
        }

        render ();
        return true;
    }
};


int main(int argc, char** argv)
{
    MainWindow window(argc, argv);
    return window.run();
}
