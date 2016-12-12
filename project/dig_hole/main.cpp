#include <iostream>

#include "HoleDigger.h"

using namespace std;

int main() {
  HoleDigger *digger = new HoleDigger();
  digger->read("geralt_remesh.obj");
  cout << "read" << endl;
  digger->sample(0.1);
  cout << "finish sampling" << endl;
  digger->digHole(0.07);
  cout << "finish digging" << endl;
  digger->save("sampled.obj");
  digger->save("sampled.off");
  return 0;
}