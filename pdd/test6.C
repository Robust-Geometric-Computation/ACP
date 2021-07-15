#include "c4d.h"

int main (int argc, char *argv[]) {
  long iLong;
  long long iLongLong;
  //cout << sizeof(iLong) << " " << sizeof(iLongLong) << endl;
  initSwapBits();
  if (argc > 1) {
    int seed = atoi(argv[1]);
    cout << "seed " << seed << endl;
    srandom(seed);
  }
  double t = getTime();
  C4D c4d;
  c4d.findIdentitiesT();
  cerr << "cpu time: " << getTime() - t << endl;
  acp::primitiveReport();
  return 0;
}
