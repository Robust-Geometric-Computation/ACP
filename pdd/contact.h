#include <iostream>
using std::istream;
using std::ostream;
#include <vector>
using std::vector;

const size_t rrr[2][4] =
  { { 1191988316, 404730216, 779230947, 247091631 },
    { 1913983098, 150084830, 1194289623, 1109051806 } };
const size_t rrrr[2][4][3] =
  { { { 1452551476, 1650559630, 1525553500 },
      { 1529282427, 268402422, 137342870 },
      { 1297000885, 512898772, 657409826 },
      { 738066950, 1032600208, 738043542 } },
    { { 432735458, 1546955785, 1615662605 },
      { 1010520012, 621957185, 1381691171 },
      { 1069490946, 1893262845, 1045285240 },
      { 679597322, 1971193484, 89789909 } } };

class ALists;

class Contact {
public:
  int/*char*/ ns[2][4];
  int is[2][4][3];

  Contact (); // { ns[0][0] = ns[1][0] = 0; }

  Contact (int ns[2][4], int is[2][4][3]);

  Contact (const Contact &f, bool mf);

  int getN () const {
    return (ns[0][0] + ns[1][0] == 0 ? 0 :
            ns[0][1] + ns[1][1] == 0 ? 1 :
            ns[0][2] + ns[1][2] == 0 ? 2 :
            ns[0][3] + ns[1][3] == 0 ? 3 : 4);
  }

  istream &read (istream &in);
  istream &read (int n, istream &in);

  ostream &write (ostream &out) {
    out << getN() << " ";
    for (int i = 0; i < 4 && ns[0][i] + ns[1][i] > 0; i++) {
      out << " ";
      for (int fm = 0; fm <= 1; fm++) {
        out << "  " << ns[fm][i];
        for (int j = 0; j < ns[fm][i]; j++)
          out << " " << is[fm][i][j];
      }
    }
    return out;
  }

  class Compare {
  public:
    bool operator() (const Contact &a, const Contact &b) const {
      for (int fm = 0; fm <= 1; fm++)
	for (int i = 0; i < 4; i++)
	  if (a.ns[fm][i] != b.ns[fm][i])
	    return a.ns[fm][i] < b.ns[fm][i];
      for (int fm = 0; fm <= 1; fm++)
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < a.ns[fm][i]; j++)
            if (a.is[fm][i][j] != b.is[fm][i][j])
              return a.is[fm][i][j] < b.is[fm][i][j];
      return false;
    }
  };

  bool operator== (const Contact &b) const {
    Compare compare;
    return !compare(*this, b) && !compare(b, *this);
  }

  bool operator< (const Contact &b) const {
    Compare compare;
    return compare(*this, b);
  }

  size_t hashCode () const {
    size_t code = 0;
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++) {
	code += rrr[fm][i] * ns[fm][i];
	for (int j = 0; j < ns[fm][i]; j++)
	  code += rrrr[fm][i][j] * is[fm][i][j];
      }
    return code;
  }

  class Equal {
  public:
    bool operator() (const Contact &a, const Contact &b) const {
      return a == b;
    }

    size_t operator() (const Contact &a) const {
      return a.hashCode();
    }
  };

  // m[1][v] is original index of (contiguous) moving vertex v.
  class VMap {
  public:
    int m[2][12];
    VMap () {
      for (int fm = 0; fm <= 1; fm++)
        for (int v = 0; v < 12; v++)
          m[fm][v] = 0;
    }

    VMap inverse (int n[2]) {
      VMap inv;
      for (int fm = 0; fm <= 1; fm++)
	for (int i = 0; i < n[fm]; i++)
	  inv.m[fm][m[fm][i]] = i;
      return inv;
    }

    void randomize (int n) {
      for (int fm = 0; fm <= 1; fm++)
        for (int v = 0; v < 12; v++) {
	  int w;
	  do {
	    m[fm][v] = random() % n;
	    w = 0;
	    while (m[fm][w] != m[fm][v])
	      w++;
	  } while (w < v);
	}
    }
  };

  void reindex (VMap &vm);

  Contact (ALists &);

  Contact getRep (VMap &) const;

  Contact getKey (VMap &) const;

  int getSwaps () const;

  vector<Contact> getOtherKeys (vector<VMap> &vms);

  int compareFaces (int h, int i);

  void swapFaces (int i, int j);

  void normalize ();

  void replace (VMap &vm);

  static void readFactors ();

  vector<Contact> getFactors () const;
};

static Contact nullContact;
static Contact::VMap nullVMap;

inline Contact::Contact () { *this = nullContact; }

inline istream &operator>> (istream &in, Contact &c) { return c.read(in); }
inline ostream &operator<< (ostream &out, Contact &c) { return c.write(out); }

void initSwapBits();
