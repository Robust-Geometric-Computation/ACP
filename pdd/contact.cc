#include "contact.h"
#include <assert.h>
#include <map>
using std::map;
#include <algorithm>
#include <fstream>
#include <unordered_map>
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;
using std::unordered_map;

Contact::Contact (int ns[2][4], int is[2][4][3]) {
  *this = nullContact;
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++) {
      this->ns[fm][i] = ns[fm][i];
      for (int j = 0; j < ns[fm][i]; j++)
        this->is[fm][i][j] = is[fm][i][j];
    }
}

vector<Contact> keys;
vector<vector<Contact> > kfactorss;
vector<vector<Contact> > kfactorss2;
unordered_map<Contact, int, Contact::Equal, Contact::Equal> kindex;

void Contact::readFactors () {
  ifstream kin;
  // kin.open("key-class3.txt");
  kin.open("min-class3.txt");
  ifstream fin;
  // fin.open("key-classes3.txt");
  fin.open("eqsub-classes3.txt");
  Contact key;
  while (kin >> key) {
    kindex[key] = keys.size();
    keys.push_back(key);
    int n = 0;
    fin >> n;
    vector<Contact> factors;
    for (int i = 0; i < n; i++) {
      int code;
      fin >> code;
      Contact factor;
      if (code != -99 && code != -999)
        factor.read(code, fin);
      else {
        factor.ns[0][0] = code;
        factor.ns[1][0] = kfactorss2.size();
        vector<Contact> factors2;
        int nn;
        fin >> nn;
        for (int j = 0; j < nn; j++) {
          Contact factor2;
          fin >> factor2;
          factors2.push_back(factor2);
        }
        kfactorss2.push_back(factors2);
      }
      factors.push_back(factor);
    }
    kfactorss.push_back(factors);
  }
}

vector<Contact> Contact::getFactors () const {
  VMap vm;
  // Contact key = getKey(vm);
  Contact key = getRep(vm);
  assert(kindex.find(key) != kindex.end());
  int index = kindex[key];
  assert(index < kfactorss.size());
  vector<Contact> &kfactors = kfactorss[index];
  vector<Contact> factors;
  for (int i = 0; i < kfactors.size(); i++) {
    Contact kfactor = kfactors[i];
    if (kfactor.ns[0][0] + kfactor.ns[1][0] == 0)
      continue;
    if (kfactor.ns[0][0] != -99 && kfactor.ns[0][0] != -999) {
      kfactor.replace(vm);
      kfactor.normalize();
      factors.push_back(kfactor);
    }
    else {
      assert(kfactor.ns[1][0] < kfactorss2.size());
      const vector<Contact> &factors2 = kfactorss2[kfactor.ns[1][0]];
      Contact factor = factors2[0];
      factor.replace(vm);
      factor.normalize();
      for (int i2 = 1; i2 < factors2.size(); i2++) {
        Contact factor2 = factors2[i2];
        factor2.replace(vm);
        factor2.normalize();
        if (factor2 < factor)
          factor = factor2;
      }
      factors.push_back(factor);
    }
  }
  return factors;
}

istream &Contact::read (istream &in) {
  *this = nullContact;
  int n = 0;
  in >> n;
  for (int i = 0; i < n; i++)
    for (int fm = 0; fm <= 1; fm++) {
      in >> ns[fm][i];
      for (int j = 0; j < ns[fm][i]; j++)
        in >> is[fm][i][j];
    }
  return in;
}

istream &Contact::read (int n, istream &in) {
  *this = nullContact;
  for (int i = 0; i < n; i++)
    for (int fm = 0; fm <= 1; fm++) {
      in >> ns[fm][i];
      for (int j = 0; j < ns[fm][i]; j++)
        in >> is[fm][i][j];
    }
  return in;
}

void Contact::reindex (VMap &vm) {
  static char old2new[1000000]; // WARNING!!!! ????  Fixed size array!!
  for (int fm = 0; fm <= 1; fm++) {
    int touched[16];
    int ntouched = 0;
    int n = 0;

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++) {
	int old = is[fm][i][j];
        if (old2new[old] == 0) {
          old2new[old] = n + 1;
	  touched[ntouched++] = old;
          vm.m[fm][n++] = old;
        }
        is[fm][i][j] = old2new[old] - 1;
      }

    for (int i = 0; i < ntouched; i++)
      old2new[touched[i]] = 0;
  }
}

#ifdef BLEEN
void Contact::reindex (VMap &vm) {
  for (int fm = 0; fm <= 1; fm++) {
    map<int, int> old2new;
    int n = 0;
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++) {
        if (old2new.find(is[fm][i][j]) == old2new.end()) {
          old2new[is[fm][i][j]] = n;
          vm.m[fm][n++] = is[fm][i][j];
        }
        is[fm][i][j] = old2new[is[fm][i][j]];
      }
  }
}
#endif

typedef int SWAPS[2];
namespace contact_swaps {
//   1
int n0012 = 1, sw0012[][2] = {{0, 1}};
//   01
int n0112 = 1, sw0112[][2] = {{1, 2}};
//   11
int n0001 = 5, sw0001[][2] = {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}};
//   001
int n0122 = 1, sw0122[][2] = {{2, 3}};
//   101
int n0011 = 3, sw0011[][2] = {{0, 1}, {2, 3}, {0, 1}};
//   011
int n0111 = 5, sw0111[][2] = {{1, 2}, {2, 3}, {1, 2}, {2, 3}, {1, 2}};
//   111
int n0000 = 23, sw0000[][2] = {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {1, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {1, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {0, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}};

int nswapss[] = { 0, n0012, n0112, n0001, n0122, n0011, n0111, n0000 };
SWAPS *swapss[] = { 0, sw0012, sw0112, sw0001, sw0122, sw0011, sw0111, sw0000 };

// swapBits[i][j][v] is the result of swapping bits 3-i and 3-j in v.
// 0 <= i < j < 4, 0 <= v < 16
int swapBits[4][4][16];
};
using namespace contact_swaps;

void initSwapBits () {
  for (int i = 0; i < 4; i++) {
    int bi = 1 << (3-i);
    for (int j = i+1; j < 4; j++) {
      int bj = 1 << (3-j);
      for (int v = 0; v < 16; v++) {
        char w = v;
        if ((v & bi) != 0)
          w |= bj;
        else
          w &= ~bj;
        if ((v & bj) != 0)
          w |= bi;
        else
          w &= ~bi;
        swapBits[i][j][v] = w;
      }
    }
  }

  assert(swapBits[0][1][0] == 0);
  assert(swapBits[0][1][4] == 8);
  assert(swapBits[0][1][8] == 4);
  assert(swapBits[0][1][12] == 12);
}

class Perm {
public:
  int p[4];
  
  Perm () {
    for (int i = 0; i < 4; i++)
      p[i] = i;
  }

  void swap (int i, int j) {
    int temp = p[i];
    p[i] = p[j];
    p[j] = temp;
  }

  void print (ostream &out) {
    out << p[0] << " " << p[1] << " " << p[2] << " " << p[3];
  }
};

class ALists {
public:
  int n[2];
  char a[2][12];

  // c must have contiguous vertices
  ALists (Contact &c) {
    for (int fm = 0; fm <= 1; fm++) {
      for (int i = 0; i < 12; i++)
        a[fm][i] = 0;
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < c.ns[fm][i]; j++)
          set(fm, c.is[fm][i][j], i);
      for (n[fm] = 12; n[fm] > 0 && a[fm][n[fm]-1] == 0; n[fm]--);
    }
  }

  void set (int fm, int v, int i) { a[fm][v] |= (1 << (4 - 1 - i)); }
  bool get (int fm, int v, int i) { return (a[fm][v] & (1 << (4 - 1 - i))) != 0; }

  void sort (Contact::VMap &vm) {
    sort(a[0], vm.m[0], n[0]);
    sort(a[1], vm.m[1], n[1]);
  }
      
  static void sort (char a[12], int m[12], int n) {
    for (int i = 1; i < n; i++) {
      int c = a[i];
      int d = m[i];
      int j = i;
      while (--j >= 0 && a[j] < c) {
        a[j+1] = a[j];
        m[j+1] = m[j];
      }
      a[j+1] = c;
      m[j+1] = d;
    }
  }

  void swap (int i, int j) {
    int *swap = swapBits[i][j];
    for (int fm = 0; fm <= 1; fm++)
      for (int v = 0; v < n[fm]; v++)
        a[fm][v] = swap[a[fm][v]];
  }

#ifdef BLEEN
  void permute (Perm &perm) {
    for (int k = 0; k < nF + nM; k++) {
      int ak = 0;
      for (int i = 0; i < 4; i++)
	if ((a[k] & (1 << (4 - 1 - i))) != 0)
	  ak |= (1 << (4 - 1 - perm.p[i]));
      a[k] = ak;
    }
    sort();
  }
#endif

  bool operator< (const ALists &that) const {
    for (int fm = 0; fm <= 1; fm++) {
      for (int i = 0; i < n[fm]; i++)
        if (a[fm][i] != that.a[fm][i])
          return a[fm][i] > that.a[fm][i];
    }
    return false;
  }

#ifdef BLEEN
  ALists minimum (int nSwaps, SWAPS *swaps, Perm &pMin) {
    ALists alMin = *this;
    Perm perm;
    for (int i = 0; i < nSwaps; i++) {
      swap(swaps[i][0], swaps[i][1]);
      perm.swap(swaps[i][0], swaps[i][1]);
      if (*this < alMin) {
        alMin = *this;
        pMin = perm;
        // cout << " swapped ";
      }
    }
    return alMin;
  }

  void getPerms (Graph &graph, int nSwaps, SWAPS *swaps);

  void print (ostream &out) {
    for (int i = 0; i < nF + nM; i++) {
      if (i != 0)
        out << " ";
      for (int j = 0; j < n; j++)
        if ((a[i] & (1 << (4 - j - 1))) != 0)
          out << j;
    }
  }
#endif

  class Compare {
  public:
    bool operator() (const ALists &a, const ALists &b) {
      if (a.n[0] != b.n[0])
        return a.n[0] < b.n[0];
      if (a.n[1] != b.n[1])
        return a.n[1] < b.n[1];
      return a < b;
    }
  };
};

Contact::Contact (ALists &lists) {
  *this = nullContact;
  for (int fm = 0; fm <= 1; fm++)
    for (int v = 0; v < lists.n[fm]; v++)
      for (int i = 0; i < 4; i++)
        if (lists.get(fm, v, i))
          is[fm][i][ns[fm][i]++] = v;
}

int Contact::getSwaps () const {
  int eqNext = 0;
  int bit = 1;
  int n = getN();
  for (int i = 1; i < n; i++) {
    if (ns[0][i-1] == ns[0][i] && ns[1][i-1] == ns[1][i])
      eqNext |= bit;
    bit *= 2;
  }
  return eqNext;
}

Contact Contact::getKey (Contact::VMap &vm) const {
  Contact contig = *this;
  contig.reindex(vm);
  ALists aLists(contig);
  aLists.sort(vm);
  return Contact(aLists);
}

Contact Contact::getRep (Contact::VMap &vm) const {
  int iSwaps = getSwaps();
  int nSwaps = nswapss[iSwaps];
  SWAPS *swaps = swapss[iSwaps];

  Contact contig = *this;
  contig.reindex(vm);

  ALists aLists(contig);
  aLists.sort(vm); // ??????
  ALists aListsMin = aLists;
  VMap vmMin = vm;

  for (int iSwap = 0; iSwap < nSwaps; iSwap++) {
    int i = swaps[iSwap][0];
    int j = swaps[iSwap][1];
    
    aLists.swap(i, j);
    aLists.sort(vm);

    if (aLists < aListsMin) {
      aListsMin = aLists;
      vmMin = vm;
    }
  }

  vm = vmMin;
  return Contact(aListsMin);
}

vector<Contact> Contact::getOtherKeys (vector<VMap> &vms) {
  int iSwaps = getSwaps();
  int nSwaps = nswapss[iSwaps];
  SWAPS *swaps = swapss[iSwaps];

  Contact contig = *this;
  VMap vm;
  contig.reindex(vm);
  assert(contig == *this);

  ALists aLists(contig);
  vector<Contact> keys;

  for (int iSwap = 0; iSwap < nSwaps; iSwap++) {
    int i = swaps[iSwap][0];
    int j = swaps[iSwap][1];
    
    aLists.swap(i, j);
    aLists.sort(vm);

    Contact key(aLists);
    if (key == *this)
      continue;

    int k = 0;
    for (k = 0; k < keys.size() && !(key == keys[k]); k++);
    if (k < keys.size())
      continue;

    {
      int n = key.getN();
      for (int j = 1; j < n; j++) {
	int i = j-1;
	assert(key.ns[0][i] + key.ns[1][i] < key.ns[0][j] + key.ns[1][j] ||
	       (key.ns[0][i] + key.ns[1][i] == key.ns[0][j] + key.ns[1][j] &&
		key.ns[0][i] <= key.ns[0][j]));
      }
    }	       

    keys.push_back(key);
    vms.push_back(vm.inverse(aLists.n));
    VMap vm0 = vms[0];
    vm0 = vms[0];
    assert(keys.size() == vms.size());
  }

  return keys;
}

int Contact::compareFaces (int h, int i) {
  int cmp; 
  if ((cmp = ns[0][h] + ns[1][h] - ns[0][i] - ns[1][i]) != 0)
    return cmp;
  if ((cmp = ns[0][h] - ns[0][i]) != 0)
    return cmp;
  for (int fm = 0; fm <= 1; fm++)
    for (int j = 0; j < ns[fm][h]; j++)
      if ((cmp = is[fm][h][j] - is[fm][i][j]) != 0)
        return cmp;
  return 0;
}

inline void swap (char &a, char &b) {
  int temp = a;
  a = b;
  b = temp;
}

inline void swap (int &a, int &b) {
  int temp = a;
  a = b;
  b = temp;
}

void swap (int a[3], int b[3]) {
  for (int i = 0; i < 3; i++)
    swap(a[i], b[i]);
}

void Contact::swapFaces (int i, int j) {
  swap(ns[0][i], ns[0][j]);
  swap(ns[1][i], ns[1][j]);
  swap(is[0][i], is[0][j]);
  swap(is[1][i], is[1][j]);
}

void Contact::normalize () {
  int n = getN();
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < n; i++)
      std::sort(is[fm][i], is[fm][i] + ns[fm][i]);
  for (int h = 0; h < n; h++)
    for (int i = h+1; i < n && ns[0][i] + ns[1][i] > 0; i++)
      if (compareFaces(h, i) > 0)
        swapFaces(h, i);
}

void Contact::replace (VMap &vm) {
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
	is[fm][i][j] = vm.m[fm][is[fm][i][j]];
}

Contact key;

int main1 (int argc, char *argv[]) {
  srandom(4479799);
  for (int i = 0; i < 32; i++)
    cout << random() << endl;

  initSwapBits();

  ifstream reps;
  reps.open("min-class3.txt");
  ifstream rfactors;
  rfactors.open("eqsub-classes3.txt");

  ofstream keys;
  keys.open("key-class3.txt");
  ofstream kfactors;
  kfactors.open("key-classes3.txt");

  int nrep = 0;
  int nkey = 0;
  Contact rep;
  while (reps >> rep) {
    nrep++;
    int nfactors;
    rfactors >> nfactors;
    vector<vector<Contact> > factorss;
    vector<int> codes;
    for (int i = 0; i < nfactors; i++) {
      vector<Contact> factors;
      int n;
      rfactors >> n;
      codes.push_back(n);
      if (n != -99 && n != -999) {
        Contact factor;
        factor.read(n, rfactors);
        factors.push_back(factor);
      }
      else {
        int nn;
        rfactors >> nn;
        for (int j = 0; j < nn; j++) {
          Contact factor;
          rfactors >> factor;
          factors.push_back(factor);
        }
      }
      factorss.push_back(factors);
    }

    nkey++;
    keys << rep << endl;
    kfactors << nfactors << " ";
    for (int i = 0; i < nfactors; i++) {
      vector<Contact> &factors = factorss[i];
      if (factors.size() == 1)
        kfactors << "   " << factors[0];
      else {
        kfactors << " " << codes[i] << " " << factors.size() << " ";
        for (int j = 0; j < factors.size(); j++)
          kfactors << "   " << factors[j];
      }
    }
    kfactors << endl;

    vector<Contact::VMap> vms;
    vector<Contact> otherKeys = rep.getOtherKeys(vms);
    assert(otherKeys.size() == vms.size());
    for (int i = 0; i < otherKeys.size(); i++) {
      Contact::VMap vm = vms[i];
      key = otherKeys[i];
      Contact::VMap vm2 = vms[i];
      nkey++;
      keys << key << endl;
      kfactors << nfactors << " ";
      for (int j = 0; j < nfactors; j++) {
        vector<Contact> &factors = factorss[j];
        if (factors.size() == 1) {
          Contact factor = factors[0];
          factor.replace(vm);
          factor.normalize();
          kfactors << "   " << factor;
        }
        else {
          kfactors << " " << codes[j] << " " << factors.size() << " ";
          for (int k = 0; k < factors.size(); k++) {
            Contact factor = factors[k];
            factor.replace(vm);
            factor.normalize();
            kfactors << "   " << factor;
          }
        }
      }
      kfactors << endl;
    }
  }

  return 0;
}

int hist[10];

int main2 (int argc, char *argv[]) {
  int a[] = { 3, 1, 4, 1, 5, 9, 2, 6 };
  std::sort(a, a+8);

  Contact::readFactors();

  for (int count = 0; count < 1000; count++) {
    Contact::VMap vm;
    vm.randomize(100);
    Contact c = keys[random() % keys.size()];
    c.replace(vm);
    c.normalize();
    vector<Contact> factors = c.getFactors();
    int nfactors = factors.size();
    hist[nfactors]++;
  }

  cout << "hist";
  for (int i = 0; i < 10; i++)
    cout << " " << hist[i];
  cout << endl;    

  return 0;
}
