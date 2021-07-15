#ifndef HEAP
#define HEAP

template<typename T>
class Heap {
 public:
  unsigned int n, nmax;
  T *elts;

  Heap () : n(0u), nmax(1u), elts(new T [2]) {}

  ~Heap () { delete [] elts; }

  Heap (const Heap &h) : n(h.n), nmax(h.nmax), elts(new T [h.nmax+1]) {
    for (unsigned int i = 0u; i < nmax + 1u; ++i)
      elts[i] = h.elts[i];
  }

  Heap & operator= (const Heap &h) {
    n = h.n;
    nmax = h.nmax;
    delete [] elts;
    elts = new T [h.nmax+1];
    for (unsigned int i = 0u; i < nmax + 1u; ++i)
      elts[i] = h.elts[i];
    return *this;
  }
  
  bool empty () const { return n == 0u; }

  T min () { return elts[1]; }

  void swap (unsigned int i, unsigned int j) {
    T temp = elts[i];
    elts[i] = elts[j];
    elts[j] = temp;
  }
    
  void expand () {
    nmax *= 2u;
    T *nelts = new T [nmax+1u];
    for (unsigned int i = 1u; i <= n; ++i)
      nelts[i] = elts[i];
    delete [] elts;
    elts = nelts;
  }

  void insert (T x) {
    if (n == nmax)
      expand();
    ++n;
    elts[n] = x;
    unsigned int i = n;
    while (i > 1u && elts[i] < elts[i/2u]) {
      swap(i, i/2u);
      i = i/2u;
    }
  }

  void removeMin () {
    elts[1u] = elts[n];
    --n;
    down(1u);
  }

  void down (unsigned int i) {
    while (i <= n/2u) {
      if (elts[i] < elts[2u*i] && 
	  (2u*i + 1u > n || elts[i] < elts[2u*i+1u]))
	return;
      if (2u*i + 1u > n || elts[2u*i] < elts[2u*i+1u]) {
	swap(i, 2*i);
	i = 2*i;
      }
      else {
	swap(i, 2*i+1);
	i = 2u*i + 1;
      }
    }
  }

  Heap (unsigned int n, T *x)
    : n(n), nmax(n), elts(new T [n+1u]) {
    for (unsigned int i = 0u; i < n; ++i)
      elts[i+1u] = x[i];
    for (unsigned int i = n/2u; i > 0u; --i)
      down(i);
  }
};

#endif
