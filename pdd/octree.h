#ifndef OCTREE
#define OCTREE

#include <vector>
#include <iostream>

using namespace std;

template <class T>
class Octree {
 public:
  class Ocelt {
  public:
    Ocelt (T d) : d(d), label(0u) {}
    bool compatible (const Ocelt &e) const { return (label & e.label) == 0u; }
    
    void expandBBox (double s) {
      double bb[6];
      d->getBBox(bb);
      for (int i = 0; i < 3; ++i) {
	bb[2*i] -= s;
	bb[2*i+1] += s;
      }
    }

    T d;
    unsigned int label;
  };

  typedef vector<Ocelt> Ocelts;
  
  typedef unsigned int ID;
  
  bool bboxOverlap (const double *a, const double *b, double s = 0.0) const
  {
    for (int i = 0; i < 3; ++i)
      if (a[2*i+1] + s < b[2*i] || b[2*i+1] + s < a[2*i])
	return false;
    return true;
  }
  Octree () : elts(elts), n(0), v(0), l(0), r(0) {}

  Octree (const Ocelts &elts) : elts(elts), n(elts.size()), v(0), l(0), r(0) {}
  
  Octree (int n, int c, double v, Octree *l, Octree *r)
    : n(n), c(c), v(v), l(l), r(r) {}

  ~Octree () { if (l) { delete l; delete r; } }

  static Octree * octree (const vector<T> &data, const double *bbox,
			  double s = 0.0, int nmax = 40, int dmax = 20) {
    Ocelts edata;
    for (int i = 0; i < data.size(); ++i)
      edata.push_back(data[i]);
    if (s != 0.0)
      for (int i = 0; i < data.size(); ++i)
	edata[i].expandBBox(s);
    return octree(edata, bbox, nmax, dmax, 0u);
  }
  
  static Octree * octree (const Ocelts &data, const double *bbox, int nmax,
			  int dmax, ID d) {
    int n = data.size();
    if (n <= nmax || dmax == 0)
      return new Octree(data);
    int c = bbox[1] - bbox[0] < bbox[3] - bbox[2] ? 2 : 1;
    if (bbox[2*c-1] - bbox[2*c-2] < bbox[5] - bbox[4])
      c = 3;
    double v = 0.5*(bbox[2*c-1] + bbox[2*c-2]);
    Ocelts data1, data2;
    partition(data, c, v, d, data1, data2);
    if (data1.size() == n && data2.size() == n)
      return new Octree(data);
    double bbox1[6], bbox2[6];
    for (int i = 0; i < 6; ++i)
      bbox1[i] = bbox2[i] = bbox[i];
    bbox1[2*c-1] = bbox2[2*c-2] = v;
    Octree *l = octree(data1, bbox1, nmax, dmax - 1, d + 1u),
      *r = octree(data2, bbox2, nmax, dmax - 1, d + 1u);
    return new Octree(n, c, v, l, r);
  }

  static void partition (const Ocelts &data, int c, double v, unsigned int d,
			 Ocelts &data1, Ocelts &data2)
  {
    int j = 2*(c - 1);
    for (int i = 0; i < data.size(); ++i) {
      const Ocelt &a = data[i];
      double bba[6];
      a.d->getBBox(bba);
      if (bba[j+1] < v)
	data1.push_back(a);
      else if (v < bba[j])
	data2.push_back(a);
      else {
	data1.push_back(a);
	Ocelt b(a);
	b.label += 1u << d;
	data2.push_back(b);
      }
    }
  }

  void list (vector<T> &res, unsigned int d = 0u, unsigned int label = 0u) const {
    if (l) {
      l->list(res, d + 1u, label);
      unsigned int rlabel = label + (1u << d);
      r->list(res, d + 1u, rlabel);
    }
    else 
      for (int i = 0; i < elts.size(); ++i) {
	const Ocelt &f = elts[i];
	if ((f.label & label) == 0u)
	  res.push_back(f.d);
      }
  }

  void pairs (vector<pair<T, T>> &res, double s = 0.0) const {
    if (l) {
      l->pairs(res, s);
      r->pairs(res, s);
    }
    else
      for (int i = 0; i < elts.size(); ++i) {
	double bbi[6];
	elts[i].d->getBBox(bbi);
	for (int j = i + 1; j < elts.size(); ++j)
	  if (elts[i].compatible(elts[j])) {
	    double bbj[6];
	    elts[j].d->getBBox(bbj);
	    if (bboxOverlap(bbi, bbj, s))
	      res.push_back(pair<T, T>(elts[i].d, elts[j].d));
	  }
      }
  }
  
  void find (double *eb, vector<T> &res, double s = 0.0, unsigned int d = 0u, 
	     unsigned int label = 0u) const {
    if (l) {
      if (eb[2*c-1] + s < v)
	l->find(eb, res, s, d + 1u, label);
      else if (v + s < eb[2*c-2])
	r->find(eb, res, s, d + 1u, label);
      else {
	l->find(eb, res, s, d + 1u, label);
	unsigned int rlabel = label + (1u << d);
	r->find(eb, res, s, d + 1u, rlabel);
      }
    }
    else 
      for (int i = 0; i < elts.size(); ++i) {
	const Ocelt &f = elts[i];
	double bbf[6];
	f.d->getBBox(bbf);
	if ((f.label & label) == 0u && bboxOverlap(eb, bbf, s))
	  res.push_back(f.d);
      }
  }

  void insert (const Ocelt &e, int nmax = 40, int dmax = 20, unsigned int d = 0u) {
    if (l) {
      double bb[6];
      e.d->getBBox(bb);
      if (bb[2*c-1] < v)
	l->insert(e, nmax, dmax - 1, d + 1u);
      else if (v < bb[2*c-2])
	r->insert(e, nmax, dmax - 1, d + 1u);
      else {
	l->insert(e, nmax, dmax - 1, d + 1u);
	Ocelt er(e);
	er.label += 1u << d;
	r->insert(er, nmax, dmax - 1, d + 1u);
      }
    }
    else {
      for (int i = 0; i < elts.size(); ++i)
	if (elts[i].d == e.d)
	  return;
      elts.push_back(e);
    }
  }

  void remove (const Ocelt &e) {
    if (l) {
      double bb[6];
      e.d->getBBox(bb);
      if (bb[2*c-1] < v)
	l->remove(e);
      else if (v < bb[2*c-2])
	r->remove(e);
      else {
	l->remove(e);
	r->remove(e);
      }
    }
    else {
      int i = 0;
      while (i < elts.size())
	if (elts[i].d == e.d) {
	  elts[i] = *elts.rbegin();
	  elts.pop_back();
	}
	else
	  ++i;
    }
  }
  
  void expand (int m, vector<Octree *> &res) {
    ID nmax = n/m;
    vector<Octree *> st;
    st.push_back(this);
    while (!st.empty()) {
      Octree *o = *st.rbegin();
      st.pop_back();
      if (!o->l || o->n <= nmax)
	res.push_back(o);
      else {
	st.push_back(o->l);
	st.push_back(o->r);
      }
    }
  }
  
  void stats (vector<int> &ne) const {
    if (l) {
      l->stats(ne);
      r->stats(ne);
    }
    else
      ne.push_back(n);
  }

  void stats () const {
    vector<int> ne;
    stats(ne);
    sort(ne.begin(), ne.end());
    int n = ne.size();
    double m = 0.0;
    for (int i = 0; i < n; ++i)
      m += ne[i];
    m /= n;
    cerr << "min " << ne[0] << "; max " << ne[n-1] << "; med = " << ne[n/2]
	 << "; mean " << m << endl;
  }

  Ocelts elts;
  int n, c;
  double v;
  Octree *l, *r;	  
};

#endif
