#pragma once

template <class S>
void MakeSet(S* x) {
  x->size() = 1;
  x->parent() = x;
}

template <class S>
S* Find(S* x) {
  S* root = x;
  while (root->parent() != root) root = root->parent();
  while (x->parent() != root) {
    S* parent = x->parent();
    x->parent() = root;
    x = parent;
  }
  return root;
}

template <class S>
void Union(S* x, S* y) {
  x = Find(x);
  y = Find(y);
  if (x == y) return;
  if (x->size() < y->size()) {
    x->parent() = y;
    y->size() += x->size();
  } else {
    y->parent() = x;
    x->size() += y->size();
  }
}
