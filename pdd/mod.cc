/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

unsigned int inverse (unsigned int a, unsigned int n)
{
  long t = 0, nt = 1, r = n, nr = a;
  while (nr != 0) {
    long q = r / nr, pt = t, pr = r;
    t = nt;
    nt = pt - q * nt;
    r = nr;
    nr = pr - q * nr;
  }
  return t < 0 ? t + n : t;
}

const int Modder::eShift = 53;
const int Modder::eMax = 1024 - eShift;
const int Modder::eMin = -1073 - eShift;
unsigned int Mods::ps[NMods];
Modder Mods::modder[NMods];
pthread_mutex_t Mods::mutex = PTHREAD_MUTEX_INITIALIZER;

Mods initializeMods(0, 0, 0);

Mods::Mods (const MValue& l, const MValue& u)
{
  if (l != u)
    for (int i = 0; i < NMods; i++) mod[i] = Mod(random() % ps[i], ps[i]);
  else {
    mpz_t m, pz, mp, two, ez, ep;
    mpz_init(m);
    mpz_init(pz);
    mpz_init(mp);
    mpz_init_set_si(two, 2);
    mpz_init(ep);
    int e = mpfr_get_z_2exp(m, l.m);
    if (e > 0) {
      cerr << "internal error in Mods::Mods" << endl;
      exit(0);
    }
    mpz_init_set_si(ez, -e);
    for (int i = 0; i < NMods; i++) {
      mpz_set_ui(pz, ps[i]);
      mpz_mod(mp, m, pz);
      mpz_powm_sec(ep, two, ez, pz);
      long mpl = mpz_get_si(mp), epl = mpz_get_si(ep);
      mod[i] = Mod(mpl, ps[i]) / Mod(epl, ps[i]);
    }
    mpz_clear(m);
    mpz_clear(pz);
    mpz_clear(mp);
    mpz_clear(two);
    mpz_clear(ez);
    mpz_clear(ep);
  }
}

void Mods::changePrime (unsigned int p)
{
  pthread_mutex_lock(&mutex);
  cout << "changePrime " << p << endl;
  for (int i = 0; i < NMods; i++)
    if (p == ps[i]) {
      ps[i] = random32bitPrime();
      cout << "changeprime " << ps[i] << endl;
      modder[i] = Modder(ps[i]);
      pthread_mutex_unlock(&mutex);
      return;
    }
  cerr << "internal error in Mods::changePrime" << endl;
  exit(0);
}

unsigned long mulmod (unsigned long a, unsigned long b, unsigned long mod)
{
  unsigned long x = 0, y = a % mod;
  while (b > 0) {
    if (b%2 == 1) {
      x = (x + y)%mod;
    }
    y = (y * 2)%mod;
    b /= 2;
  }
  return x%mod;
}

unsigned long modulo (unsigned long base, unsigned long exponent,
		      unsigned long mod)
{
  unsigned long x = 1, y = base;
  while (exponent > 0) {
    if (exponent % 2 == 1)
      x = (x*y)%mod;
    y = (y * y) % mod;
    exponent = exponent/2;
  }
  return x%mod;
}

int Miller (unsigned long p, int iteration)
{
  int i;
  unsigned long s;
  if (p < 2)
    return 0;
  if (p != 2 && p % 2 == 0)
    return 0;
  s = p - 1;
  while (s%2 == 0)
    s /= 2;
  for (i = 0; i < iteration; i++) {
    unsigned long a = random() % (p - 1) + 1, temp = s;
    unsigned long mod = modulo(a, temp, p);
    while (temp != p - 1 && mod != 1 && mod != p - 1) {
      mod = mulmod(mod, mod, p);
      temp *= 2;
    }
    if (mod != p - 1 && temp%2 == 0)
      return 0;
  }
  return 1;
}

bool isPrime (unsigned int p) { return Miller(p, 40); }

unsigned int random32bitPrime ()
{
  while (true) {
    unsigned int p = random();
    p |= (1 << 31);
    if (!isPrime(p))
      continue;
    return p;
  }
}

