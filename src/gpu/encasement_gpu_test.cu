#include "encasement_gpu.h"

int sub, iterations;

void _calculate(PTR<Object<PPoly2>> f, PTR<Object<PPoly2>> g) {

	Parameter::disable();

	auto t1 = std::chrono::high_resolution_clock::now();

	//find all the regions
	vector<Rectangle *> regions = get_regions(f, g, -2, -2, 2, 2, sub, iterations);

	Parameter::enable();

	vector<Root2PTR> roots;

	vector<PTR<Object<PPoly2>> > curves;
	curves.push_back(f);
	curves.push_back(g);

	//run encasement on all small regions and collect intersections
	for (int i = 0; i < regions.size(); i++) {
		PTR<Object<PV2>> ll = new InputPoint(regions[i]->box.x.lb(),
				regions[i]->box.y.lb());
		PTR<Object<PV2>> ur = new InputPoint(regions[i]->box.x.ub(),
				regions[i]->box.y.ub());

		Encasement * encasement = new Encasement(ll, ur, 27);

		vector<Root2PTR> lroots;

		encasement->init();

		encasement->calculate(curves);

		lroots = encasement->getRoots();

		roots.insert(roots.end(), lroots.begin(), lroots.end());

		printf("enc[%d] faces: %ld\n", i, encasement->faces.size());

	}

	for (int i = 0; i < roots.size(); i++) {
		PV2 p = roots[i]->getApprox();
		printf("r[%d]: [%.16g, %.16g] [%.16g, %.16g]\n", i, p.x.lb(), p.x.ub(),
				p.y.lb(), p.y.ub());
	}

	printf("roots size: %ld\n", roots.size());

	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << std::chrono::duration_cast < std::chrono::milliseconds
			> (t2 - t1).count() << " ms\n";

	Parameter::disable();

}

void calculate(vector<PTR<Object<PPoly2>> > curves) {

	for (int i = 0; i < curves.size(); i++) {
		for (int j = i + 1; j < curves.size(); j++) {
			_calculate(curves[i], curves[j]);
		}
	}

}

int main(int argc, char ** argv) {

	if (argc < 2) {
		printf("usage: ./encasement_gpu <file_name> [subdivision dimension = 100] [iterations = 1]\n");
		exit(1);
	}

	if(argc > 2) {
		sub = atoi(argv[2]);
	} else {
		sub = 100;
	}

	if(argc > 3) {
		iterations = atoi(argv[3]);
	} else {
		iterations = 1;
	}

	vector<PTR<Object<PPoly2>> > curves = loadCurves(argv[1]);

	try {
		calculate(curves);
	} catch (PrecisionException e) {
		printf("Precision Exception\n");
	} catch (SignException e) {
		printf("Sign Exception\n");
	} catch (IdentityException e) {
		printf("Identity Exception\n");
	}

}

