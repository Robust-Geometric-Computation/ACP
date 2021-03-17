#include "encasement_gpu.h"

vector<PTR<Object<PPoly2>> > loadCurves(const char* fname) {

	vector<PTR<Object<PPoly2>> > curves;

	ifstream inFile(fname);
	string line;

	Parameter::enable();

	while (getline(inFile, line)) {

		istringstream buf(line);
		istream_iterator<std::string> beg(buf), end;

		vector<Parameter> coef;
		vector<int> pow;

		for (; beg != end;) {
			double c = stod((*beg));
			beg++;
			int i = stoi((*beg));
			beg++;
			int j = stoi((*beg));
			beg++;

			coef.push_back(Parameter::input(c));
			pow.push_back(i);
			pow.push_back(j);

		}

		curves.push_back(new Ellipse(PPoly2(coef.size(), &coef[0], &pow[0])));

	}

	Parameter::disable();

	return curves;
}

thrust::host_vector<i_pair<Number> > get_regions_gpu(PTR<Object<PPoly2>> f_in,
		PTR<Object<PPoly2>> g_in, double xl, double yl, double xu, double yu, int sub, int iterations) {

	//two device vectors, swapping buffers
	//first has the m intervals
	//subdivide so that second has floor(N / m) boxes for each of the m intervals (give iterators so no dead threads)
	//run the copy_if on that vector and keep the output iterator to know how many
	//swap buffers

	// input size
	int N = sub * sub;

	// define some types
	typedef thrust::device_vector<i_pair<Number> > d_vector_interval;
	typedef thrust::device_vector<int> d_vector_int;

	typedef d_vector_interval::iterator d_vector_interval_it;
	typedef d_vector_int::iterator d_vector_int_it;

	PPoly2 fp = f_in->get();
	PPoly2 fpx = fp.derX();
	PPoly2 fpy = fp.derY();
	PPoly2 fpxx = fpx.derX();
	PPoly2 fpxy = fpx.derY();
	PPoly2 fpyy = fpy.derY();

	PPoly2 gp = g_in->get();
	PPoly2 gpx = gp.derX();
	PPoly2 gpy = gp.derY();
	PPoly2 gpxx = gpx.derX();
	PPoly2 gpxy = gpx.derY();
	PPoly2 gpyy = gpy.derY();

	interval_gpu<Number> fcoef[fp.nt];
	interval_gpu<Number> fxcoef[fp.nt];
	interval_gpu<Number> fycoef[fp.nt];
	interval_gpu<Number> fxxcoef[fp.nt];
	interval_gpu<Number> fxycoef[fp.nt];
	interval_gpu<Number> fyycoef[fp.nt];

	interval_gpu<Number> gcoef[fp.nt];
	interval_gpu<Number> gxcoef[fp.nt];
	interval_gpu<Number> gycoef[fp.nt];
	interval_gpu<Number> gxxcoef[fp.nt];
	interval_gpu<Number> gxycoef[fp.nt];
	interval_gpu<Number> gyycoef[fp.nt];

	for (int i = 0; i < fp.nt; i++) {
		fcoef[i] = interval_gpu<T>(T(fp.a[i].lb()), T(fp.a[i].ub()));
	}
	for (int i = 0; i < fpx.nt; i++) {
		fxcoef[i] = interval_gpu<T>(T(fpx.a[i].lb()), T(fpx.a[i].ub()));
	}
	for (int i = 0; i < fpy.nt; i++) {
		fycoef[i] = interval_gpu<T>(T(fpy.a[i].lb()), T(fpy.a[i].ub()));
	}
	for (int i = 0; i < fpxx.nt; i++) {
		fxxcoef[i] = interval_gpu<T>(T(fpxx.a[i].lb()), T(fpxx.a[i].ub()));
	}
	for (int i = 0; i < fpxy.nt; i++) {
		fxycoef[i] = interval_gpu<T>(T(fpxy.a[i].lb()), T(fpxy.a[i].ub()));
	}
	for (int i = 0; i < fpyy.nt; i++) {
		fyycoef[i] = interval_gpu<T>(T(fpyy.a[i].lb()), T(fpyy.a[i].ub()));
	}

	for (int i = 0; i < gp.nt; i++) {
		gcoef[i] = interval_gpu<T>(T(gp.a[i].lb()), T(gp.a[i].ub()));
	}
	for (int i = 0; i < gpx.nt; i++) {
		gxcoef[i] = interval_gpu<T>(T(gpx.a[i].lb()), T(gpx.a[i].ub()));
	}
	for (int i = 0; i < gpy.nt; i++) {
		gycoef[i] = interval_gpu<T>(T(gpy.a[i].lb()), T(gpy.a[i].ub()));
	}
	for (int i = 0; i < gpxx.nt; i++) {
		gxxcoef[i] = interval_gpu<T>(T(gpxx.a[i].lb()), T(gpxx.a[i].ub()));
	}
	for (int i = 0; i < gpxy.nt; i++) {
		gxycoef[i] = interval_gpu<T>(T(gpxy.a[i].lb()), T(gpxy.a[i].ub()));
	}
	for (int i = 0; i < gpyy.nt; i++) {
		gyycoef[i] = interval_gpu<T>(T(gpyy.a[i].lb()), T(gpyy.a[i].ub()));
	}

	cuda_poly2<Number> hf(fcoef, fp.m, fp.nt);
	cuda_poly2<Number> hfx(fxcoef, fpx.m, fpx.nt);
	cuda_poly2<Number> hfy(fycoef, fpy.m, fpy.nt);
	cuda_poly2<Number> hfxx(fxxcoef, fpxx.m, fpxx.nt);
	cuda_poly2<Number> hfxy(fxycoef, fpxy.m, fpxy.nt);
	cuda_poly2<Number> hfyy(fyycoef, fpyy.m, fpyy.nt);

	cuda_poly2<Number> hg(gcoef, gp.m, gp.nt);
	cuda_poly2<Number> hgx(gxcoef, gpx.m, gpx.nt);
	cuda_poly2<Number> hgy(gycoef, gpy.m, gpy.nt);
	cuda_poly2<Number> hgxx(gxxcoef, gpxx.m, gpxx.nt);
	cuda_poly2<Number> hgxy(gxycoef, gpxy.m, gpxy.nt);
	cuda_poly2<Number> hgyy(gyycoef, gpyy.m, gpyy.nt);

	cuda_poly2<Number> * f = hf.copy_device();
	cuda_poly2<Number> * fx = hfx.copy_device();
	cuda_poly2<Number> * fy = hfy.copy_device();
	cuda_poly2<Number> * fxx = hfxx.copy_device();
	cuda_poly2<Number> * fxy = hfxy.copy_device();
	cuda_poly2<Number> * fyy = hfyy.copy_device();

	cuda_poly2<Number> * g = hg.copy_device();
	cuda_poly2<Number> * gx = hgx.copy_device();
	cuda_poly2<Number> * gy = hgy.copy_device();
	cuda_poly2<Number> * gxx = hgxx.copy_device();
	cuda_poly2<Number> * gxy = hgxy.copy_device();
	cuda_poly2<Number> * gyy = hgyy.copy_device();

	// allocate storage for array
	d_vector_int values(N);

	// initialize array to [0, 1, 2, ... ]
	thrust::sequence(values.begin(), values.end());

	// allocate output storage, here we conservatively assume all values will be copied
	d_vector_interval sub_div_intervals(values.size());
	d_vector_interval predicate_intervals(values.size());

	d_vector_interval_it predicate_intervals_end = predicate_intervals.begin()
			+ 1;
	// setup initial interval
	i_pair<Number> initial;
	initial.x = interval_gpu<Number>(xl, xu);
	initial.y = interval_gpu<Number>(yl, yu);

	sub_div_intervals[0] = initial;
	predicate_intervals[0] = initial;

	int num_intervals = 1;

	for (int index = 0; index < iterations; index++) {

		// subdivide the interval into subintervals

		int subs_per_thread = N / num_intervals;

		//subs_per_thread = max(4, subs_per_thread);

		cout << "subs_per: " << subs_per_thread << endl;

		int rows = int(sqrtf(subs_per_thread));
		int cols = subs_per_thread / rows;

		cout << "rows: " << rows << endl;
		cout << "cols: " << cols << endl;

		subs_per_thread = rows * cols;

		cout << "subs_per: " << subs_per_thread << endl;

		int total_threads = subs_per_thread * num_intervals;

		cout << "total_threads: " << total_threads << endl;

		struct subdivide<Number> subdivide_functor(total_threads,
				thrust::raw_pointer_cast(&predicate_intervals[0]),
				num_intervals);

		thrust::transform(values.begin(), values.begin() + total_threads,
				sub_div_intervals.begin(), subdivide_functor);

		// copy all intervals that ARE ambiguous
		//struct is_ambiguous<Number> ambiguous_functor(f, g);
		struct is_ambiguous_quad<Number> ambiguous_functor(f, fx, fy, fxx, fxy,
				fyy, g, gx, gy, gxx, gxy, gyy);

		predicate_intervals_end = thrust::copy_if(sub_div_intervals.begin(),
				sub_div_intervals.begin() + total_threads,
				predicate_intervals.begin(), ambiguous_functor);

		num_intervals = int(
				predicate_intervals_end - predicate_intervals.begin());

		//FOR TESTING

		struct function_value<Number> function_value_functor(f, g);

		d_vector_interval function_values(num_intervals);

		thrust::transform(predicate_intervals.begin(), predicate_intervals_end,
				function_values.begin(), function_value_functor);

		thrust::host_vector<i_pair<Number> > hfunction_values = function_values;
		thrust::host_vector<i_pair<Number> > hintervals(
				predicate_intervals.begin(), predicate_intervals_end);

		printf("removed intervals size: %d\n", num_intervals);

		//for(int i = 0; i < hintervals.size(); i++) {
		//    printf("f([%f, %f][%f, %f]) -> [%f, %f]\n", hintervals[i].x.lower(), hintervals[i].x.upper(), hintervals[i].y.lower(), hintervals[i].y.upper(), hfunction_values[i].x.lower(), hfunction_values[i].x.upper());
		//    printf("g([%f, %f][%f, %f]) -> [%f, %f]\n\n", hintervals[i].x.lower(), hintervals[i].x.upper(), hintervals[i].y.lower(), hintervals[i].y.upper(), hfunction_values[i].y.lower(), hfunction_values[i].y.upper());
		//}

	}

	return thrust::host_vector < i_pair<Number>
			> (predicate_intervals.begin(), predicate_intervals_end);
}

vector<Rectangle *> get_all_rects(PTR<Object<PPoly2>> f, PTR<Object<PPoly2>> g, double xl,
		double yl, double xu, double yu, int sub, int iterations) {

	thrust::host_vector<i_pair<Number> > hintervals = get_regions_gpu(f, g, xl,
			yl, xu, yu, sub, iterations);

	vector<Rectangle *> all_rects;

	for (int i = 0; i < hintervals.size(); i++) {
		all_rects.push_back(
				new Rectangle(hintervals[i].x.lower(), hintervals[i].y.lower(),
						hintervals[i].x.upper(), hintervals[i].y.upper()));
	}

	return all_rects;

}
vector<Rectangle *> get_regions(PTR<Object<PPoly2>> f, PTR<Object<PPoly2>> g, double xl,
		double yl, double xu, double yu, int sub, int iterations) {

	cudaGetLastError();
	cudaDeviceReset();
	cudaSetDevice(0);

	vector<Rectangle *> all_rects = get_all_rects(f, g, xl, yl, xu, yu, sub, iterations);

	vector<Rectangle *> bounding_boxes;

	map<int, set<Rectangle *> > m = connected_components(all_rects);

	for (map<int, set<Rectangle *> >::iterator it = m.begin(); it != m.end();
			it++) {

		set<Rectangle *> s = (*it).second;

		double rxl, ryl, rxu, ryu;
		bool first = true;

		for (set<Rectangle *>::iterator rectit = s.begin(); rectit != s.end();
				rectit++) {
			Rectangle * rect = *(rectit);
			if (first) {
				rxl = rect->box.x.lb();
				ryl = rect->box.y.lb();
				rxu = rect->box.x.ub();
				ryu = rect->box.y.ub();
				first = false;
			} else {
				rxl = fmin(rxl, rect->box.x.lb());
				ryl = fmin(ryl, rect->box.y.lb());
				rxu = fmax(rxu, rect->box.x.ub());
				ryu = fmax(ryu, rect->box.y.ub());
			}
		}

		bounding_boxes.push_back(new Rectangle(rxl, ryl, rxu, ryu));

	}

	return bounding_boxes;

}
