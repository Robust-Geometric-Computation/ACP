CFLAGS = -O3 -std=c++11 #-DEXACTRAT

COMPILE = g++ $(CFLAGS) -c

LINK = g++ $(CFLAGS)

LFLAGS = -lmpfr -lpthread -lgmp -lgmpxx

ACPH =	par.h ppar.h mpar.h mod.h acp.h

ACP = par.h par.cc ppar.h ppar.cc mpar.h mpar.cc mod.h mod.cc acp.h acp.cc	

all: test1 test2 test3 test4 test57a test57b test6 test8

test1:	pack.o mink.o polyhedron.o triangulate.o io.o acp.o
	$(LINK) pack.o mink.o polyhedron.o triangulate.o io.o \
	acp.o $(LFLAGS) -otest1

test2:	test2.C cdt.o delaunay.o mink.o polyhedron.o triangulate.o io.o acp.o
	$(LINK) test2.C cdt.o delaunay.o mink.o polyhedron.o triangulate.o io.o acp.o \
	$(LFLAGS) -otest2

test3:	test3.C cdt.o delaunay.o mink.o polyhedron.o triangulate.o io.o acp.o
	$(LINK) test3.C cdt.o delaunay.o mink.o polyhedron.o triangulate.o io.o acp.o \
	$(LFLAGS) -otest3

test4:	mesh.o delaunay.o polyhedron.o triangulate.o io.o acp.o
	$(LINK) mesh.o delaunay.o polyhedron.o triangulate.o io.o acp.o \
	$(LFLAGS) -otest4

test57a: freespace-graph_test1.o freespace-graph.o path3d.o hull.o \
	polyhedron.o freespace.o triangulate.o io.o acp.o mink.o
	$(LINK) freespace-graph_test1.o freespace-graph.o hull.o path3d.o \
	polyhedron.o freespace.o triangulate.o io.o acp.o mink.o \
	$(LFLAGS) -otest57a

test57b: freespace-graph_test2.o freespace-graph.o path3d.o hull.o \
	polyhedron.o freespace.o triangulate.o io.o acp.o mink.o
	$(LINK) freespace-graph_test2.o freespace-graph.o hull.o path3d.o \
	polyhedron.o freespace.o triangulate.o io.o acp.o mink.o \
	$(LFLAGS) -otest57b

test6: test6.C c4d.o acp.o angle.o poly.o root.o contact.o polym.o
	$(LINK) test6.C c4d.o acp.o angle.o poly.o root.o contact.o polym.o \
	$(LFLAGS) -otest6

test8:	test8.C cspace.o root.o poly.o mink.o polyhedron.o \
	triangulate.o io.o acp.o
	$(LINK) test8.C cspace.o root.o poly.o mink.o \
	 polyhedron.o triangulate.o io.o acp.o $(LFLAGS) -otest8

test8r2: csres.o cspace.o polym.o root.o poly.o mink.o \
	polyhedronr2.o triangulate.o io.o acp.o 
	$(LINK) csres.o cspace.o polym.o root.o poly.o mink.o \
	 polyhedronr2.o triangulate.o io.o acp.o $(LFLAGS) -otest8r2

test8r3: csres.o cspace.o polym.o root.o poly.o mink.o \
	polyhedronr3.o triangulate.o io.o acp.o 
	$(LINK) csres.o cspace.o polym.o root.o poly.o mink.o \
	 polyhedronr3.o triangulate.o io.o acp.o $(LFLAGS) -otest8r3

test8r4: csres.o cspace.o polym.o root.o poly.o mink.o \
	polyhedronr4.o triangulate.o io.o acp.o 
	$(LINK) csres.o cspace.o polym.o root.o poly.o mink.o \
	 polyhedronr4.o triangulate.o io.o acp.o $(LFLAGS) -otest8r4

csres.o: csres.C cspace.h polym.h root.h poly.h mink.h polyhedron.h io.h octree.h
	$(COMPILE) csres.C

polyhedronr2.o: polyhedronr.C polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) -DRESTEST2 polyhedronr.C -opolyhedronr2.o

polyhedronr3.o: polyhedronr.C polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) -DRESTEST3 polyhedronr.C -opolyhedronr3.o

polyhedronr4.o: polyhedronr.C polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) -DRESTEST4 polyhedronr.C -opolyhedronr4.o

polyhedron.o: polyhedron.C polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) polyhedron.C

io.o:	io.C io.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) io.C

triangulate.o: triangulate.C triangulate.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) triangulate.C

mink.o: mink.C mink.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) mink.C

hull.o: hull.C hull.h polyhedron.h pv.h $(ACPH)
	$(COMPILE) hull.C

pack.o:	pack.C pack.h mink.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) pack.C

cspace.o: cspace.C cspace.h root.h poly.h mink.h polyhedron.h io.h octree.h \
	pv.h $(ACPH)
	$(COMPILE) cspace.C

acp.o:	$(ACP)
	$(COMPILE) acp.cc

hull:	hull.o io.o triangulate.o polyhedron.o acp.o
	$(LINK) hull.o io.o triangulate.o polyhedron.o acp.o \
	$(LFLAGS) -o hull

delaunay.o: delaunay.C delaunay.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) delaunay.C

cdt.o:	cdt.C cdt.h delaunay.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) cdt.C

mesh.o:	mesh.C mesh.h heap.h delaunay.h polyhedron.h octree.h pv.h $(ACPH)
	$(COMPILE) mesh.C

poly.o:	poly.cc poly.h pv.h $(ACPH)
	$(COMPILE) poly.cc

root.o:	root.cc root.h poly.h pv.h $(ACPH)
	$(COMPILE) root.cc

polym.o: polym.cc polym.h poly.h pv.h acp.h par.h mod.h ppar.h mpar.h
	$(COMPILE) polym.cc

freespace-graph.o: freespace-graph.cc freespace-graph.h path3d.h hull.h \
	polyhedron.h octree.h acp.h geometry3d.h
	$(COMPILE) freespace-graph.cc

freespace-graph_test1.o: freespace-graph_test.cc freespace-graph.h path3d.h hull.h \
	polyhedron.h octree.h acp.h geometry3d.h
	$(COMPILE) -DTWO_ROOM_EXAMPLE freespace-graph_test.cc -ofreespace-graph_test1.o

freespace-graph_test2.o: freespace-graph_test.cc freespace-graph.h path3d.h hull.h \
	polyhedron.h octree.h acp.h geometry3d.h
	$(COMPILE) -DLATTICE_EXAMPLE freespace-graph_test.cc -ofreespace-graph_test2.o

path3d.o: path3d.cc path3d.h polyhedron.h octree.h acp.h geometry3d.h io.h
	$(COMPILE) path3d.cc

freespace.o: freespace.cc freespace.h hull.h polyhedron.h octree.h acp.h \
	geometry3d.h
	$(COMPILE) freespace.cc

c4d.o:	c4d.cc c4d.h angle.h pv.h root.h poly2.h poly.h acp.h
	$(COMPILE) c4d.cc

angle.o: angle.cc angle.h pv.h root.h poly2.h poly.h acp.h polym.h
	$(COMPILE) angle.cc

contact.o: contact.h contact.cc
	$(COMPILE) contact.cc
clean: 
	rm -f *.o *~ test1 test2 test3 test4 test57a test57b test6 test8 \
	test8r2 test8r3 test8r4
