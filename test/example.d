#!/usr/bin/env dub
/+ dub.sdl:
    name        "example"
    targetType  "executable"

    dependency "rtree" version="*" path=".."
+/

import rtree;

//
// This is a direct port to D of the C++ version of the RTree test program.
//

alias ValueType = int;

struct Rect
{
	this(int a_minX, int a_minY, int a_maxX, int a_maxY)
	{
		min[0] = a_minX;
		min[1] = a_minY;

		max[0] = a_maxX;
		max[1] = a_maxY;
	}


	int[2] min;
	int[2] max;
}

Rect[] rects =
[
	Rect(0, 0, 2, 2), // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
	Rect(5, 5, 7, 7),
	Rect(8, 5, 9, 6),
	Rect(7, 1, 9, 2),
];

Rect search_rect = Rect(6, 4, 10, 6); // search will find above rects that this one overlaps


bool MySearchCallback(ValueType id)
{
	import std.stdio : writeln;
	writeln("Hit data rect ",  id);
	return true; // keep going
}


int test(Tree)(ref Tree tree)
{
	int i, nhits;
	auto nrects = rects.length;
	import std.stdio : printf;
	printf("nrects = %ld\n", nrects);

	for(i=0; i<nrects; i++)
	{
		tree.insert(rects[i].min, rects[i].max, i); // Note, all values including zero are fine in this version
	}

	import std.functional : toDelegate;
	auto result = tree.search(search_rect.min, search_rect.max);

	printf("Search resulted in [");
	foreach(e; result[])
		printf("%d, ", e);
	printf("]\n");

	import std.algorithm : equal;
	int[2] etalon = [1, 2];
	assert(result[].equal(etalon[]));

	// Iterator test
	//int itIndex = 0;
	//MyTree::Iterator it;
	//for( tree.GetFirst(it);
	//!tree.IsNull(it);
	//tree.GetNext(it) )
	//{
	//int value = tree.GetAt(it);

	//int boundsMin[2] = {0,0};
	//int boundsMax[2] = {0,0};
	//it.GetBounds(boundsMin, boundsMax);
	//cout << "it[" << itIndex++ << "] " << value << " = (" << boundsMin[0] << "," << boundsMin[1] << "," << boundsMax[0] << "," << boundsMax[1] << ")\n";
	//}

	//// Iterator test, alternate syntax
	//itIndex = 0;
	//tree.GetFirst(it);
	//while( !it.IsNull() )
	//{
	//int value = *it;
	//++it;
	//cout << "it[" << itIndex++ << "] " << value << "\n";
	//}

	return 0;

	// Output:
	//
	// nrects = 4
	// Hit data rect 1
	// Hit data rect 2
	// Search resulted in 2 hits
	// it[0] 0 = (0,0,2,2)
	// it[1] 1 = (5,5,7,7)
	// it[2] 2 = (8,5,9,6)
	// it[3] 3 = (7,1,9,2)
	// it[0] 0
	// it[1] 1
	// it[2] 2
	// it[3] 3
}

int main()
{
	() {
		import std.experimental.allocator.gc_allocator : GCAllocator;
		auto tree = RTree!(GCAllocator, ValueType, int, 2, float).make();
		test(tree);
	} ();

	() @nogc {
		import std.experimental.allocator.mallocator : Mallocator;
		auto tree = RTree!(Mallocator, ValueType, int, 2, float).make();
		test(tree);
	} ();

	() @nogc {
		import std.experimental.allocator.building_blocks.region;
		alias R = InSituRegion!(1 * 1024, 16);
		R r;
		auto tree = RTree!(R, ValueType, int, 2, float).make(r);
		test(tree);
	} ();
	
	return 0;
}
