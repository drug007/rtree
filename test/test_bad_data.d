#!/usr/bin/env dub
/+ dub.sdl:
    name        "test_bad_data"
    targetType  "executable"

    dependency "rtree" version="*" path=".."
+/

import rtree;

alias ValueType = int;
alias CoordType = long;

struct Rect
{
	this(CoordType a_minX, CoordType a_minY, CoordType a_maxX, CoordType a_maxY)
	{
		min[0] = a_minX;
		min[1] = a_minY;

		max[0] = a_maxX;
		max[1] = a_maxY;
	}


	CoordType[2] min;
	CoordType[2] max;
}

int main(string[] args)
{
	import std.stdio : writeln;

	if ( args.length < 2 ) {
		writeln("Usage: ", args[0], " inFile");
		return -1;
	}

	import std.container.array : Array;
	alias RectVector = Array!Rect;
	RectVector rectVector;

	// read the data
	{
		import std.stdio : File;
		auto f = File(args[1]);
		foreach (e; f.byRecord!(CoordType, CoordType, CoordType, CoordType)("%s %s %s %s"))
		{
			rectVector ~= Rect(e[0], e[1], e[0]+e[2], e[1]+e[3]);
		}
	}

	import std.experimental.allocator.mallocator : Mallocator;
	auto tree = RTree!(Mallocator, ValueType, CoordType, 2, float).make();

	writeln("number of rectangles is ", rectVector.length);

	foreach(i; 0..cast(int)rectVector.length)
	{
		writeln(i, ": ", rectVector[i]);
		tree.insert(rectVector[i].min, rectVector[i].max, i); // Note, all values including zero are fine in this version
	}

	auto search_rect = Rect(6, 4, 10, 6);
	auto result = tree.search(search_rect.min, search_rect.max);

	writeln("Search resulted in ", result.length, " hits");
	writeln("Hit data rect ", result[]);

	//// Iterator test
	//int itIndex = 0;
	//MyTree::Iterator it;
	//for( tree.GetFirst(it);
	//		 !tree.IsNull(it);
	//		 tree.GetNext(it) )
	//{
	//	int value = tree.GetAt(it);

	//	CoordType boundsMin[2] = {0,0};
	//	CoordType boundsMax[2] = {0,0};
	//	it.GetBounds(boundsMin, boundsMax);
	//	cout << "it[" << itIndex++ << "] " << value << " = (" << boundsMin[0] << "," << boundsMin[1] << "," << boundsMax[0] << "," << boundsMax[1] << ")\n";
	//}

	//// Iterator test, alternate syntax
	//itIndex = 0;
	//tree.GetFirst(it);
	//while( !it.IsNull() )
	//{
	//	CoordType value = *it;
	//	++it;
	//	cout << "it[" << itIndex++ << "] " << value << "\n";
	//}

	return 0;
}

