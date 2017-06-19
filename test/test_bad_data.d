//
// TestBadData.d
//

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


bool MySearchCallback(ValueType id, void* arg)
{
	import std.stdio : writeln;
	writeln("Hit data rect ", id);
	return true; // keep going
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

	alias MyTree = RTree!(ValueType, CoordType, 2, float);
	auto tree = new MyTree();

	int i, nhits;
	writeln("number of rectangles is ", rectVector.length);

	for(i=0; i<rectVector.length(); i++)
  {
		writeln(i, ": ", rectVector[i]);
    tree.Insert(rectVector[i].min, rectVector[i].max, i); // Note, all values including zero are fine in this version
  }

	auto search_rect = Rect(6, 4, 10, 6);
	nhits = tree.Search(search_rect.min, search_rect.max, &MySearchCallback, null);

	writeln("Search resulted in ", nhits, " hits");

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

