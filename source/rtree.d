module rtree;

/// \class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree!(Object, float, 3) myTree;
///
/// This is ported to D C++ version by Yariv Barkan (https://github.com/nushoin/RTree)
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than (void*).sizeof and simple type
/// ELEMTYPE Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3
/// ELEMTYPEREAL Type of element that allows fractional and large values such as float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed size allocator for efficiency.
///        Instead of using a callback function for returned results, I recommend use efficient pre-sized, grow-only memory
///        array similar to MFC CArray or STL Vector for returning search query result.
///
class RTree(DATATYPE, ELEMTYPE, alias NUMDIMS,
	ELEMTYPEREAL, alias int MAXNODES = 8, alias int MINNODES = MAXNODES / 2)
{
public:

	// Precomputed volumes of the unit spheres for the first few dimensions
	enum float[] UnitSphereVolume = [
		0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
		4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
		5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
		3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
		1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
		0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
		0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
	];

	/// Unit sphere constant for required number of dimensions
	enum unitSphereVolume = cast(ELEMTYPEREAL) UnitSphereVolume[NUMDIMS];

	alias Callback = bool function(DATATYPE, void*);

	this()
	{
		static assert(MAXNODES > MINNODES);
		static assert(MINNODES > 0);

		m_root = AllocNode();
		m_root.m_level = 0;
	}

	~this()
	{
		Reset(); // Free, or reset node memory
	}

	/// Insert entry
	/// \param a_min Min of bounding rect
	/// \param a_max Max of bounding rect
	/// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
	void insert(ref const ELEMTYPE[NUMDIMS] a_min, ref const ELEMTYPE[NUMDIMS] a_max, ref const(DATATYPE) a_dataId)
	{
		debug
		{
			for(int index=0; index<NUMDIMS; ++index)
			{
				assert(a_min[index] <= a_max[index]);
			}
		}

		Branch branch;
		branch.m_data = a_dataId;
		branch.m_child = null;

		for(int axis=0; axis<NUMDIMS; ++axis)
		{
			branch.m_rect.m_min[axis] = a_min[axis];
			branch.m_rect.m_max[axis] = a_max[axis];
		}

		InsertRect(branch, &m_root, 0);
	}

	/// Remove entry
	/// \param a_min Min of bounding rect
	/// \param a_max Max of bounding rect
	/// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
	void remove(const ELEMTYPE[NUMDIMS] a_min, const ELEMTYPE[NUMDIMS] a_max, ref const(DATATYPE) a_dataId)
	{
		debug
		{
			for(int index=0; index<NUMDIMS; ++index)
			{
				assert(a_min[index] <= a_max[index]);
			}
		}

		Rect rect;

		for(int axis=0; axis<NUMDIMS; ++axis)
		{
			rect.m_min[axis] = a_min[axis];
			rect.m_max[axis] = a_max[axis];
		}

		RemoveRect(&rect, a_dataId, &m_root);
	}

	/// Find all within search rectangle
	/// \param a_min Min of search bounding rect
	/// \param a_max Max of search bounding rect
	/// \param a_searchResult Search result array.  Caller should set grow size. Function will reset, not append to array.
	/// \param a_resultCallback Callback function to return result.  Callback should return 'true' to continue searching
	/// \param a_context User context to pass as parameter to a_resultCallback
	/// \return Returns the number of entries found
	int search(const ELEMTYPE[NUMDIMS] a_min, const ELEMTYPE[NUMDIMS] a_max, Callback a_resultCallback, void* a_context)
	{
		debug
		{
			for(int index=0; index<NUMDIMS; ++index)
			{
				assert(a_min[index] <= a_max[index]);
			}
		}

		Rect rect;

		for(int axis=0; axis<NUMDIMS; ++axis)
		{
			rect.m_min[axis] = a_min[axis];
			rect.m_max[axis] = a_max[axis];
		}

		// NOTE: May want to return search result another way, perhaps returning the number of found elements here.

		int foundCount = 0;
		Search(m_root, &rect, foundCount, a_resultCallback, a_context);

		return foundCount;
	}


	/// Remove all entries from tree
	void removeAll()
	{
		// Delete all existing nodes
		Reset();

		m_root = AllocNode();
		m_root.m_level = 0;
	}

	/// Count the data elements in this container.  This is slow as no internal counter is maintained.
	int count()
	{
		int count = 0;
		CountRec(m_root, count);

		return count;
	}

protected:

	/// Minimal bounding rectangle (n-dimensional)
	struct Rect
	{
		ELEMTYPE[NUMDIMS] m_min;                      ///< Min dimensions of bounding box
		ELEMTYPE[NUMDIMS] m_max;                      ///< Max dimensions of bounding box
	}

	/// May be data or may be another subtree
	/// The parents level determines this.
	/// If the parents level is 0, then this is data
	struct Branch
	{
		Rect m_rect;                                  ///< Bounds
		Node* m_child;                                ///< Child node
		DATATYPE m_data;                              ///< Data Id
	}

	/// Node for each branch level
	struct Node
	{
		bool IsInternalNode()                         { return (m_level > 0); } // Not a leaf, but a internal node
		bool IsLeaf()                                 { return (m_level == 0); } // A leaf, contains data

		int m_count;                                  ///< Count
		int m_level;                                  ///< Leaf is zero, others positive
		Branch[MAXNODES] m_branch;                    ///< Branch
	}

	/// A link list of nodes for reinsertion after a delete operation
	struct ListNode
	{
		ListNode* m_next;                             ///< Next in list
		Node* m_node;                                 ///< Node
	}

	/// Variables for finding a split partition
	struct PartitionVars
	{
		enum { NOT_TAKEN = -1 } // indicates that position

		int[MAXNODES+1]    m_partition;
		int                m_total;
		int                m_minFill;
		int[2]             m_count;
		Rect[2]            m_cover;
		ELEMTYPEREAL[2]    m_area;

		Branch[MAXNODES+1] m_branchBuf;
		int                m_branchCount;
		Rect               m_coverSplit;
		ELEMTYPEREAL       m_coverSplitArea;
	}

	Node* AllocNode()
	{
		auto newNode = new Node;
		InitNode(newNode);
		return newNode;
	}

	void FreeNode(Node* a_node)
	{
		assert(a_node);
		destroy(a_node);
		a_node = null;
	}

	void InitNode(Node* a_node)
	{
		a_node.m_count = 0;
		a_node.m_level = -1;
	}

	void InitRect(Rect* a_rect)
	{
		for(int index = 0; index < NUMDIMS; ++index)
		{
			a_rect.m_min[index] = cast(ELEMTYPE) 0;
			a_rect.m_max[index] = cast(ELEMTYPE) 0;
		}
	}

	bool InsertRectRec(ref Branch a_branch, Node* a_node, Node** a_newNode, int a_level)
	{
		assert(a_node && a_newNode);
		assert(a_level >= 0 && a_level <= a_node.m_level);

		// recurse until we reach the correct level for the new record. data records
		// will always be called with a_level == 0 (leaf)
		if(a_node.m_level > a_level)
		{
			// Still above level for insertion, go down tree recursively
			Node* otherNode;

			// find the optimal branch for this record
			int index = PickBranch(&a_branch.m_rect, a_node);

			// recursively insert this record into the picked branch
			bool childWasSplit = InsertRectRec(a_branch, a_node.m_branch[index].m_child, &otherNode, a_level);

			if (!childWasSplit)
			{
				// Child was not split. Merge the bounding box of the new record with the
				// existing bounding box
				a_node.m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node.m_branch[index].m_rect));
				return false;
			}
			else
			{
				// Child was split. The old branches are now re-partitioned to two nodes
				// so we have to re-calculate the bounding boxes of each node
				a_node.m_branch[index].m_rect = NodeCover(a_node.m_branch[index].m_child);
				Branch branch;
				branch.m_child = otherNode;
				branch.m_rect = NodeCover(otherNode);

				// The old node is already a child of a_node. Now add the newly-created
				// node to a_node as well. a_node might be split because of that.
				return AddBranch(&branch, a_node, a_newNode);
			}
		}
		else if(a_node.m_level == a_level)
		{
			// We have reached level for insertion. Add rect, split if necessary
			return AddBranch(&a_branch, a_node, a_newNode);
		}
		else
		{
			// Should never occur
			debug
			{
				assert(0);
			}
			else
			{
				return false;
			}
		}
	}

	bool InsertRect(ref Branch a_branch, Node** a_root, int a_level)
	{
		assert(a_root);
		assert(a_level >= 0 && a_level <= (*a_root).m_level);
		debug
		{
			for(int index=0; index < NUMDIMS; ++index)
			{
				assert(a_branch.m_rect.m_min[index] <= a_branch.m_rect.m_max[index]);
			}
		}

		Node* newNode;

		if (InsertRectRec(a_branch, *a_root, &newNode, a_level))  // Root split
		{
			// Grow tree taller and new root
			Node* newRoot = AllocNode();
			newRoot.m_level = (*a_root).m_level + 1;

			Branch branch;

			// add old root node as a child of the new root
			branch.m_rect = NodeCover(*a_root);
			branch.m_child = *a_root;
			AddBranch(&branch, newRoot, null);

			// add the split node as a child of the new root
			branch.m_rect = NodeCover(newNode);
			branch.m_child = newNode;
			AddBranch(&branch, newRoot, null);

			// set the new root as the root node
			*a_root = newRoot;

			return true;
		}

		return false;
	}

	Rect NodeCover(Node* a_node)
	{
		assert(a_node);

		Rect rect = a_node.m_branch[0].m_rect;
		for(int index = 1; index < a_node.m_count; ++index)
		{
			rect = CombineRect(&rect, &(a_node.m_branch[index].m_rect));
		}

		return rect;
	}

	bool AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode)
	{
		assert(a_branch);
		assert(a_node);

		if(a_node.m_count < MAXNODES)  // Split won't be necessary
		{
			a_node.m_branch[a_node.m_count] = *a_branch;
			++a_node.m_count;

			return false;
		}
		else
		{
			assert(a_newNode);

			SplitNode(a_node, a_branch, a_newNode);
			return true;
		}
	}

	void DisconnectBranch(Node* a_node, int a_index)
	{
		assert(a_node && (a_index >= 0) && (a_index < MAXNODES));
		assert(a_node.m_count > 0);

		// Remove element by swapping with the last element to prevent gaps in array
		a_node.m_branch[a_index] = a_node.m_branch[a_node.m_count - 1];

		--a_node.m_count;
	}

	int PickBranch(const Rect* a_rect, Node* a_node)
	{
		assert(a_rect && a_node);

		bool firstTime = true;
		ELEMTYPEREAL increase;
		ELEMTYPEREAL bestIncr = cast(ELEMTYPEREAL) -1;
		ELEMTYPEREAL area;
		ELEMTYPEREAL bestArea;
		int best;
		Rect tempRect;

		for(int index=0; index < a_node.m_count; ++index)
		{
			Rect* curRect = &a_node.m_branch[index].m_rect;
			area = CalcRectVolume(curRect);
			tempRect = CombineRect(a_rect, curRect);
			increase = CalcRectVolume(&tempRect) - area;
			if((increase < bestIncr) || firstTime)
			{
				best = index;
				bestArea = area;
				bestIncr = increase;
				firstTime = false;
			}
			else if((increase == bestIncr) && (area < bestArea))
			{
				best = index;
				bestArea = area;
				bestIncr = increase;
			}
		}
		return best;
	}

	Rect CombineRect(const Rect* a_rectA, const Rect* a_rectB)
	{
		assert(a_rectA && a_rectB);

		Rect newRect;

		for(int index = 0; index < NUMDIMS; ++index)
		{
			import std.algorithm : min, max;
			newRect.m_min[index] = min(a_rectA.m_min[index], a_rectB.m_min[index]);
			newRect.m_max[index] = max(a_rectA.m_max[index], a_rectB.m_max[index]);
		}

		return newRect;
	}

	void SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode)
	{
		assert(a_node);
		assert(a_branch);

		// Could just use local here, but member or external is faster since it is reused
		PartitionVars localVars;
		PartitionVars* parVars = &localVars;

		// Load all the branches into a buffer, initialize old node
		GetBranches(a_node, a_branch, parVars);

		// Find partition
		ChoosePartition(parVars, MINNODES);

		// Create a new node to hold (about) half of the branches
		*a_newNode = AllocNode();
		(*a_newNode).m_level = a_node.m_level;

		// Put branches from buffer into 2 nodes according to the chosen partition
		a_node.m_count = 0;
		LoadNodes(a_node, *a_newNode, parVars);

		assert((a_node.m_count + (*a_newNode).m_count) == parVars.m_total);
	}

	ELEMTYPEREAL RectSphericalVolume(Rect* a_rect)
	{
		assert(a_rect);

		ELEMTYPEREAL sumOfSquares = cast(ELEMTYPEREAL) 0;
		ELEMTYPEREAL radius;

		for(int index=0; index < NUMDIMS; ++index)
		{
			ELEMTYPEREAL halfExtent = (cast(ELEMTYPEREAL) a_rect.m_max[index] - cast(ELEMTYPEREAL) a_rect.m_min[index]) * 0.5f;
			sumOfSquares += halfExtent * halfExtent;
		}

		import std.math : sqrt;
		radius = cast(ELEMTYPEREAL) sqrt(sumOfSquares);

		// Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
		static if(NUMDIMS == 3)
		{
			return (radius * radius * radius * unitSphereVolume);
		}
		else static if(NUMDIMS == 2)
		{
			return (radius * radius * unitSphereVolume);
		}
		else
		{
			return cast(ELEMTYPEREAL) (pow(radius, NUMDIMS) * unitSphereVolume);
		}
	}

	ELEMTYPEREAL RectVolume(Rect* a_rect)
	{
		assert(a_rect);

		ELEMTYPEREAL volume = cast(ELEMTYPEREAL) 1;

		for(int index=0; index<NUMDIMS; ++index)
		{
			volume *= a_rect.m_max[index] - a_rect.m_min[index];
		}

		assert(volume >= cast(ELEMTYPEREAL) 0);

		return volume;
	}

	ELEMTYPEREAL CalcRectVolume(Rect* a_rect)
	{
		version(RTREE_USE_SPHERICAL_VOLUME)
		{
			return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
		}
		else
		{
			return RectVolume(a_rect); // Faster but can cause poor merges
		}
	}

	void GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars)
	{
		assert(a_node);
		assert(a_branch);

		assert(a_node.m_count == MAXNODES);

		// Load the branch buffer
		for(int index=0; index < MAXNODES; ++index)
		{
			a_parVars.m_branchBuf[index] = a_node.m_branch[index];
		}
		a_parVars.m_branchBuf[MAXNODES] = *a_branch;
		a_parVars.m_branchCount = MAXNODES + 1;

		// Calculate rect containing all in the set
		a_parVars.m_coverSplit = a_parVars.m_branchBuf[0].m_rect;
		for(int index=1; index < MAXNODES+1; ++index)
		{
			a_parVars.m_coverSplit = CombineRect(&a_parVars.m_coverSplit, &a_parVars.m_branchBuf[index].m_rect);
		}
		a_parVars.m_coverSplitArea = CalcRectVolume(&a_parVars.m_coverSplit);
	}

	void ChoosePartition(PartitionVars* a_parVars, int a_minFill)
	{
		assert(a_parVars);

		ELEMTYPEREAL biggestDiff;
		int group, chosen, betterGroup;

		InitParVars(a_parVars, a_parVars.m_branchCount, a_minFill);
		PickSeeds(a_parVars);

		while (((a_parVars.m_count[0] + a_parVars.m_count[1]) < a_parVars.m_total)
			&& (a_parVars.m_count[0] < (a_parVars.m_total - a_parVars.m_minFill))
			&& (a_parVars.m_count[1] < (a_parVars.m_total - a_parVars.m_minFill)))
		{
			biggestDiff = cast(ELEMTYPEREAL) -1;
			for(int index=0; index<a_parVars.m_total; ++index)
			{
				if(PartitionVars.NOT_TAKEN == a_parVars.m_partition[index])
				{
					Rect* curRect = &a_parVars.m_branchBuf[index].m_rect;
					Rect rect0 = CombineRect(curRect, &a_parVars.m_cover[0]);
					Rect rect1 = CombineRect(curRect, &a_parVars.m_cover[1]);
					ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars.m_area[0];
					ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars.m_area[1];
					ELEMTYPEREAL diff = growth1 - growth0;
					if(diff >= 0)
					{
						group = 0;
					}
					else
					{
						group = 1;
						diff = -diff;
					}

					if(diff > biggestDiff)
					{
						biggestDiff = diff;
						chosen = index;
						betterGroup = group;
					}
					else if((diff == biggestDiff) && (a_parVars.m_count[group] < a_parVars.m_count[betterGroup]))
					{
						chosen = index;
						betterGroup = group;
					}
				}
			}
			Classify(chosen, betterGroup, a_parVars);
		}

		// If one group too full, put remaining rects in the other
		if((a_parVars.m_count[0] + a_parVars.m_count[1]) < a_parVars.m_total)
		{
			if(a_parVars.m_count[0] >= a_parVars.m_total - a_parVars.m_minFill)
			{
				group = 1;
			}
			else
			{
				group = 0;
			}

			for(int index=0; index<a_parVars.m_total; ++index)
			{
				if(PartitionVars.NOT_TAKEN == a_parVars.m_partition[index])
				{
					Classify(index, group, a_parVars);
				}
			}
		}

		assert((a_parVars.m_count[0] + a_parVars.m_count[1]) == a_parVars.m_total);
		assert((a_parVars.m_count[0] >= a_parVars.m_minFill) &&
			(a_parVars.m_count[1] >= a_parVars.m_minFill));
	}

	void LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
	{
		assert(a_nodeA);
		assert(a_nodeB);
		assert(a_parVars);

		for(int index=0; index < a_parVars.m_total; ++index)
		{
			assert(a_parVars.m_partition[index] == 0 || a_parVars.m_partition[index] == 1);

			int targetNodeIndex = a_parVars.m_partition[index];
			Node*[2] targetNodes = [a_nodeA, a_nodeB];

			// It is assured that AddBranch here will not cause a node split.
			bool nodeWasSplit = AddBranch(&a_parVars.m_branchBuf[index], targetNodes[targetNodeIndex], null);
			assert(!nodeWasSplit);
		}
	}

	void InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
	{
		assert(a_parVars);

		a_parVars.m_count[0] = a_parVars.m_count[1] = 0;
		a_parVars.m_area[0] = a_parVars.m_area[1] = cast(ELEMTYPEREAL) 0;
		a_parVars.m_total = a_maxRects;
		a_parVars.m_minFill = a_minFill;
		for(int index=0; index < a_maxRects; ++index)
		{
			a_parVars.m_partition[index] = PartitionVars.NOT_TAKEN;
		}
	}

	void PickSeeds(PartitionVars* a_parVars)
	{
		int seed0 = 0;
		int seed1 = seed0 + 1;
		ELEMTYPEREAL worst, waste;
		ELEMTYPEREAL[MAXNODES+1] area;

		for(int index=0; index<a_parVars.m_total; ++index)
		{
			area[index] = CalcRectVolume(&a_parVars.m_branchBuf[index].m_rect);
		}

		worst = -a_parVars.m_coverSplitArea - 1;
		for(int indexA = seed0; indexA < a_parVars.m_total-1; ++indexA)
		{
			for(int indexB = seed1; indexB < a_parVars.m_total; ++indexB)
			{
				Rect oneRect = CombineRect(&a_parVars.m_branchBuf[indexA].m_rect, &a_parVars.m_branchBuf[indexB].m_rect);
				waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
				if(waste > worst)
				{
					worst = waste;
					seed0 = indexA;
					seed1 = indexB;
				}
			}
		}

		Classify(seed0, 0, a_parVars);
		Classify(seed1, 1, a_parVars);
	}

	void Classify(int a_index, int a_group, PartitionVars* a_parVars)
	{
		assert(a_parVars);
		assert(PartitionVars.NOT_TAKEN == a_parVars.m_partition[a_index]);

		a_parVars.m_partition[a_index] = a_group;

		// Calculate combined rect
		if (a_parVars.m_count[a_group] == 0)
		{
			a_parVars.m_cover[a_group] = a_parVars.m_branchBuf[a_index].m_rect;
		}
		else
		{
			a_parVars.m_cover[a_group] = CombineRect(&a_parVars.m_branchBuf[a_index].m_rect, &a_parVars.m_cover[a_group]);
		}

		// Calculate volume of combined rect
		a_parVars.m_area[a_group] = CalcRectVolume(&a_parVars.m_cover[a_group]);

		++a_parVars.m_count[a_group];
	}

	bool RemoveRect(Rect* a_rect, ref const DATATYPE a_id, Node** a_root)
	{
		assert(a_rect && a_root);
		assert(*a_root);

		ListNode* reInsertList = null;

		if (!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))
		{
			// Found and deleted a data item
			// Reinsert any branches from eliminated nodes
			while(reInsertList)
			{
				Node* tempNode = reInsertList.m_node;

				for(int index = 0; index < tempNode.m_count; ++index)
				{
					// TODO go over this code. should I use (tempNode.m_level - 1)?
					InsertRect(tempNode.m_branch[index],
					       a_root,
					       tempNode.m_level);
				}

				ListNode* remLNode = reInsertList;
				reInsertList = reInsertList.m_next;

				FreeNode(remLNode.m_node);
				FreeListNode(remLNode);
			}

			// Check for redundant root (not leaf, 1 child) and eliminate TODO replace
			// if with while? In case there is a whole branch of redundant roots...
			if((*a_root).m_count == 1 && (*a_root).IsInternalNode())
			{
				Node* tempNode = (*a_root).m_branch[0].m_child;

				assert(tempNode);
				FreeNode(*a_root);
				*a_root = tempNode;
			}
			return false;
		}
		else
		{
			return true;
		}
	}

	bool RemoveRectRec(Rect* a_rect, ref const DATATYPE a_id, Node* a_node, ListNode** a_listNode)
	{
		assert(a_rect && a_node && a_listNode);
		assert(a_node.m_level >= 0);

		if(a_node.IsInternalNode())  // not a leaf node
		{
			for(int index = 0; index < a_node.m_count; ++index)
			{
				if(Overlap(a_rect, &(a_node.m_branch[index].m_rect)))
				{
					if(!RemoveRectRec(a_rect, a_id, a_node.m_branch[index].m_child, a_listNode))
					{
						if(a_node.m_branch[index].m_child.m_count >= MINNODES)
						{
							// child removed, just resize parent rect
							a_node.m_branch[index].m_rect = NodeCover(a_node.m_branch[index].m_child);
						}
						else
						{
							// child removed, not enough entries in node, eliminate node
							ReInsert(a_node.m_branch[index].m_child, a_listNode);
							DisconnectBranch(a_node, index); // Must return after this call as count has changed
						}
						return false;
					}
				}
			}
			return true;
		}
		else // A leaf node
		{
			for(int index = 0; index < a_node.m_count; ++index)
			{
				if(a_node.m_branch[index].m_data == a_id)
				{
					DisconnectBranch(a_node, index); // Must return after this call as count has changed
					return false;
				}
			}
			return true;
		}
	}

	ListNode* AllocListNode()
	{
		return new ListNode;
	}

	void FreeListNode(ListNode* a_listNode)
	{
		destroy(a_listNode);
	}

	bool Overlap(Rect* a_rectA, Rect* a_rectB)
	{
		assert(a_rectA && a_rectB);

		for(int index=0; index < NUMDIMS; ++index)
		{
			if (a_rectA.m_min[index] > a_rectB.m_max[index] ||
				a_rectB.m_min[index] > a_rectA.m_max[index])
			{
				return false;
			}
		}
		return true;
	}

	void ReInsert(Node* a_node, ListNode** a_listNode)
	{
		ListNode* newListNode;

		newListNode = AllocListNode();
		newListNode.m_node = a_node;
		newListNode.m_next = *a_listNode;
		*a_listNode = newListNode;
	}

	bool Search(Node* a_node, Rect* a_rect, ref int a_foundCount, Callback a_resultCallback, void* a_context)
	{
		assert(a_node);
		assert(a_node.m_level >= 0);
		assert(a_rect);

		if(a_node.IsInternalNode())
		{
			// This is an internal node in the tree
			for(int index=0; index < a_node.m_count; ++index)
			{
				if(Overlap(a_rect, &a_node.m_branch[index].m_rect))
				{
					if(!Search(a_node.m_branch[index].m_child, a_rect, a_foundCount, a_resultCallback, a_context))
					{
						// The callback indicated to stop searching
						return false;
					}
				}
			}
		}
		else
		{
			// This is a leaf node
			for(int index=0; index < a_node.m_count; ++index)
			{
				if(Overlap(a_rect, &a_node.m_branch[index].m_rect))
				{
					DATATYPE* id = &a_node.m_branch[index].m_data;
					++a_foundCount;

					// NOTE: There are different ways to return results.  Here's where to modify
					if(a_resultCallback)
					{
						if(!a_resultCallback(*id, a_context))
						{
							return false; // Don't continue searching
						}
					}
				}
			}
		}

		return true; // Continue searching
	}

	void RemoveAllRec(Node* a_node)
	{
		assert(a_node);
		assert(a_node.m_level >= 0);

		if(a_node.IsInternalNode()) // This is an internal node in the tree
		{
			for(int index=0; index < a_node.m_count; ++index)
			{
				RemoveAllRec(a_node.m_branch[index].m_child);
			}
		}
		FreeNode(a_node);
	}

	void Reset()
	{
		RemoveAllRec(m_root);
	}

	void CountRec(Node* a_node, ref int a_count)
	{
		if(a_node.IsInternalNode())  // not a leaf node
		{
			for(int index = 0; index < a_node.m_count; ++index)
			{
				CountRec(a_node.m_branch[index].m_child, a_count);
			}
		}
		else // A leaf node
		{
			a_count += a_node.m_count;
		}
	}

	Node* m_root;                                    ///< Root of tree
}
