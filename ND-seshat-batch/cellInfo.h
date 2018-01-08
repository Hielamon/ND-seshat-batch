#pragma once

//Structure to handle coordinates
struct coo {
	int x, y, s, t;

	coo() {
		x = y = s = t = -1;
	}

	coo(int a, int b, int c, int d) {
		x = a; y = b; s = c; t = d;
	}

	bool operator==(coo &R) {
		return x == R.x && y == R.y && s == R.s && t == R.t;
	}
};

inline bool operator<(const coo &A, const coo &B)
{
	if (A.x < B.x) return true;
	if (A.x == B.x) {
		if (A.y < B.y) return true;
		if (A.y == B.y) {
			if (A.s < B.s) return true;
			if (A.s == B.s)
				if (A.t < B.t) return true;
		}
	}
	return false;
}

//The class for warping the useful infomation in a cell
class CellInfo
{
public:
	CellInfo(coo &box_, std::vector<bool> &segMask_, int talla_)
		: box(box_), segN(segMask_.size()), segMask(segMask_), talla(talla_)
	{}
	
	CellInfo() {}

	~CellInfo() {}

	//The bounding box of the cell
	coo box;

	//total number of segments
	int segN;

	//the mask for determining whether the segment is convered by this cell
	std::vector<bool> segMask;

	//total number of segments covered in this cell
	int talla; 
};

