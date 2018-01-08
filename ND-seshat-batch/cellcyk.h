/*Copyright 2014 Francisco Alvaro

 This file is part of SESHAT.

    SESHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SESHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SESHAT.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _CELLCYK_
#define _CELLCYK_
#include <cstdio>
#include <vector>
#include "hypothesis.h"
//#include "sample.h"
#include "cellInfo.h"


class CellCYK{
public:
	//Cell info , include bounding-box, segMask, and talla
	std::shared_ptr<CellInfo> pCInfo;

	//total number of all non-terminals
	int nnt;
	//Hypotheses for every non-terminals
	std::vector<std::shared_ptr<Hypothesis>> vNoTerm;

	//Next cell in linked list (CYK table of same size)
	std::shared_ptr<CellCYK> sig;

	//Methods
	CellCYK(int n, int segN)
	{
		sig = NULL;
		nnt = n;

		pCInfo = std::make_shared<CellInfo>();
		pCInfo->segMask.resize(segN);
		pCInfo->segN = segN;
		pCInfo->talla = 0;

		//Create (empty) hypotheses
		vNoTerm.resize(nnt);
	}
	~CellCYK() {}

	void setRegion(coo &box, const int &segIdx)
	{
		pCInfo->segMask[segIdx] = true;
		pCInfo->box = box;
	}

	void setRegion(Sample &M, const int &segIdx)
	{
		pCInfo->segMask[segIdx] = true;
		SegUnit &seg = M.getSegUnit(segIdx);
		pCInfo->box.x = seg.ROI.x;
		pCInfo->box.y = seg.ROI.y;
		pCInfo->box.s = seg.ROI.br().x;
		pCInfo->box.t = seg.ROI.br().y;
	}

	//Comparison operator for logspace ordering
	bool operator<(const CellCYK &C)
	{
		coo &boxC = C.pCInfo->box;
		coo &box = pCInfo->box;
		if (box.x < boxC.x)
			return true;
		if (box.x == boxC.x) {
			if (box.y < boxC.y)
				return true;
			if (box.y == boxC.y) {
				if (box.s < boxC.s)
					return true;
				if (box.s == boxC.s)
					if (box.t < boxC.t)
						return true;
			}
		}
		return false;
	}

	//Set the covered strokes to the union of cells A and B
	void ccUnion(CellCYK &A, CellCYK &B)
	{
		for (int i = 0; i < pCInfo->segN; i++)
			pCInfo->segMask[i] = (A.pCInfo->segMask[i] || B.pCInfo->segMask[i]) ? true : false;
	}

	//Check if cell H covers the same strokes that this
	bool ccEqual(std::shared_ptr<CellCYK> &pCell)
	{
		if (pCInfo->talla != pCell->pCInfo->talla)
			return false;

		for (int i = 0; i<pCInfo->segN; i++)
			if (pCInfo->segMask[i] != pCell->pCInfo->segMask[i])
				return false;

		return true;
	}

	//Check if the intersection between the strokes of this cell and H is empty
	bool compatible(std::shared_ptr<CellCYK> &pCell)
	{
		for (int i = 0; i<pCInfo->segN; i++)
			if (pCInfo->segMask[i] && pCell->pCInfo->segMask[i])
				return false;

		return true;
	}

  
};

inline void MergeRegionsCenter(std::shared_ptr<CellCYK>& pCA, int ntIDA,
						 std::shared_ptr<CellCYK>& pCB, int ntIDB,
						 std::shared_ptr<CellCYK>& pCS, int ntIDS,
						 char merge_cen)
{
	std::shared_ptr<Hypothesis> &a = pCA->vNoTerm[ntIDA];
	std::shared_ptr<Hypothesis> &b = pCB->vNoTerm[ntIDB];
	std::shared_ptr<Hypothesis> &s = pCS->vNoTerm[ntIDS];
	switch (merge_cen) {
	case 'A': //Data Hypothesis a
		s->lcen = a->lcen;
		s->rcen = a->rcen;
		break;
	case 'B': //Data Hypothesis b
		s->lcen = b->lcen;
		s->rcen = b->rcen;
		break;
	case 'C': //Center point
		s->lcen = (pCA->pCInfo->box.y + pCA->pCInfo->box.t) / 2;
		s->rcen = (pCB->pCInfo->box.y + pCB->pCInfo->box.t) / 2;
		break;
	case 'M': //Mean of both centers
		s->lcen = (a->lcen + b->lcen) / 2; //a->lcen;
		s->rcen = (a->rcen + b->rcen) / 2; //b->rcen;
		break;
	default:
		HL_CERR("Error: Unrecognized option " << merge_cen << " in merge regions");
	}
}


#endif
