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
#ifndef _TABLECYK_
#define _TABLECYK_

#include <cstdio>
#include <map>
#include <list>
#include <climits>
#include <cfloat>
#include <memory>
#include "cellcyk.h"
#include "hypothesis.h"




class TableCYK{
	std::vector<std::shared_ptr<CellCYK>> T;
	std::vector<std::map<coo, std::shared_ptr<CellCYK>>> TS;
	int N, K;

	//Hypothesis that accounts for the target (input) math expression
	std::shared_ptr<Hypothesis> Target;
	std::shared_ptr<CellCYK> targetC;
	int targetHIndex;

	//Percentage of strokes covered by the most likely hypothesis (target)
	int pm_comps;
 
public:
	TableCYK(int n, int k)
	{
		N = n;
		K = k;

		//Target = NULL;
		Target =std::make_shared<Hypothesis>(-1, -FLT_MAX, std::shared_ptr<CellInfo>());
		targetHIndex = -1;
		pm_comps = 0.0;

		T.resize(N);
		TS.resize(N);
	}
	~TableCYK()
	{
	}

	void getMLInfo(std::shared_ptr<Hypothesis> &H, std::shared_ptr<CellCYK> &C, int &HIndex)
	{
		H = Target;
		HIndex = targetHIndex;
		C = targetC;
	}

	std::shared_ptr<CellCYK>& get(int n)
	{
		return T[n - 1];
	}

	int size(int n)
	{
		return TS[n - 1].size();
	}

	void updateTarget(std::shared_ptr<CellCYK>& pCell, int HIndex)
	{

		if (HIndex < 0 || HIndex >= pCell->nnt)
			HL_CERR("the H index is invalid!");

		int pcomps = 0;

		for (int i = 0; i < pCell->pCInfo->segN; i++)
			if (pCell->pCInfo->segMask[i])
				pcomps++;

		if (pcomps > pm_comps || (pcomps == pm_comps && pCell->vNoTerm[HIndex]->pr > Target->pr)) {
			pm_comps = pcomps;
			Target = pCell->vNoTerm[HIndex];
			targetHIndex = HIndex;
			targetC = pCell;
		}
	}

	void add(int n, std::shared_ptr<CellCYK>& pCell, int noterm_id, std::vector<bool> &esinit)
	{
		coo key(pCell->pCInfo->box.x, pCell->pCInfo->box.y, pCell->pCInfo->box.s, pCell->pCInfo->box.t);
		std::map<coo, std::shared_ptr<CellCYK>>::iterator it = TS[n - 1].find(key);

		pCell->pCInfo->talla = n;

		if (it == TS[n - 1].end()) {
			//Link as head of size  n
			pCell->sig = T[n - 1];
			T[n - 1] = pCell;
			TS[n - 1][key] = pCell;

			if (noterm_id >= 0) {
				if (esinit[noterm_id])
					updateTarget(pCell, noterm_id);
			}
			else {
				for (int nt = 0; nt<pCell->nnt; nt++)
					if (pCell->vNoTerm[nt].use_count() && esinit[nt])
						updateTarget(pCell, nt);
			}
		}
		else { //Maximize probability avoiding duplicates

			int VA, VB;
			if (noterm_id < 0) {
				VA = 0;
				VB = pCell->nnt;
			}
			else {
				VA = noterm_id;
				VB = VA + 1;
			}

			std::shared_ptr<CellCYK> r = it->second;

			if (!pCell->ccEqual(r)) {
				//The cells cover the same region with a different set of strokes

				float maxpr_c = -FLT_MAX;
				for (int i = VA; i<VB; i++)
					if (pCell->vNoTerm[i].use_count() && pCell->vNoTerm[i]->pr > maxpr_c)
						maxpr_c = pCell->vNoTerm[i]->pr;

				float maxpr_r = -FLT_MAX;
				for (int i = 0; i<r->nnt; i++)
					if (r->vNoTerm[i].use_count() && r->vNoTerm[i]->pr > maxpr_r)
						maxpr_r = r->vNoTerm[i]->pr;

				//If the new cell contains the most likely hypothesis, replace the hypotheses
				if (maxpr_c > maxpr_r) {

					//Copy the new set of strokes
					for (int i = 0; i<pCell->pCInfo->segN; i++)
						r->pCInfo->segMask[i] = pCell->pCInfo->segMask[i];

					//Replace the hypotheses for each non-terminal
					for (int i = 0; i<pCell->nnt; i++)
						if (pCell->vNoTerm[i].use_count()) {

							if (r->vNoTerm[i].use_count()) {
								r->vNoTerm[i]->copy(pCell->vNoTerm[i]);

								if (esinit[i])
									updateTarget(r, i);
							}
							else {
								r->vNoTerm[i] = pCell->vNoTerm[i];

								//Set to NULL such that the "delete celda" doesn't delete the hypothesis
								pCell->vNoTerm[i].reset();

								if (esinit[i])
									updateTarget(r, i);
							}

						}
						else if (r->vNoTerm[i].use_count()) {
							r->vNoTerm[i].reset();
						}

				}

				//Finished
			}
			else
			{
				for (int i = VA; i<VB; i++) {

					if (pCell->vNoTerm[i].use_count()) {
						if (r->vNoTerm[i].use_count()) {

							if (pCell->vNoTerm[i]->pr > r->vNoTerm[i]->pr) {
								//Maximize probability (replace)
								r->vNoTerm[i]->copy(pCell->vNoTerm[i]);

								if (esinit[i])
									updateTarget(r, i);
							}
						}
						else {
							r->vNoTerm[i] = pCell->vNoTerm[i];

							//Set to NULL such that the "delete celda" doesn't delete the hypothesis
							pCell->vNoTerm[i].reset();

							if (esinit[i])
								updateTarget(r, i);
						}
					}
				}
			}

			pCell.reset();
		}
	}

	
};


#endif
