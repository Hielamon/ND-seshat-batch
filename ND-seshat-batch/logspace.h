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
#ifndef _LOGSPACE_
#define _LOGSPACE_

#include <cstdio>
#include <list>
#include <algorithm>
#include "cellcyk.h"

class LogSpace {
	int N;
	int RX, RY;
	std::vector<std::shared_ptr<CellCYK>> data;

	void quicksort(std::vector<std::shared_ptr<CellCYK>> &vec, int ini, int fin)
	{
		if (ini < fin) {
			int piv = partition(vec, ini, fin);
			quicksort(vec, ini, piv);
			quicksort(vec, piv + 1, fin);
		}
	}
	int partition(std::vector<std::shared_ptr<CellCYK>> &vec, int ini, int fin)
	{
		int piv = vec[ini]->pCInfo->box.x;
		int i = ini - 1, j = fin + 1;

		do {
			do {
				j--;
			} while (vec[j]->pCInfo->box.x > piv);
			do {
				i++;
			} while (vec[i]->pCInfo->box.x < piv);

			if (i<j) {
				std::shared_ptr<CellCYK> aux = vec[i];
				vec[i] = vec[j];
				vec[j] = aux;
			}
		} while (i<j);

		return j;
	}
	void bsearch(int sx, int sy, int ss, int st, std::list<std::shared_ptr<CellCYK>> &set)
	{
		//Binary search of "sx"
		int i, j;
		for (i = 0, j = N; i<j; ) {
			int m = (i + j) / 2;

			if (sx <= data[m]->pCInfo->box.x)
				j = m;
			else
				i = m + 1;
		}

		//Retrieve the compatible regions
		while (i<N && data[i]->pCInfo->box.x <= ss) {
			if (data[i]->pCInfo->box.y <= st && data[i]->pCInfo->box.t >= sy) {
				set.push_back(data[i]);
			}
			i++;
		}
	}
	void bsearchStv(int sx, int sy, int ss, int st, std::list<std::shared_ptr<CellCYK>> &set, bool U_V, std::shared_ptr<CellCYK> &cd)
	{
		//Binary search of "sx"
		int i, j;
		for (i = 0, j = N; i<j; ) {
			int m = (i + j) / 2;

			if (sx <= data[m]->pCInfo->box.x)
				j = m;
			else
				i = m + 1;
		}

		//Retrieve the compatible regions
		if (U_V) { //Direction 'Up' (U)
			while (i<N && data[i]->pCInfo->box.x <= ss) {
				if (data[i]->pCInfo->box.t <= st && data[i]->pCInfo->box.t >= sy && data[i]->pCInfo->box.s <= ss) {

					//For checking in the x-coordinate
					bool checkRange = data[i]->pCInfo->box.s < cd->pCInfo->box.s &&
						data[i]->pCInfo->box.x > cd->pCInfo->box.s;

					if (data[i]->pCInfo->box.t < cd->pCInfo->box.y && checkRange)
						sy = std::max(std::max(data[i]->pCInfo->box.y, data[i]->pCInfo->box.t - RY), sy);
					set.push_back(data[i]);
				}
				i++;
			}
		}
		else { //Direction 'Down' (V)
			while (i<N && data[i]->pCInfo->box.x <= ss) {
				if (data[i]->pCInfo->box.y <= st && data[i]->pCInfo->box.y >= sy && data[i]->pCInfo->box.s <= ss) {

					//For checking in the x-coordinate
					bool checkRange = data[i]->pCInfo->box.s < cd->pCInfo->box.s &&
						data[i]->pCInfo->box.x > cd->pCInfo->box.s;

					if (data[i]->pCInfo->box.y > cd->pCInfo->box.t && checkRange)
						st = std::min(std::min(data[i]->pCInfo->box.t, data[i]->pCInfo->box.y + RY), st);
					set.push_back(data[i]);
				}
				i++;
			}
		}
	}
	void bsearchHBP(int sx, int sy, int ss, int st, std::list<std::shared_ptr<CellCYK>> &set, std::shared_ptr<CellCYK> &cd)
	{
		//Binary search of "sx"
		int i, j;
		for (i = 0, j = N; i<j; ) {
			int m = (i + j) / 2;

			if (sx <= data[m]->pCInfo->box.x)
				j = m;
			else
				i = m + 1;
		}

		//Retrieve the compatible regions
		while (i<N && data[i]->pCInfo->box.x <= ss) {
			if (data[i]->pCInfo->box.y <= st && data[i]->pCInfo->box.t >= sy) {

				//for checking the overlap in y-coordinate
				int tly = std::max(data[i]->pCInfo->box.y, cd->pCInfo->box.y);
				int bry = std::min(data[i]->pCInfo->box.t, cd->pCInfo->box.t);

				if (data[i]->pCInfo->box.x > cd->pCInfo->box.s && bry > tly)
				{
					//ss = std::min(data[i]->pCInfo->box.s, ss);
					ss = std::min(std::min(data[i]->pCInfo->box.s, data[i]->pCInfo->box.x + RX), ss);
				}
				set.push_back(data[i]);
			}
			i++;
		}
	}

public:
	LogSpace(std::shared_ptr<CellCYK> &pCell, int nr, int dx, int dy)
	{
		//List length
		N = nr;
		//Size of the "reference symbol"
		RX = dx;
		RY = dy;

		//Create a vector to store the regions
		data.resize(N);
		int i = 0;
		for (std::shared_ptr<CellCYK> r = pCell; r.use_count(); r = r->sig)
			data[i++] = r;

		//Sort the regions
		quicksort(data, 0, N - 1);
	}
	~LogSpace() {}

	void getH(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = std::max(pCell->pCInfo->box.x + 1, pCell->pCInfo->box.s - (int)(RX * 2));  // (sx,sy)------
		ss = pCell->pCInfo->box.s + RX * 10;                    //  ------------
		sy = pCell->pCInfo->box.y - 2 * RY;                      //  ------------
		st = pCell->pCInfo->box.t + 2 * RY;                      //  ------(ss,st)

																 //Retrieve the regions
		bsearchHBP(sx, sy, ss, st, set, pCell);
	}

	//Below region
	void getV(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = pCell->pCInfo->box.x - 2 * RX;
		ss = pCell->pCInfo->box.s + 2 * RX;
		sy = std::max(pCell->pCInfo->box.t - RY, pCell->pCInfo->box.y + 1);
		st = pCell->pCInfo->box.t + RY * 3;

		//Retrieve the regions
		bsearchStv(sx, sy, ss, st, set, false, pCell);
	}

	//Above region. Although only Below is really considered this is necessary to
	//solve the problem of the case | aaa|
	//                              |bbbb|
	//such that "a" would never find "b" because its 'sx' would start before "b.x"
	void getU(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = pCell->pCInfo->box.x - 2 * RX;
		ss = pCell->pCInfo->box.s + 2 * RX;
		sy = pCell->pCInfo->box.y - RY * 3;
		st = std::min(pCell->pCInfo->box.y + RY, pCell->pCInfo->box.t - 1);

		//Retrieve the regions
		bsearchStv(sx, sy, ss, st, set, true, pCell);
	}

	//Inside region (sqrt)
	void getI(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = pCell->pCInfo->box.x + 1;  // (sx,sy)------
		ss = pCell->pCInfo->box.s + RX; //  ------------
										//ss = pCell->pCInfo->box.s + 1; //  ------------
		sy = pCell->pCInfo->box.y + 1;  //  ------------
										//sy = pCell->pCInfo->box.y - 2*RY;  //  ------------
		st = pCell->pCInfo->box.t + RY; //  ------(ss,st)

										//Retrieve the regions
		bsearch(sx, sy, ss, st, set);
	}

	//Mroot region (n-th sqrt)
	void getM(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = pCell->pCInfo->box.x - 2 * RX;            // (sx,sy)------
		ss = std::min(pCell->pCInfo->box.x + 2 * RX, pCell->pCInfo->box.s); //  ------------
		sy = pCell->pCInfo->box.y - RY;              //  ------------
		st = std::min(pCell->pCInfo->box.y + 2 * RY, pCell->pCInfo->box.t); //  ------(ss,st)

																			//Retrieve the regions
		bsearch(sx, sy, ss, st, set);
	}

	//SubSupScript regions
	void getS(std::shared_ptr<CellCYK> &pCell, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = pCell->pCInfo->box.x - 1;      // (sx,sy)------
		ss = pCell->pCInfo->box.x + 1;      //  ------------
		sy = pCell->pCInfo->box.y - RY;   //  ------------
		st = pCell->pCInfo->box.t + RY;   //  ------(ss,st)

		bsearch(sx, sy, ss, st, set);
	}

	void getS(coo &box, std::list<std::shared_ptr<CellCYK>> &set)
	{
		int sx, sy, ss, st;

		//Set the region to search
		sx = box.x - 1;      // (sx,sy)------
		ss = box.x + 1;      //  ------------
		sy = box.y - RY;   //  ------------
		st = box.t + RY;   //  ------(ss,st)

		bsearch(sx, sy, ss, st, set);
	}
};

#endif
