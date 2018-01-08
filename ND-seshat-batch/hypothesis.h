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
#ifndef _HYPOTHESIS_
#define _HYPOTHESIS_

//class ProductionB;
//class ProductionT;
//
//struct CellCYK;
//struct Grammar;

#include <cstdio>
#include <cmath>
#include <list>
#include <memory>
#include "production.h"
#include "cellInfo.h"
//#include "cellcyk.h"
//#include "grammar.h"



struct Hypothesis{
	int clase; //If the hypothesis encodes a terminal symbols this is the class id (-1 otherwise)
	double pr; //log-probability

	//Cell info , include bounding-box, segMask, and talla
	std::shared_ptr<CellInfo> pCInfo;
	
	//References to left-child (hi) and right-child (hd) to create the derivation tree
	std::shared_ptr<Hypothesis> hleft, hright;

	//The production used to create this hypothesis (either Binary or terminal)
	std::shared_ptr<ProductionB> prod;
	std::shared_ptr<ProductionT> pt;
	//int ptID;

	//Auxiliar var to retrieve the used production in the special SSE treatment
	std::shared_ptr<ProductionB> prod_sse;

	//Vertical center left and right
	int lcen, rcen;

	//CellCYK *parent; //Parent cell
	//int ntid;        //Nonterminal ID in parent

					 //Methods
	Hypothesis(int c, double p, std::shared_ptr<CellInfo> pCInfo_)
	{
		pCInfo = pCInfo_;
		clase = c;
		pr = p;
		//hright = hleft = NULL;
		lcen = rcen = 0;
	}
	~Hypothesis() {}

	void copy(std::shared_ptr<Hypothesis> &H)
	{
		clase = H->clase;
		pr = H->pr;
		hleft = H->hleft;
		hright = H->hright;
		lcen = H->lcen;
		rcen = H->rcen;
		pCInfo = H->pCInfo;
		prod = H->prod;
		pt = H->pt;
		prod_sse = H->prod_sse;
	}

	//Check if the intersection between the strokes of this cell and H is empty
	bool compatible(std::shared_ptr<Hypothesis> &pH)
	{
		for (int i = 0; i<pCInfo->segN; i++)
			if (pCInfo->segMask[i] && pH->pCInfo->segMask[i])
				return false;

		return true;
	}
};

#endif
