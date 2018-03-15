#pragma once
#include <algorithm>
#include <sstream>
#include "symSet.h"
#include "grammar.h"
#include "sample.h"
#include "tablecyk.h"
#include "logspace.h"
#include "sparel.h"
#include "relationSet.h"

#define NB 1

//#define SHOW_PIPELINE
//#define LOG_COMPUTE

inline void PrintSymSeg(std::shared_ptr<Hypothesis> &H)
{
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count()) {
		PrintSymSeg(H->hleft);
		PrintSymSeg(H->hright);
	}
	else {
		std::string clatex = H->pt->getTeX(H->clase);
		std::cout << clatex << std::endl;
	}
}

inline void PrintLatex(std::shared_ptr<Hypothesis> &H, std::shared_ptr<Grammar> &pG) {
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count())
	{
		//H->prod->printOut(pG, H);
		std::string &outStr = H->prod->get_outstr();
		if (!outStr.empty()) {

			int pd1 = outStr.find("$1");
			int pd2 = outStr.find("$2");

			int i = 0;
			if (pd2 >= 0 && pd1 >= 0 && pd2 < pd1) 
			{
				while (outStr[i] != '$' || outStr[i + 1] != '2') 
				{
					putchar(outStr[i]);
					i++;
				}
				i += 2;

				PrintLatex(H->hright, pG);
				if (H->hright->clase < 0)

				while (outStr[i] != '$' || outStr[i + 1] != '1')
				{
					putchar(outStr[i]);
					i++;
				}
				i += 2;

				PrintLatex(H->hleft, pG);
			}
			else
			{
				if (pd1 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '1')
					{
						putchar(outStr[i]);
						i++;
					}
					i += 2;

					PrintLatex(H->hleft, pG);
				}
				if (pd2 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '2')
					{
						putchar(outStr[i]);
						i++;
					}
					i += 2;

					PrintLatex(H->hright, pG);
				}
			}

			while (outStr[i]) 
			{
				putchar(outStr[i]);
				i++;
			}
		}
		
	}
	else {
		std::string clatex = H->pt->getTeX(H->clase);
		std::cout << clatex;
	}
	//std::cout << std::endl;
}

inline void PrintLatexToString(std::shared_ptr<Hypothesis> &H, std::shared_ptr<Grammar> &pG, std::stringstream &result) {
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count())
	{
		//H->prod->printOut(pG, H);
		std::string &outStr = H->prod->get_outstr();
		if (!outStr.empty()) {

			int pd1 = outStr.find("$1");
			int pd2 = outStr.find("$2");

			int i = 0;
			if (pd2 >= 0 && pd1 >= 0 && pd2 < pd1)
			{
				while (outStr[i] != '$' || outStr[i + 1] != '2')
				{
					result << outStr[i];
					i++;
				}
				i += 2;

				PrintLatexToString(H->hright, pG, result);
				if (H->hright->clase < 0)
					while (outStr[i] != '$' || outStr[i + 1] != '1')
					{
						result << outStr[i];
						i++;
					}
				i += 2;

				PrintLatexToString(H->hleft, pG, result);
			}
			else
			{
				if (pd1 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '1')
					{
						result << outStr[i];
						i++;
					}
					i += 2;

					PrintLatexToString(H->hleft, pG, result);
				}
				if (pd2 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '2')
					{
						result << outStr[i];
						i++;
					}
					i += 2;

					PrintLatexToString(H->hright, pG, result);
				}
			}

			while (outStr[i])
			{
				result << outStr[i];
				i++;
			}
		}

	}
	else {
		std::string clatex = H->pt->getTeX(H->clase);
		result << clatex;
	}
	//std::cout << std::endl;
}

inline void drawCell(cv::Mat &src, std::shared_ptr<CellCYK> &pCell)
{
	coo &box = pCell->pCInfo->box;
	cv::Rect roi(box.x, box.y, box.s - box.x, box.t - box.y);
	cv::rectangle(src, roi, RandomColor(), 1);
}

inline void drawCellWithColor(cv::Mat &src, std::shared_ptr<CellCYK> &pCell, cv::Scalar &color)
{
	coo &box = pCell->pCInfo->box;
	cv::Rect roi(box.x, box.y, box.s - box.x, box.t - box.y);
	cv::rectangle(src, roi, color, 3);
}

inline void PrintRelationSet(std::shared_ptr<Hypothesis> &H, std::shared_ptr<RelationSet> &pRelSet,
							 std::shared_ptr<Sample> &M)
{
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count())
	{
		char type = H->prod->tipo();
		RelationType relType;
		std::shared_ptr<Hypothesis> &Hleft = H->hleft, Hright = H->hright;
		std::shared_ptr<Hypothesis> h1, h2;
		switch (type)
		{
		case 'H':
			relType = RelationType::H;
			h1 = Hleft->hRight.use_count() == 0 ? Hleft : Hleft->hRight;
			h2 = Hright->hLeft.use_count() == 0 ? Hright : Hright->hLeft;
			break;
		case 'B':
			relType = RelationType::SUB;
			h1 = Hleft->hRight.use_count() == 0 ? Hleft : Hleft->hRight;
			h2 = Hright->hLeft.use_count() == 0 ? Hright : Hright->hLeft;
			break;
		case 'P':
			relType = RelationType::SUP;
			h1 = Hleft->hRight.use_count() == 0 ? Hleft : Hleft->hRight;
			h2 = Hright->hLeft.use_count() == 0 ? Hright : Hright->hLeft;
			break;
		case 'V':
			relType = RelationType::V;
			h1 = Hleft->hBottom.use_count() == 0 ? Hleft : Hleft->hBottom;
			h2 = Hright->hTop.use_count() == 0 ? Hright : Hright->hTop;
			break;
		case 'e':
			relType = RelationType::V;
			h1 = Hleft->hBottom.use_count() == 0 ? Hleft : Hleft->hBottom;
			h2 = Hright->hTop.use_count() == 0 ? Hright : Hright->hTop;
			break;
		case 'I':
			relType = RelationType::INS;
			h1 = Hleft->hRight.use_count() == 0 ? Hleft : Hleft->hRight;
			h2 = Hright->hRight.use_count() == 0 ? Hright : Hright->hRight;
			break;
		case 'M':
			relType = RelationType::MROOT;
			h1 = Hleft->hRight.use_count() == 0 ? Hleft : Hleft->hRight;
			h2 = Hright->hLeft.use_count() == 0 ? Hright : Hright->hLeft;
			break;
		case 'S':
			relType = RelationType::SUP;
			h1 = Hleft->hleft->hRight.use_count() == 0 ? Hleft->hleft : Hleft->hleft->hRight;
			h2 = Hright->hLeft.use_count() == 0 ? Hright : Hright->hLeft;
			break;
		default:
			break;
		}

		if (h1->pt.use_count() == 0 || h2->pt.use_count() == 0)
			HL_CERR("The impossible case happened, the relation input must be single symbol");

		coo &box1 = h1->pCInfo->box, &box2 = h2->pCInfo->box;
		RelationUnit relUnit;
		relUnit.relType = relType;
		relUnit.ROI1 = cv::Rect(box1.x, box1.y, box1.s - box1.x, box1.t - box1.y);
		relUnit.ROI2 = cv::Rect(box2.x, box2.y, box2.s - box2.x, box2.t - box2.y);
		M->getSegPoints(relUnit.ROI1, relUnit.vPoint1);
		M->getSegPoints(relUnit.ROI2, relUnit.vPoint2);

		pRelSet->addRelation(relUnit);
		PrintRelationSet(Hleft, pRelSet, M);
		PrintRelationSet(Hright, pRelSet, M);
	}
	else {
		return;
	}
	//std::cout << std::endl;
}

class MeParser
{
public:
	MeParser(const std::shared_ptr<SymSet> &pSymSet_,
			 const std::shared_ptr<Grammar> &pG_,
			 const std::shared_ptr<GMM> &pGMM_)
		: pSymSet(pSymSet_), pG(pG_), pGMM(pGMM_) 
	{
		//ptfactor = 0.37846575;
		ptfactor = 0.0;
		//pbfactor = 0.14864657;
		pbfactor = 1.0;

		//clusterF = 0.15680604;
		clusterF = 1;
		//rfactor = 0.63225349;
		dfactor = 1.0;
		rfactor = 1.0;

		qfactor = 1.88593577;
		//InsPen = 2.11917745;
		InsPen = 1;
	}
	MeParser() {}
	~MeParser()	{}

	void reSetup(const std::shared_ptr<SymSet> &pSymSet_,
				 const std::shared_ptr<Grammar> &pG_,
				 const std::shared_ptr<GMM> &pGMM_)
	{
		pSymSet = pSymSet_;
		pG = pG_;
		pGMM = pGMM_;
	}

	std::string parse(std::shared_ptr<Sample> &M)
	{
		parseCore(M);

		std::shared_ptr<Hypothesis> mlh;
		std::shared_ptr<CellCYK> mlc;
		int mlhIdx;
		tcyk->getMLInfo(mlh, mlc, mlhIdx);

		assert(mlh == mlc->vNoTerm[mlhIdx]);

		if (!mlh.use_count() || !mlc.use_count())
			HL_CERR("\nNo hypothesis found!!");

		std::cout << "\nMost Likely Hypothesis " << mlc->pCInfo->talla << " Segmentations" << std::endl;

		std::cout << "Math Symbols : " << std::endl;
		PrintSymSeg(mlh);
		std::cout << std::endl;

		std::cout << "Latex : " << std::endl;
		//PrintLatex(mlh, pG);
		//std::cout << std::endl;
		std::string latexStr;
		std::stringstream ioStr;
		PrintLatexToString(mlh, pG, ioStr);
		latexStr = ioStr.str();
		std::cout << latexStr << std::endl;
		return latexStr;
	}

	std::shared_ptr<RelationSet> getRelationSet(std::shared_ptr<Sample> &M)
	{
		if (tcyk.use_count() == 0) return nullptr;

		std::shared_ptr<RelationSet> relSet = std::make_shared<RelationSet>();

		std::shared_ptr<Hypothesis> mlh;
		std::shared_ptr<CellCYK> mlc;
		int mlhIdx;
		tcyk->getMLInfo(mlh, mlc, mlhIdx);

		assert(mlh == mlc->vNoTerm[mlhIdx]);

		if (!mlh.use_count() || !mlc.use_count())
			HL_CERR("\nNo hypothesis found!!");

		PrintRelationSet(mlh, relSet, M);
		return relSet;
	}

private:
	std::shared_ptr<SymSet> pSymSet;
	std::shared_ptr<Grammar> pG;
	std::shared_ptr<GMM> pGMM;
	std::shared_ptr<SpaRel> pSPR;

	std::shared_ptr<TableCYK> tcyk;

	//Penalty parameters for distance
	//ClusterF 
	float clusterF, dfactor;

	//ProductionTSF, ProductionBSF, RelationSF
	float ptfactor, pbfactor, rfactor;

	//SymbolSF, InsPenalty
	float qfactor, InsPen;

	std::shared_ptr<TableCYK> parseCore(std::shared_ptr<Sample> &M)
	{
		//Compute the normalized size of a symbol for sample M
		M->detRefSymbol();

		M->computeSegDistance(M->RX, M->RY);

		int N = M->getSegUnitSize();
		int K = pG->noTerminales.size();

		//Cocke-Younger-Kasami (CYK) algorithm for 2D-SCFG
		tcyk = std::make_shared<TableCYK>(N, K);

		std::cout << "CYK table initialization:" << std::endl;
		initCYKterms(M, *tcyk, N, K);

		std::vector<std::shared_ptr<LogSpace>> logspace(N + 1);
		std::list<std::shared_ptr<CellCYK>> c1setH, c1setV, c1setU, c1setI, c1setM, c1setS;
		//SpaRel SPR();
		pSPR = std::make_shared<SpaRel>(pGMM, M);

		//Init spatial space for size 1
		logspace[1] = std::make_shared<LogSpace>(tcyk->get(1), tcyk->size(1), M->RX, M->RY);

		std::cout << "\nCYK parsing algorithm" << std::endl;
		std::cout << "Size 1: Generated " << tcyk->size(1) << std::endl;

		//CYK algorithm main loop
		for (int talla = 2; talla <= N; talla++)
		{
			for (int a = 1; a < talla; a++)
			{
				int b = talla - a;
				for (std::shared_ptr<CellCYK> c1 = tcyk->get(a); c1.use_count(); c1 = c1->sig)
				{
					//Clear lists
					c1setH.clear();
					c1setV.clear();
					c1setU.clear();
					c1setI.clear();
					c1setM.clear();
					c1setS.clear();


					//Get the subset of regions close to c1 according to different spatial relations
					logspace[b]->getH(c1, c1setH); //Horizontal (right)
					logspace[b]->getV(c1, c1setV); //Vertical (down)
					logspace[b]->getU(c1, c1setU); //Vertical (up)
					logspace[b]->getI(c1, c1setI); //Inside (sqrt)
					logspace[b]->getM(c1, c1setM); //mroot (sqrt[i])

#ifdef SHOW_PIPELINE
					int time = 1;
					cv::Mat tmpImg1 = M->getRGBImg(), tmpImg2;
					drawCellWithColor(tmpImg1, c1, cv::Scalar(0, 0, 255));
					int winx = 10, winy = 10;
					tmpImg2 = tmpImg1.clone();
					cv::imshow("Current cell", tmpImg2);
					cv::moveWindow("Current cell", winx, winy);
					cv::waitKey(time);
#endif // SHOW_PIPELINE

					for (auto c2 = c1setH.begin(); c2 != c1setH.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE

						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsH, Grammar::H);
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsSup, Grammar::SUP);
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsSub, Grammar::SUB);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setH traverse", tmpImg2);
					cv::moveWindow("c1setH traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setV.begin(); c2 != c1setV.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsV, Grammar::V);
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsVe, Grammar::VE);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setV traverse", tmpImg2);
					cv::moveWindow("c1setV traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setU.begin(); c2 != c1setU.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(*c2, c1, M, *tcyk, talla, pG->prodsV, Grammar::V);
						traverseProductionB(*c2, c1, M, *tcyk, talla, pG->prodsVe, Grammar::VE);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setU traverse", tmpImg2);
					cv::moveWindow("c1setU traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setI.begin(); c2 != c1setI.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsIns, Grammar::INS);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setI traverse", tmpImg2);
					cv::moveWindow("c1setI traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setM.begin(); c2 != c1setM.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, *tcyk, talla, pG->prodsMrt, Grammar::MRT);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setM traverse", tmpImg2);
					cv::moveWindow("c1setM traverse", winx, winy);
					cv::waitKey(time);
#endif // SHOW_PIPELINE

					//Look for combining {x_subs} y {x^sups} in {x_subs^sups}
					//TODO :
					for (int pps = 0; pps<c1->nnt; pps++) {

						//If c1->noterm[pa] is a Hypothesis of a subscript (parent_son)
						if (c1->vNoTerm[pps] && c1->vNoTerm[pps]->prod && c1->vNoTerm[pps]->prod->tipo() == 'B') {

							logspace[b + c1->vNoTerm[pps]->hleft->pCInfo->talla]->getS(c1, c1setS); //sup/sub-scripts union

							for (auto c2 = c1setS.begin(); c2 != c1setS.end(); c2++) {

								if ((*c2)->pCInfo->box.x == c1->pCInfo->box.x && c1 != *c2) {

									traverseProductionSSE(c1, *c2, M, *tcyk, talla, pG->prodsSSE);
								}

							}//end for c2 in c1setS

							c1setS.clear();
						}
					}//end for(int pps=0; pps<c1->nnt; pps++)


				}//for (std::shared_ptr<CellCYK> c1 = tcyk.get(a); c1.use_count(); c1 = c1->sig)
			}//for (int a = 1; a < talla; a++)

			std::cout << "Size " << talla << ": Generated " << tcyk->size(talla) << std::endl;

			if (talla <= N)
			{
				//Create new logspace structure of size "talla"
				logspace[talla] = std::make_shared<LogSpace>(tcyk->get(talla), tcyk->size(talla), M->RX, M->RY);

#ifdef SHOW_PIPELINE
				int time = 0;
				cv::Mat tmpImg1 = M->getRGBImg(), tmpImg2;
				int winx = 10, winy = 10;
				for (std::shared_ptr<CellCYK> c = tcyk.get(talla); c.use_count(); c = c->sig)
				{
					drawCellWithColor(tmpImg1, c, RandomColor());
				}
				std::stringstream ioStr;
				ioStr << talla;
				cv::imshow("talla" + ioStr.str(), tmpImg1);
				cv::moveWindow("talla" + ioStr.str(), winx, winy);
				cv::waitKey(time);
#endif // SHOW_PIPELINE
			}


		}//for (int talla = 2; talla <= N; talla++)

		return tcyk;
	}

	void initCYKterms(std::shared_ptr<Sample> &M, TableCYK &tcyk, int N, int K)
	{
		for (int i = 0; i<M->getSegUnitSize(); i++) {

			int cmy, asc, des;

			std::cout << "Segment " << i << std::endl;

			std::vector<int> clase;
			std::vector<float> pr;

			M->getSegUnitInfo(i, NB, clase, pr, cmy, asc, des);

			std::shared_ptr<CellCYK> pCell = std::make_shared<CellCYK>(K, N);

			pCell->setRegion(*M, i);

			int actualNB = clase.size();

			bool insertar = false;
			for (size_t i = 0; i < pG->prodTerms.size(); i++)
			{
				std::shared_ptr<ProductionT> &prod = pG->prodTerms[i];
				int ntID = prod->getNoTerm();

				for (int k = 0; k < actualNB; k++)
				{
					if (pr[k] > 0.0 && prod->getClase(clase[k]) && prod->getPrior(clase[k]) > -FLT_MAX)
					{
						//the origin probability
						/*float prob = log(InsPen)
						+ ptfactor * prod->getPrior(clase[k])
						+ qfactor  * log(pr[k])
						+ dfactor  * log(duration->prob(clase[k], 1));*/
						float prob = log(InsPen)
							+ ptfactor * prod->getPrior(clase[k])
							+ qfactor * log(pr[k]);

						if (pCell->vNoTerm[ntID].use_count())
						{
							if (pCell->vNoTerm[ntID]->pr > prob + prod->getPrior(clase[k]))
								continue;
							else
								pCell->vNoTerm[ntID].reset();
						}

						insertar = true;

						//Create new symbol
						pCell->vNoTerm[ntID] = std::make_shared<Hypothesis>(clase[k], prob, pCell->pCInfo);
						pCell->vNoTerm[ntID]->pt = prod;

						//Compute the vertical centroid according to the type of symbol
						StdSymInfo sSymInfo = pSymSet->stdInfoClase(clase[k]);
						double rel_h = sSymInfo.rel_t - sSymInfo.rel_y;
						double ratio = (pCell->pCInfo->box.t - pCell->pCInfo->box.y) / rel_h;
						double cen = pCell->pCInfo->box.y + (STD_CENTER_Y - sSymInfo.rel_y) * ratio;
						//Vertical center

						pCell->vNoTerm[ntID]->lcen = cen;
						pCell->vNoTerm[ntID]->rcen = cen;

						pCell->vNoTerm[ntID]->lineTop = pCell->pCInfo->box.y;
						pCell->vNoTerm[ntID]->lineBottom = pCell->pCInfo->box.t;
						pCell->vNoTerm[ntID]->totalSymWidth = pCell->pCInfo->box.s - pCell->pCInfo->box.x;
					}
				}
			}

			if (insertar)
			{
				for (int j = 0; j < K; j++) {
					if (pCell->vNoTerm[j].use_count()) {
						std::cout << pSymSet->strClase(pCell->vNoTerm[j]->clase) << " "
							<< pG->key2str(j) << " " << exp(pCell->vNoTerm[j]->pr) << std::endl;
						/*printf("%12s [%s] %g\n", pSymSet->strClase(pCell->vNoTerm[j]->clase),
						pG->key2str(j), exp(pCell->vNoTerm[j]->pr));*/
					}
				}

				//Add to parsing table (size=1)
				tcyk.add(1, pCell, -1, pG->esInit);
			}
		}
	}

	std::shared_ptr<CellCYK> fusion(std::shared_ptr<Sample> &M, std::shared_ptr<ProductionB>& pd, std::shared_ptr<CellCYK>& pCA,
									int ntIDA, std::shared_ptr<CellCYK>& pCB, int ntIDB, double prob)
	{
		std::shared_ptr<CellCYK> pCS;

		if (!pCA->compatible(pCB) || pd->prior == -FLT_MAX)
			return pCS;

		//Penalty according to distance between Cell
		float grpen = 1.0;
		//TODO : induce the penalty function
		if (clusterF > 0.0) {

			grpen = M->group_penalty(pCA->pCInfo->segMask, pCB->pCInfo->segMask);
			//If distance is infinity -> not visible
			if (grpen >= M->INF_DIST)
				return pCS;

			/* float cax = (pCA->pCInfo->box.x + pCA->pCInfo->box.s)*0.5;
			float cay = (pCA->pCInfo->box.y + pCA->pCInfo->box.t)*0.5;
			float cbx = (pCB->pCInfo->box.x + pCB->pCInfo->box.s)*0.5;
			float cby = (pCB->pCInfo->box.y + pCB->pCInfo->box.t)*0.5;

			grpen = std::sqrt((cbx - cax)*(cbx - cax) + (cby - cay)*(cby - cay)) / M->NORMF;*/

			//Compute penalty
			grpen = std::max(-0.9f, grpen);
			grpen = 1.0 / (1.0 + grpen);
			grpen = pow(grpen, clusterF);
			//grpen = 1 * grpen;
			//std::cout << "grpen = " << grpen << std::endl;
		}

		//Get nonterminal
		int ps = pd->S;
		int N = M->getSegUnitSize();
		int K = pG->noTerminales.size();

		pCS = std::make_shared<CellCYK>(K, N);

		std::shared_ptr<Hypothesis> &A = pCA->vNoTerm[ntIDA], &B = pCB->vNoTerm[ntIDB];
		int tallaA = A->pCInfo->talla;
		int tallaB = B->pCInfo->talla;
		int tallaS = tallaA + tallaB;

		//Compute the (log)probability
		//prob = pbfactor * pd->prior + rfactor * log(prob * grpen) + A->pr + B->pr;
		double tmpProb = prob;
		double pdPrior = pd->prior;
		if (pdPrior > 0.0)
		{
			char typeC = pd->tipo();
			//char allType[7] = { 'H', 'B', 'P', 'V', 'e', 'I', 'M' };
			char allType[5] = { 'H', 'B', 'P', 'V', 'e' };
			for (size_t i = 0; i < 5; i++)
			{
				if (typeC == allType[i]) continue;
				double cdpr = pSPR->getWithType(A, B, allType[i]);
				if (cdpr > prob)
				{
					pdPrior = 0.0;
					break;
				}
			}
		}


		prob = pbfactor * pdPrior + rfactor * log(prob) + dfactor * log(grpen)/* * (1.0 / tallaS)*/ + A->pr + B->pr;

#ifdef LOG_COMPUTE
		HL_GENERAL_LOG("Production: " + pd->sS << "(" << tallaS << ") ---> " + pd->sA << "(" << tallaA << ") " + pd->sB << "(" << tallaB << ")");
		HL_GENERAL_LOG("(" << prob << ") " << "dist penalty: " << grpen << " ; A prob: " << A->pr << " ; B prob : " << B->pr);
		std::cout << std::endl;
#endif // LOG_COMPUTE

		//Copute resulting region
		pCS->pCInfo->box.x = std::min(pCA->pCInfo->box.x, pCB->pCInfo->box.x);
		pCS->pCInfo->box.y = std::min(pCA->pCInfo->box.y, pCB->pCInfo->box.y);
		pCS->pCInfo->box.s = std::max(pCA->pCInfo->box.s, pCB->pCInfo->box.s);
		pCS->pCInfo->box.t = std::max(pCA->pCInfo->box.t, pCB->pCInfo->box.t);

		//Set the Cell covered
		pCS->ccUnion(*pCA, *pCB);

		int clase = -1;
		if (!pd->check_out() && pSymSet->checkClase(pd->get_outstr()))
			clase = pSymSet->keyClase(pd->get_outstr());

		//Create hypothesis
		pCS->vNoTerm[ps] = std::make_shared<Hypothesis>(clase, prob, pCS->pCInfo);

		MergeRegionsCenter(pCA, ntIDA, pCB, ntIDB, pCS, ps, pd->merge_cen);

		pCS->vNoTerm[ps]->hleft = A;
		pCS->vNoTerm[ps]->hright = B;
		pCS->vNoTerm[ps]->prod = pd;

		if (clase >= 0)
		{
			for (size_t i = 0; i < pG->prodTerms.size(); i++)
			{
				std::shared_ptr<ProductionT> &prod = pG->prodTerms[i];

				if (prod->getClase(clase) && prod->getPrior(clase) > -FLT_MAX)
				{
					pCS->vNoTerm[ps]->pt = prod;
					break;
				}
			}
		}

		return pCS;
	}


	void traverseProductionB(std::shared_ptr<CellCYK> &c1,
							 std::shared_ptr<CellCYK> &c2,
							 std::shared_ptr<Sample> &M, TableCYK &tcyk, int talla,
							 std::vector<std::shared_ptr<ProductionB>> &vProds, Grammar::PBTYPE pType)
	{
		for (size_t i = 0; i < vProds.size(); i++)
		{
			std::shared_ptr<ProductionB> &pd = vProds[i];

			if (pd->prior == -FLT_MAX) continue;
			//Production S -> A B
			int ps = pd->S;
			int pa = pd->A;
			int pb = pd->B;


			if (c1->vNoTerm[pa].use_count() && c2->vNoTerm[pb].use_count()) {
				std::shared_ptr<Hypothesis> &ha = c1->vNoTerm[pa], hb = c2->vNoTerm[pb];

				double cdpr = 0.0;

				switch (pType)
				{
				case Grammar::H:
					cdpr = pSPR->getHorProb(ha, hb);
					break;
				case Grammar::SUP:
					cdpr = pSPR->getSupProb(ha, hb);
					break;
				case Grammar::SUB:
					cdpr = pSPR->getSubProb(ha, hb);
					break;
				case Grammar::V:
					cdpr = pSPR->getVerProb(ha, hb);
					break;
				case Grammar::VE:
					cdpr = pSPR->getVerProb(ha, hb, true);
					break;
				case Grammar::INS:
					cdpr = pSPR->getInsProb(ha, hb);
					break;
				case Grammar::MRT:
					cdpr = pSPR->getMrtProb(ha, hb);
					break;
				case Grammar::SSE:
					//cdpr = pSPR->getSupProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					break;
				default:
					break;
				}
				if (cdpr <= 0.0) continue;

#ifdef LOG_COMPUTE
				std::string typeStr = pG->strType(pType);
				std::cout << "The space type : " << typeStr << "	; space prob : " << cdpr << std::endl;
#endif // LOG_COMPUTE

				std::shared_ptr<CellCYK> pCell = fusion(M, pd, c1, pa, c2, pb, cdpr);

				if (!pCell.use_count()) continue;

				if (pCell->vNoTerm[ps].use_count()) {
					tcyk.add(talla, pCell, ps, pG->esInit); //Add to parsing table (size=talla)
				}
				else {
					tcyk.add(talla, pCell, -1, pG->esInit); //Add to parsing table
				}
			}
		}
	}

	void traverseProductionSSE(std::shared_ptr<CellCYK> &c1,
							   std::shared_ptr<CellCYK> &c2,
							   std::shared_ptr<Sample> &M, TableCYK &tcyk, int talla,
							   std::vector<std::shared_ptr<ProductionB>> &vProdSSE)
	{
		int N = M->getSegUnitSize();
		int K = pG->noTerminales.size();

		for (size_t i = 0; i < vProdSSE.size(); i++)
		{
			std::shared_ptr<ProductionB> &pd = vProdSSE[i];

			if (pd->prior == -FLT_MAX) continue;
			//Production S -> A B
			int ps = pd->S;
			int pa = pd->A;
			int pb = pd->B;

			if (c1->vNoTerm[pa] && c2->vNoTerm[pb]
				&& c1->vNoTerm[pa]->prod && c2->vNoTerm[pb]->prod
				&& c1->vNoTerm[pa]->hleft == c2->vNoTerm[pb]->hleft
				&& c1->vNoTerm[pa]->prod->tipo() == 'B'
				&& c2->vNoTerm[pb]->prod->tipo() == 'P'
				&& c1->vNoTerm[pa]->hright->compatible(c2->vNoTerm[pb]->hright)) {


				//Subscript and superscript should start almost vertically aligned
				if (std::abs(c1->vNoTerm[pa]->hright->pCInfo->box.x - c2->vNoTerm[pb]->hright->pCInfo->box.x) > 3 * M->RX) continue;
				//Subscript and superscript should not overlap
				if (std::max(solape(c1->vNoTerm[pa]->hright, c2->vNoTerm[pb]->hright),
							 solape(c2->vNoTerm[pb]->hright, c1->vNoTerm[pa]->hright)) > 0.1) continue;

				//TODO : add the score strategy for SSE situation
				float prob = c1->vNoTerm[pa]->pr + c2->vNoTerm[pb]->pr - c1->vNoTerm[pa]->hleft->pr;

				std::shared_ptr<CellCYK> pCell = std::make_shared<CellCYK>(K, N);

				pCell->pCInfo->box.x = std::min(c1->pCInfo->box.x, c2->pCInfo->box.x);
				pCell->pCInfo->box.y = std::min(c1->pCInfo->box.y, c2->pCInfo->box.y);
				pCell->pCInfo->box.s = std::max(c1->pCInfo->box.s, c2->pCInfo->box.s);
				pCell->pCInfo->box.t = std::max(c1->pCInfo->box.t, c2->pCInfo->box.t);

				pCell->vNoTerm[ps] = std::make_shared<Hypothesis>(-1, prob, pCell->pCInfo);

				pCell->vNoTerm[ps]->lcen = c1->vNoTerm[pa]->lcen;
				pCell->vNoTerm[ps]->rcen = c1->vNoTerm[pa]->rcen;
				pCell->ccUnion(*c1, *c2);

				pCell->vNoTerm[ps]->hleft = c1->vNoTerm[pa];
				pCell->vNoTerm[ps]->hright = c2->vNoTerm[pb]->hright;
				pCell->vNoTerm[ps]->prod = pd;
				//Save the production of the superscript in order to recover it when printing the used productions
				pCell->vNoTerm[ps]->prod_sse = c2->vNoTerm[pb]->prod;

				MergeRegionsCenter(c1, pa, c2, pb, pCell, ps, pd->merge_cen);

				tcyk.add(talla, pCell, ps, pG->esInit);
			}
		}
	}

};



