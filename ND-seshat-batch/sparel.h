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
#ifndef _SPAREL_
#define _SPAREL_

class CellCYK;

#include <cstdio>
#include "hypothesis.h"
#include "cellcyk.h"
#include "gmm.h"
#include "sample.h"

//Aux functions
inline std::shared_ptr<Hypothesis> leftmost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count())
		return h;

	std::shared_ptr<Hypothesis> izq = leftmost(h->hleft);
	std::shared_ptr<Hypothesis> der = leftmost(h->hright);

	return izq->pCInfo->box.x < der->pCInfo->box.x ? izq : der;
}

inline std::shared_ptr<Hypothesis> rightmost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count())
		return h;

	std::shared_ptr<Hypothesis> izq = rightmost(h->hleft);
	std::shared_ptr<Hypothesis> der = rightmost(h->hright);

	return izq->pCInfo->box.s > der->pCInfo->box.s ? izq : der;
}

//Spcecial with frac term
inline std::shared_ptr<Hypothesis> topmost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count() || (h->prod.use_count() && h->prod->tipo() == 'V'))
		return h;

	std::shared_ptr<Hypothesis> izq = topmost(h->hleft);
	std::shared_ptr<Hypothesis> der = topmost(h->hright);

	return izq->pCInfo->box.y < der->pCInfo->box.y ? izq : der;
}

//Spcecial with frac term
inline std::shared_ptr<Hypothesis> buttommost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count() || (h->prod.use_count() && h->prod->tipo() == 'V'))
		return h;

	std::shared_ptr<Hypothesis> izq = buttommost(h->hleft);
	std::shared_ptr<Hypothesis> der = buttommost(h->hright);

	return izq->pCInfo->box.t > der->pCInfo->box.t ? izq : der;
}


//Percentage of the area of region A that overlaps with region B
inline float solape(std::shared_ptr<Hypothesis> &a, std::shared_ptr<Hypothesis> &b) {
	int x = std::max(a->pCInfo->box.x, b->pCInfo->box.x);
	int y = std::max(a->pCInfo->box.y, b->pCInfo->box.y);
	int s = std::min(a->pCInfo->box.s, b->pCInfo->box.s);
	int t = std::min(a->pCInfo->box.t, b->pCInfo->box.t);

	if (s >= x && t >= y) {
		float aSolap = (s - x + 1.0)*(t - y + 1.0);
		float aTotal = (a->pCInfo->box.s - a->pCInfo->box.x + 1.0)*(a->pCInfo->box.t - a->pCInfo->box.y + 1.0);

		return aSolap / aTotal;
	}

	return 0.0;
}

inline bool checkLeftRight(std::shared_ptr<Hypothesis> &h1, std::shared_ptr<Hypothesis> &h2)
{
	//Check left-to-right order constraint in Hor/Sub/Sup relationships
	std::shared_ptr<Hypothesis> rma = rightmost(h1);
	std::shared_ptr<Hypothesis> lmb = leftmost(h2);

	int rmaw = rma->pCInfo->box.s - rma->pCInfo->box.x;
	int lmbw = lmb->pCInfo->box.s - lmb->pCInfo->box.x;
	int rminshift = rmaw * 0.2;

	int sqrtID = 58;
	if (rma->clase == sqrtID)
	{
		rminshift = rmaw - lmbw * 0.5;
	}

	if (lmb->pCInfo->box.x < rma->pCInfo->box.x + rminshift || lmb->pCInfo->box.s + rminshift <= rma->pCInfo->box.s)
		return false;
	else
		return true;
}

//Check the upon symbols
inline bool checkVerticalUpon(std::shared_ptr<Hypothesis> &ha, coo &bbox, double bcen)
{
	if (!ha.use_count())
		HL_CERR("Invalid Pointer");

	if (ha->pt.use_count())
	{
		coo &abox = ha->pCInfo->box;
		int awidth = abox.s - abox.x, aheight = abox.t - abox.y;

		if ((abox.x < bbox.x && (bbox.x - abox.x) > 0.5 * awidth) ||
			(abox.s > bbox.s && (abox.s - bbox.s) > 0.5 * awidth))
		{
			double downRatio = double(abox.t - bcen) / aheight;
			if (downRatio > 0.3) return false;
		}

		return true;
	}

	bool lcheck = checkVerticalUpon(ha->hleft, bbox, bcen);
	bool rcheck = checkVerticalUpon(ha->hright, bbox, bcen);

	return lcheck && rcheck;
}

//Check the down symbols
inline bool checkVerticalDown(std::shared_ptr<Hypothesis> &hb, coo &abox, double acen)
{
	if (!hb.use_count())
		HL_CERR("Invalid Pointer");

	if (hb->pt.use_count())
	{
		coo &bbox = hb->pCInfo->box;
		int bwidth = bbox.s - bbox.x, bheight = bbox.t - bbox.y;

		if ((bbox.x < abox.x && (abox.x - bbox.x) > 0.5 * bwidth) ||
			(bbox.s > abox.s && (bbox.s - abox.s) > 0.5 * bwidth))
		{
			double upRatio = double(acen - bbox.y) / bheight;
			if (upRatio > 0.3) return false;
		}

		return true;
	}

	bool lcheck = checkVerticalDown(hb->hleft, abox, acen);
	bool rcheck = checkVerticalDown(hb->hright, abox, acen);

	return lcheck && rcheck;
}

inline bool isDigit(std::shared_ptr<Hypothesis> &H)
{
	return H->clase >= 3 && H->clase <= 11;
}

class SpaRel {
public:
	const int NRELS = 6;
	const int NFEAT = 9;

private:
	std::shared_ptr<GMM> model;
	std::shared_ptr<Sample> mue;
	std::vector<float> probs;

	double compute_prob(std::shared_ptr<Hypothesis> &h1, std::shared_ptr<Hypothesis> &h2, int k)
	{
		//Set probabilities according to spatial constraints  
		double result = 1.0;
		if (k <= 2) {
			//Check left-to-right order constraint in Hor/Sub/Sup relationships
			std::shared_ptr<Hypothesis> rma = rightmost(h1);
			std::shared_ptr<Hypothesis> lmb = leftmost(h2);

			int rmaw = rma->pCInfo->box.s - rma->pCInfo->box.x;
			int lmbw = lmb->pCInfo->box.s - lmb->pCInfo->box.x;
			//int minshift = 0;
			int rminshift = rmaw * 0.2;
			int lminshift = lmbw * 0.2;

			if (lmb->pCInfo->box.x < rma->pCInfo->box.x + rminshift || lmb->pCInfo->box.s + lminshift <= rma->pCInfo->box.s)
				return 0.0;

			if (k == 0)
			{
				double overlap = std::max(solape(rma, lmb), solape(lmb, rma));
				result = 1.0 / (1.0 + exp((overlap - 0.85) * 10));
			}

			return 1.0;
		}

		//Compute probabilities
		std::vector<float> sample(NFEAT);

		getFeas(h1, h2, sample, mue->RY);

		//Get spatial relationships probability from the model
		model->posterior(&sample[0], &probs[0]);

		//Slightly smooth probabilities because GMM classifier can provide
		//to biased probabilities. Thsi way we give some room to the
		//language model (the 2D-SCFG grammar)
		smooth(probs);

		return probs[k];
	}

	void smooth(std::vector<float> &post)
	{
		for (int i = 0; i<NRELS; i++)
			post[i] = (post[i] + 0.02) / (1.00 + NRELS*0.02);
	}

public:
	SpaRel(std::shared_ptr<GMM> &gmm, std::shared_ptr<Sample> &m)
	{
		model = gmm;
		mue = m;
		probs.resize(NRELS);
	}

	~SpaRel() {}

	void getFeas(std::shared_ptr<Hypothesis> &a, std::shared_ptr<Hypothesis> &b, std::vector<float> &sample, int ry)
	{
		//Normalization factor: combined height
		float F = std::max(a->pCInfo->box.t, b->pCInfo->box.t) - std::min(a->pCInfo->box.y, b->pCInfo->box.y) + 1;

		sample[0] = (b->pCInfo->box.t - b->pCInfo->box.y + 1) / F;
		sample[1] = (a->rcen - b->lcen) / F;
		sample[2] = ((a->pCInfo->box.s + a->pCInfo->box.x) / 2.0 - (b->pCInfo->box.s + b->pCInfo->box.x) / 2.0) / F;
		sample[3] = (b->pCInfo->box.x - a->pCInfo->box.s) / F;
		sample[4] = (b->pCInfo->box.x - a->pCInfo->box.x) / F;
		sample[5] = (b->pCInfo->box.s - a->pCInfo->box.s) / F;
		sample[6] = (b->pCInfo->box.y - a->pCInfo->box.t) / F;
		sample[7] = (b->pCInfo->box.y - a->pCInfo->box.y) / F;
		sample[8] = (b->pCInfo->box.t - a->pCInfo->box.t) / F;
	}

	double getHorProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb, double maxScore = 1.0)
	{
		coo &abox = ha->pCInfo->box, &bbox = hb->pCInfo->box;
		int awidth = abox.s - abox.x, aheight = abox.t - abox.y;
		int bwidth = bbox.s - bbox.x, bheight = bbox.t - bbox.y;

		//Exclude the obviously invalid situation
		if (abox.t <= bbox.y || abox.y >= bbox.t || !checkLeftRight(ha, hb)) return 0.0;

		double score = 0.0;

		int cenDiff = hb->lcen - ha->rcen;
		int tDiff = hb->lineTop - ha->lineTop;
		int bDiff = hb->lineBottom - ha->lineBottom;

		double avgH = 0.5 * (aheight + bheight);

		double cenRatio = cenDiff / avgH;

		double cenScore = 1.0 - std::abs(cenRatio);
		score = cenScore;
		return std::min(maxScore, score);
	}

	double getSubProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb, double maxScore = 1.0)
	{
		coo &abox = ha->pCInfo->box, &bbox = hb->pCInfo->box;
		int awidth = abox.s - abox.x, aheight = abox.t - abox.y;
		int bwidth = bbox.s - bbox.x, bheight = bbox.t - bbox.y;

		//Exclude the obviously invalid situation
		if (abox.y > bbox.y || !checkLeftRight(ha, hb)) return 0.0;

		double score = 0.0;
		std::shared_ptr<Hypothesis> topb = topmost(hb);

		int cenDiff = topb->lcen - ha->rcen;

		double avgH = 0.5 * (aheight + topb->pCInfo->box.t - topb->pCInfo->box.y);

		double cenRatio = cenDiff / avgH;
		double cenScore = 1 / (1 + exp((-cenRatio + 0.1) * 15));;

		score = cenScore;

		// Add the penalty for distance in SUB
		double xsDiff = double(bbox.x - abox.s) - mue->RX;
		double ytDiff = double(bbox.y - abox.t) - mue->RY;
		double penRatio = std::max(xsDiff, ytDiff) / avgH;
		double verticalPen = 0.0;
		if (penRatio > 0)
			verticalPen = 1.0 / (1.0 + exp(-(penRatio - 0.2) * 15));

		score -= (verticalPen * 0.5);

		return std::min(maxScore, score);
	}

	double getSupProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb, double maxScore = 1.0)
	{
		coo &abox = ha->pCInfo->box, &bbox = hb->pCInfo->box;
		int awidth = abox.s - abox.x, aheight = abox.t - abox.y;
		int bwidth = bbox.s - bbox.x, bheight = bbox.t - bbox.y;

		//Exclude the obviously invalid situation
		if (abox.t < bbox.t || !checkLeftRight(ha, hb)) return 0.0;

		double score = 0.0;
		std::shared_ptr<Hypothesis> bottomb = buttommost(hb);

		int cenDiff = ha->rcen - bottomb->lcen;

		double avgH = 0.5 * (aheight + bottomb->pCInfo->box.t - bottomb->pCInfo->box.y);

		double cenRatio = cenDiff / avgH;
		double cenScore = 1 / (1 + exp((-cenRatio + 0.1) * 15));;

		score = cenScore;

		// Add the penalty for distance in SUP
		double xsDiff = double(bbox.x - abox.s) - mue->RX;
		double ytDiff = double(abox.y - bbox.t) - mue->RY;
		double penRatio = std::max(xsDiff, ytDiff) / avgH;
		double verticalPen = 0.0;
		if (penRatio > 0)
			verticalPen = 1.0 / (1.0 + exp(-(penRatio - 0.2) * 15));

		score -= (verticalPen * 0.5);

		return std::min(maxScore, score);
	}

	double getVerProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb, bool strict = false)
	{
		coo &abox = ha->pCInfo->box, &bbox = hb->pCInfo->box;

		//Pruning
		if (bbox.y < (abox.y + abox.t) / 2
			|| abs((abox.x + abox.s) / 2 - (bbox.x + bbox.s) / 2) > 2.5*mue->RX
			|| (bbox.x > abox.s || bbox.s < abox.x))
			return 0.0;

		//Check the position relationship of outer lines
		double bcen = (hb->lcen + hb->rcen)*0.5;
		double acen = (ha->lcen + ha->rcen)*0.5;
		if (!checkVerticalUpon(ha, bbox, bcen) || !checkVerticalDown(hb, abox, acen))
			return 0.0;

		double score = 0.0;

		score = compute_prob(ha, hb, 3);

		if (strict)
		{
			//Penalty for strict relationships
			float penalty = abs(abox.x - bbox.x) / (3.0*mue->RX)
				+ abs(abox.s - bbox.s) / (3.0*mue->RX);

			if (penalty > 0.95) penalty = 0.95;

			score *= (1 - penalty);
		}

		return score;
	}

	double getInsProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		int sqrtH = hb->pCInfo->box.t - hb->pCInfo->box.y;

		if (solape(hb, ha) < 0.4 ||
			hb->pCInfo->box.x < ha->pCInfo->box.x || hb->pCInfo->box.y < (ha->pCInfo->box.y - sqrtH*0.3))
			return 0.0;


		std::shared_ptr<Hypothesis> rma = rightmost(ha);
		std::shared_ptr<Hypothesis> rmb = rightmost(hb);

		if (rma->pCInfo->box.s <= rmb->pCInfo->box.x + 1) return 0.0;

		/*if (solape(hb, ha) < 0.5 ||
		hb->pCInfo->box.x < ha->pCInfo->box.x || hb->pCInfo->box.y < ha->pCInfo->box.y)
		return 0.0;*/

		return compute_prob(ha, hb, 4);
	}

	double getMrtProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		return compute_prob(ha, hb, 5);
	}
};

#endif
