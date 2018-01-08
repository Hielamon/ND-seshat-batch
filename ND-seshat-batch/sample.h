#pragma once
#include <OpencvCommon.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <commonMacro.h>
#include <algorithm>

#include "symSet.h"


//The unit result of segmentation or symbol recognition procedure 
struct SegUnit
{
	cv::Rect ROI;
	double score;
	int symID;
	std::string symStr;
};

inline int computeRectArea(cv::Rect r)
{
	return (r.width > 0 && r.height > 0) ? r.width * r.height : std::max(r.width, r.height);
}

//The Class include the Handwritten formula image
//and the symbol recognition result
class Sample
{
public:
	Sample(const std::shared_ptr<SymSet> &pSymSet_) : pSymSet(pSymSet_) {}
	~Sample() {}

	//Load the Sample information from the annotation file of VOC 2007
	bool LoadFromVOC2007XML(const std::string &xmlName, const std::string &imgPath,
							std::string &gtLatexID, const std::string &charMapName = "VOC2007/charmap_.txt")
	{
		if (!vSegUnits.empty()) vSegUnits.clear();

		std::fstream fs(charMapName, std::ios::in);
		if (!fs.is_open())
			HL_CERR_RETURN_FALSE("Failed to open the file " + charMapName);
		int charMapSize;
		fs >> charMapSize;

		std::vector<std::pair<std::string, int>> mCharMap;
		for (size_t i = 0; i < charMapSize; i++)
		{
			std::string charName;
			int mapID;
			fs >> charName >> mapID;

			mCharMap.push_back(std::pair<std::string, int>(charName, mapID));
		}
		fs.close();

		boost::property_tree::ptree pt, annotation;
		boost::property_tree::read_xml(xmlName, pt);
		if (pt.empty())
			HL_CERR_RETURN_FALSE("Cannot load the xml tree from " + xmlName);

		std::string imgName, latexIDStr;
		annotation = pt.get_child("annotation");
		for (auto iter = annotation.begin(); iter != annotation.end(); iter++)
		{
			if (iter->first == "filename")
			{
				imgName = iter->second.data();
				Img = cv::imread(imgPath + imgName, cv::IMREAD_GRAYSCALE);
				cv::threshold(Img, Img, 100, 255, cv::ThresholdTypes::THRESH_BINARY);
			}
			else if (iter->first == "latex")
			{
				latexIDStr = iter->second.data();
				latexIDStr.assign(latexIDStr.begin() + 1, latexIDStr.end() - 1);
			}
			else if (iter->first == "size")
			{
				W = iter->second.get<int>("width");
				H = iter->second.get<int>("height");
			}
			else if (iter->first == "object")
			{
				SegUnit seg;
				int nameID = iter->second.get<int>("name");
				if (nameID >= 0 && nameID < mCharMap.size())
				{
					seg.symID = mCharMap[nameID].second;
					if (seg.symID < 0 || seg.symID >= pSymSet->getNClases())
						HL_CERR_RETURN_FALSE("The symbol " + mCharMap[nameID].first +
								" doesn't have a valid ID (" << seg.symID << ") in the symbol set");

					seg.symStr = pSymSet->strClase(seg.symID);
				}
				else
					HL_CERR_RETURN_FALSE("The symbol ID " << nameID << " from xml file is not valid ");

				seg.score = 1.0;

				int x = iter->second.get<int>("bndbox.xmin");
				int y = iter->second.get<int>("bndbox.ymin");
				int s = iter->second.get<int>("bndbox.xmax");
				int t = iter->second.get<int>("bndbox.ymax");
				seg.ROI = cv::Rect(x, y, s - x, t - y);

				vSegUnits.push_back(seg);
			}
		}

		gtLatexID = latexIDStr;
		return true;
	}

	void ShowSample(const std::string &windowName = "Sample")
	{
		if (Img.empty())
			HL_CERR("There is not a image loaded");

		double megapix = std::max(0.5 * 1e6, double(W*H));
		double scale = std::sqrt(megapix / (W * H));
		cv::Mat showImg;
		cv::cvtColor(Img, showImg, cv::COLOR_GRAY2BGR);
		cv::resize(showImg, showImg, cv::Size(W*scale, H*scale));

		std::stringstream ioStr;
		int i = 0;
		std::for_each(vSegUnits.begin(), vSegUnits.end(), [&](SegUnit &seg)
		{
			cv::Scalar color = RandomColor();
			cv::Rect scaledROI(seg.ROI.x * scale, seg.ROI.y * scale,
							   seg.ROI.width * scale, seg.ROI.height * scale);
			cv::rectangle(showImg, scaledROI, color, 2 * scale, cv::LINE_AA);
			ioStr.str("");
			//ioStr << seg.symID << " | " << seg.symStr << " | " << seg.score;
			ioStr << seg.symID << " | " << seg.symStr << " | " << i;
			i++;
			cv::putText(showImg, ioStr.str(), scaledROI.tl() - cv::Point(0, 2 * scale), cv::HersheyFonts::FONT_HERSHEY_COMPLEX, 0.8, color);
		});
		resizeShow(windowName, showImg);
		cv::waitKey(0);
	}

	void  detRefSymbol()
	{
		std::vector<int> vmedx, vmedy;
		int nregs = 0, lAr;
		float mAr = 0;
		RX = 0, RY = 0;

		//Compute reference symbol for normalization
		for (int i = 0; i<vSegUnits.size(); i++) {
			int ancho = vSegUnits[i].ROI.width + 1;
			int alto = vSegUnits[i].ROI.height + 1;
			float aspectratio = (float)ancho / alto;
			int area = ancho*alto;

			vmedx.push_back(ancho);
			vmedy.push_back(alto);

			mAr += area;
			if (aspectratio >= 0.25 && aspectratio <= 4.0) {
				RX += ancho;
				RY += alto;
				nregs++;
			}
		}

		//Average area
		mAr /= vmedx.size();
		lAr = (int)(sqrt(mAr) + 0.5);
		lAr *= 0.9;

		if (nregs > 0) {
			RX /= nregs;
			RY /= nregs;
		}
		else {
			for (int i = 0; i<vSegUnits.size(); i++) {
				int ancho = vSegUnits[i].ROI.width + 1;
				int alto = vSegUnits[i].ROI.height + 1;

				RX += ancho;
				RY += alto;
				nregs++;
			}
			RX /= nregs;
			RY /= nregs;
		}

		//Compute median
		sort(vmedx.begin(), vmedx.end());
		sort(vmedy.begin(), vmedy.end());

		//Reference is the average of (mean,median,avg_area)
		RX = (RX + vmedx[vmedx.size() / 2] + lAr) / 3.0;
		RY = (RY + vmedy[vmedy.size() / 2] + lAr) / 3.0;
	}

	void computeSegDistance(int rx, int ry)
	{
		int N = vSegUnits.size();
		vvSegDist.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			vvSegDist[i].resize(N, 0.0);
		}

		float aux_x = rx;
		float aux_y = ry;
		NORMF = sqrt(aux_x*aux_x + aux_y*aux_y);
		INF_DIST = FLT_MAX / NORMF;

		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = i + 1; j < N; j++)
			{
				vvSegDist[i][j] = segDistance(i, j) / NORMF;
				vvSegDist[j][i] = vvSegDist[i][j];
			}
		}
	}

	int getSegUnitSize()
	{
		return vSegUnits.size();
	}

	void getSegUnitInfo(int segIdx, int NBest, std::vector<int> &vSymIdx,
						std::vector<float> &vProb, int &cmy, int &asc, int &des)
	{
		if (NBest != 1)
			HL_CERR("Request teh NB = 1 currently!!!");

		vSymIdx.resize(NBest);
		vProb.resize(NBest);

		SegUnit &curSeg = vSegUnits[segIdx];

		//Get from the Ground truth, just only one
		vSymIdx[0] = curSeg.symID;
		vProb[0] = 1.0;

		cv::Mat segMat = Img(curSeg.ROI);

		int n = 0;
		cmy = 0;
		float asc_ = 0, des_ = 0;
		float wasc = 0.1, wdes = 1.9;
		float paso = 1.8 / curSeg.ROI.height;
		float sumasc = 0, sumdes = 0;

		for (size_t i = 0; i < curSeg.ROI.height; i++)
		{
			uchar *pSegRow = segMat.ptr(i);
			for (size_t j = 0; j < curSeg.ROI.width; j++)
			{
				if (pSegRow[j] > 100)
				{
					int y = i + curSeg.ROI.y;
					sumasc += wasc;
					asc_ += y*wasc;

					n++;
					cmy += y;

					sumdes += wdes;
					des_ += y*wdes;
				}
			}

			wasc += paso;
			wdes -= paso;
		}

		asc = asc_ / sumasc;
		cmy = cmy / n;
		des = des_ / sumdes;
	}

	SegUnit &getSegUnit(int segIdx)
	{
		return vSegUnits[segIdx];
	}



	cv::Mat getRGBImg()
	{
		cv::Mat result;
		cv::cvtColor(Img, result, cv::COLOR_GRAY2BGR);
		return result;
	}

	float segDistanceOld(int si, int sj) {

		cv::Rect &segIROI = vSegUnits[si].ROI, &segJROI = vSegUnits[sj].ROI;

		cv::Point tl, br;
		cv::Point br1 = segIROI.br(), br2 = segJROI.br();
		tl.x = segIROI.x > segJROI.x ? segIROI.x : segJROI.x;
		tl.y = segIROI.y > segJROI.y ? segIROI.y : segJROI.y;
		br.x = br1.x < br2.x ? br1.x : br2.x;
		br.y = br1.y < br2.y ? br1.y : br2.y;

		float dist = FLT_MAX;

		if (tl.x <= br.x && tl.y <= br.y)
		{
			//The overlap is valid
			dist = -std::sqrt((br.x - tl.x)*(br.y - tl.y));
		}
		else
		{
			//When the overlap is invalid
			float dist1 = tl.x - br.x;
			float dist2 = tl.y - br.y;
			bool main_x = dist1 > dist2;
			dist = std::max(dist1, dist2);

			//Get Union ROI
			cv::Point utl, ubr;
			cv::Rect uROI = GetUnionRoi(segIROI, segJROI);
			utl = uROI.tl();
			ubr = uROI.br();

			cv::Rect emptyROI;
			emptyROI.x = tl.x < br.x ? tl.x : br.x;
			emptyROI.y = tl.y < br.y ? tl.y : br.y;
			emptyROI.width = std::abs(tl.x - br.x);
			emptyROI.height = std::abs(tl.y - br.y);

			//Check the visible if the overlap is invalid
			for (size_t i = 0; i < vSegUnits.size(); i++)
			{
				if (i == si || i == sj)continue;
				cv::Point stl = vSegUnits[i].ROI.tl(), sbr = vSegUnits[i].ROI.br();

				/*if ((stl.x >= utl.x && stl.y >= utl.y && stl.x < ubr.x && stl.y < ubr.y) ||
				(sbr.x >= utl.x && sbr.y >= utl.y && sbr.x < ubr.x && sbr.y < ubr.y))*/
				{
					cv::Rect olROI;
					if (GetOverlapRoi(emptyROI, vSegUnits[i].ROI, olROI))
					{
						float r = main_x ? (olROI.height / (emptyROI.height + 0.00001)) : (olROI.width / (emptyROI.width + 0.00001));
						if (r > 0.25)
						{
							//invisible dist
							dist = FLT_MAX;
							break;
						}
					}
				}
			}

		}

		return dist;
	}

	float segDistance(int si, int sj) {

		cv::Rect &segIROI = vSegUnits[si].ROI, &segJROI = vSegUnits[sj].ROI;

		cv::Point tl, br;
		cv::Point bri = segIROI.br(), brj = segJROI.br();
		tl.x = segIROI.x > segJROI.x ? segIROI.x : segJROI.x;
		tl.y = segIROI.y > segJROI.y ? segIROI.y : segJROI.y;
		br.x = bri.x < brj.x ? bri.x : brj.x;
		br.y = bri.y < brj.y ? bri.y : brj.y;

		float dist = FLT_MAX;

		if (tl.x <= br.x && tl.y <= br.y)
		{
			//The overlap is valid
			//dist = -std::sqrt((br.x - tl.x)*(br.y - tl.y));
			dist = 0;
		}
		else
		{
			//When the overlap is invalid, get the min distance of two box
			cv::Point cenPti((bri.x + segIROI.x)*0.5, (bri.y + segIROI.y)*0.5);
			cv::Point cenPtj((brj.x + segJROI.x)*0.5, (brj.y + segJROI.y)*0.5);
			cv::Point vcenIJ = cenPti - cenPtj;
			float cw = std::abs(vcenIJ.x), ch = std::abs(vcenIJ.y);
			float cenDist = std::sqrt(cw*cw + ch*ch);

			float k = ch / cw;
			float cutRatioI = (k * segIROI.width) < segIROI.height ? (segIROI.width * 0.5) / cw : (segIROI.height * 0.5) / ch;
			float cutRatioJ = (k * segJROI.width) < segJROI.height ? (segJROI.width * 0.5) / cw : (segJROI.height * 0.5) / ch;

			dist = (1.0 - cutRatioI - cutRatioJ) * cenDist;

			//Check the visible if the overlap is invalid

			//Get Union ROI
			cv::Point utl, ubr;
			cv::Rect uROI = GetUnionRoi(segIROI, segJROI);
			utl = uROI.tl();
			ubr = uROI.br();

			cv::Rect emptyROI;
			emptyROI.x = tl.x < br.x ? tl.x : br.x;
			emptyROI.y = tl.y < br.y ? tl.y : br.y;
			emptyROI.width = std::abs(tl.x - br.x);
			emptyROI.height = std::abs(tl.y - br.y);

			bool checkXSpace = tl.x > br.x, checkYSpace = tl.y > br.y;
			for (size_t i = 0; i < vSegUnits.size(); i++)
			{
				if (i == si || i == sj)continue;
				cv::Point stl = vSegUnits[i].ROI.tl(), sbr = vSegUnits[i].ROI.br();

				bool isInXSpace = checkXSpace &&
					((stl.x >= emptyROI.x && stl.x <= emptyROI.br().x) /*||
																	   (sbr.x >= emptyROI.x && sbr.x <= emptyROI.br().x)*/
					 );
				bool isInYSpace = checkYSpace &&
					((stl.y >= emptyROI.y && stl.y <= emptyROI.br().y) /*||
																	   (sbr.y >= emptyROI.y && sbr.y <= emptyROI.br().y)*/
					 );

				/*if ((stl.x >= utl.x && stl.y >= utl.y && stl.x < ubr.x && stl.y < ubr.y) &&
				(sbr.x >= utl.x && sbr.y >= utl.y && sbr.x < ubr.x && sbr.y < ubr.y))*/
				if (isInXSpace || isInYSpace)
				{
					cv::Rect olROI;
					if (GetOverlapRoi(emptyROI, vSegUnits[i].ROI, olROI))
					{
						float olArea = computeRectArea(olROI);
						float emptyArea = computeRectArea(emptyROI);
						float r = std::sqrt(olArea / (emptyArea + 0.000001));
						if (r > 0.05)
						{
							//invisible dist
							dist = FLT_MAX;
							break;
						}
					}
				}

			}

		}
		return dist;
	}

	float getDist(int si, int sj) {
		if (si<0 || sj<0 || si >= vvSegDist.size() || sj >= vvSegDist[si].size()) 
			HL_CERR("ERROR: segment id out of range in getDist(" << si << ", " << sj << ")");
		return vvSegDist[si][sj];
	}

	float group_penalty(std::vector<bool> &cccA, std::vector<bool> &cccB)
	{
		//Minimum or single-linkage clustering
		float dmin = FLT_MAX;
		for (int i = 0; i < cccA.size(); i++)
		{
			if (cccA[i]) 
			{
				for (int j = 0; j < cccB.size(); j++)
				{
					if (cccB[j] && j != i && getDist(i, j) < dmin)
					{
						dmin = getDist(i, j);
					}
				}
			}
		}

		return dmin;
	}

	//Normalized reference symbol size
	int RX, RY;

	//Move to CellCYK class
	/*void Sample::setRegion(CellCYK *c, int segIdx) {
	c->ccc[segIdx] = true;

	c->x = vSegUnits[segIdx].ROI.x;
	c->y = vSegUnits[segIdx].ROI.y;
	c->s = vSegUnits[segIdx].ROI.br().x;
	c->t = vSegUnits[segIdx].ROI.br().y;
	}*/

public:
	float INF_DIST;  //Infinite distance value (visibility)
	float NORMF;     //Normalization factor for distances

private:
	cv::Mat Img;
	std::vector<SegUnit> vSegUnits;
	std::vector<std::vector<float>> vvSegDist;
	int W, H;

	std::shared_ptr<SymSet> pSymSet;
};
