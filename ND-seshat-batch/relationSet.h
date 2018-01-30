#pragma once
#include <OpencvCommon.h>
#include <commonMacro.h>
#include <iostream>

enum RelationType
{
	H, SUB, SUP, V, INS, MROOT
};

struct RelationUnit
{
	cv::Rect ROI1, ROI2;
	RelationType relType;
};

inline cv::Rect paddingROI(cv::Rect &ROI, double scale, int padding)
{
	return cv::Rect(ROI.x * scale - padding, ROI.y * scale - padding, ROI.width * scale + 2 * padding, ROI.height * scale + 2 * padding);
}

inline cv::Point ROICenter(cv::Rect &ROI)
{
	return cv::Point(ROI.x + ROI.width * 0.5, ROI.y + ROI.height * 0.5);
}

class RelationSet
{
public:
	RelationSet() {}
	~RelationSet() {}

	void addRelation(RelationUnit &ru)
	{
		switch (ru.relType)
		{
		case H:
			vRelH.push_back(ru);
			break;
		case SUB:
			vRelSub.push_back(ru);
			break;
		case SUP:
			vRelSup.push_back(ru);
			break;
		case V:
			vRelV.push_back(ru);
			break;
		case INS:
			vRelIns.push_back(ru);
			break;
		case MROOT:
			vRelMroot.push_back(ru);
			break;
		default:
			HL_CERR("Invalid space relation type : " << ru.relType);
			break;
		}
	}

	bool SaveToFile(const std::string &fName)
	{

	}

	bool ReadFromFile(const std::string &fName)
	{

	}

	void DrawInImage(const cv::Mat &Img, cv::Mat &out)
	{
		int W = Img.cols, H = Img.rows;
		double megapix = std::max(0.5 * 1e6, double(W*H));
		double scale = std::sqrt(megapix / (W * H));
		cv::Mat showImg = Img.clone();
		cv::resize(showImg, showImg, cv::Size(W*scale, H*scale));

		cv::Scalar color;
		std::string tag;
		int padding;
		color = cv::Scalar(0, 165, 255);
		tag = "H";
		padding = 0;
		drawRelationUnits(showImg, vRelH, color, tag, padding, scale);

		color = cv::Scalar(180, 105, 255);
		tag = "SUB";
		padding = 2;
		drawRelationUnits(showImg, vRelSub, color, tag, padding, scale);

		color = cv::Scalar(0, 0, 255);
		tag = "SUP";
		padding = -2;
		drawRelationUnits(showImg, vRelSup, color, tag, padding, scale);

		color = cv::Scalar(255, 255, 0);
		tag = "V";
		padding = -4;
		drawRelationUnits(showImg, vRelV, color, tag, padding, scale);

		color = cv::Scalar(255, 48, 155);  
		tag = "INS";
		padding = 4;
		drawRelationUnits(showImg, vRelIns, color, tag, padding, scale);

		color = cv::Scalar(43, 90, 139);
		tag = "MROOT";
		padding = 6;
		drawRelationUnits(showImg, vRelMroot, color, tag, padding, scale);

		out = showImg;
	}

private:
	std::vector<RelationUnit> vRelH, vRelSub, vRelSup, vRelV, vRelIns, vRelMroot;

	void drawRelationUnits(cv::Mat img, std::vector<RelationUnit> vRelUnits, cv::Scalar color, std::string tag, int padding, double scale)
	{
		for (size_t i = 0; i < vRelUnits.size(); i++)
		{
			RelationUnit &ru = vRelUnits[i];
			cv::Rect pROI1 = paddingROI(ru.ROI1, scale, padding);
			cv::Rect pROI2 = paddingROI(ru.ROI2, scale, padding);
			cv::rectangle(img, pROI1, color, 1, cv::LINE_AA);
			cv::rectangle(img, pROI2, color, 1, cv::LINE_AA);
			cv::Point ptS = ROICenter(pROI1), ptE = ROICenter(pROI2);
			if (tag == "MROOT" || tag == "INS")
			{
				ptS = pROI1.tl();
				ptE = pROI2.tl();
			}
			cv::arrowedLine(img, ptS, ptE, color, 1, cv::LINE_AA);
			cv::putText(img, tag, (ptS + ptE)*0.5, cv::HersheyFonts::FONT_HERSHEY_COMPLEX, 0.8, color);
		}
	}
};
