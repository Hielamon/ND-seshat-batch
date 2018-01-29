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

inline cv::Point ROICenter(cv::Rect &ROI, double scale)
{
	return cv::Point(ROI.x + ROI.width * 0.5, ROI.y + ROI.height * 0.5) * scale;
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

	int ShowInImage(cv::Mat &Img, const std::string &winName)
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
		padding = 1;
		drawRelationUnits(showImg, vRelSub, color, tag, padding, scale);

		color = cv::Scalar(219, 112, 147);
		tag = "SUP";
		padding = -1;
		drawRelationUnits(showImg, vRelSup, color, tag, padding, scale);

		color = cv::Scalar(255, 255, 0);
		tag = "V";
		padding = -2;
		drawRelationUnits(showImg, vRelV, color, tag, padding, scale);

		color = cv::Scalar(255, 48, 155);  
		tag = "INS";
		padding = 2;
		drawRelationUnits(showImg, vRelIns, color, tag, padding, scale);

		color = cv::Scalar(43, 90, 139);
		tag = "MROOT";
		padding = -3;
		drawRelationUnits(showImg, vRelMroot, color, tag, padding, scale);

		//cv::imshow("Show the Relation Sets", showImg);
		cv::imshow(winName, showImg);
		return cv::waitKey(0);
	}

private:
	std::vector<RelationUnit> vRelH, vRelSub, vRelSup, vRelV, vRelIns, vRelMroot;

	void drawRelationUnits(cv::Mat img, std::vector<RelationUnit> vRelUnits, cv::Scalar color, std::string tag, int padding, double scale)
	{
		for (size_t i = 0; i < vRelUnits.size(); i++)
		{
			RelationUnit &ru = vRelUnits[i];
			cv::rectangle(img, paddingROI(ru.ROI1, scale, padding), color, 1, cv::LINE_AA);
			cv::rectangle(img, paddingROI(ru.ROI2, scale, padding), color, 1, cv::LINE_AA);
			cv::Point ptS = ROICenter(ru.ROI1, scale), ptE = ROICenter(ru.ROI2, scale);
			cv::arrowedLine(img, ptS, ptE, color, 1, cv::LINE_AA);
			cv::putText(img, tag, (ptS + ptE)*0.5, cv::HersheyFonts::FONT_HERSHEY_COMPLEX, 0.8, color);
		}
	}
};
