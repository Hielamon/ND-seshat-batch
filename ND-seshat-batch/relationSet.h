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
	std::vector<cv::Point> vPoint1, vPoint2;
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

inline std::string stringRelType(RelationType &relType)
{
	std::string result;
	switch (relType)
	{
	case H:
		result = "Horizontal(0)";
		break;
	case SUB:
		result = "Subscript(1)";
		break;
	case SUP:
		result = "Superscript(2)";
		break;
	case V:
		result = "Vertical(3)";
		break;
	case INS:
		result = "Inside(4)";
		break;
	case MROOT:
		result = "mth root(5)";
		break;
	default:
		break;
	}

	return result;
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
		std::ofstream fs(fName, std::ios::out);
		if (!fs.is_open())
			HL_CERR_RETURN_FALSE("Failed to open file " + fName);

		fs << totalSize() << std::endl;

		saveRelationUnits(fs, vRelH);
		saveRelationUnits(fs, vRelSub);
		saveRelationUnits(fs, vRelSup);
		saveRelationUnits(fs, vRelV);
		saveRelationUnits(fs, vRelIns);
		saveRelationUnits(fs, vRelMroot);

		fs.close();
	}

	bool ReadFromFile(const std::string &fName, std::vector<RelationUnit> &vrelUnit)
	{
		std::ifstream fs(fName, std::ios::in);
		if (!fs.is_open())
			HL_CERR_RETURN_FALSE("Failed to open file " + fName);

		int totalNum;
		fs >> totalNum;
		char typeC;
		
		for (size_t i = 0; i < totalNum; i++)
		{
			RelationUnit relUnit;
			loadRelationUnit(fs, typeC, relUnit);
			vrelUnit.push_back(relUnit);
			addRelationUnit(relUnit);
		}

		fs.close();

		return true;
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

	int totalSize()
	{
		return vRelH.size() + vRelSub.size() + vRelSup.size() + vRelV.size() + vRelIns.size() + vRelMroot.size();
	}

	void drawRelationUnits(cv::Mat& img, std::vector<RelationUnit> &vRelUnits,
						   cv::Scalar color, std::string tag, int padding, double scale)
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

	void saveRelationUnits(std::ofstream &fs, std::vector<RelationUnit> &vRelUnits)
	{
		for (size_t i = 0; i < vRelUnits.size(); i++)
		{
			RelationUnit &relUnit = vRelUnits[i];
			fs << getTypeC(relUnit.relType) << std::endl;
			fs << relUnit.ROI1.x << " " << relUnit.ROI1.y << " "
				<< relUnit.ROI1.width << " " << relUnit.ROI1.height << " "
				<< relUnit.ROI2.x << " " << relUnit.ROI2.y << " "
				<< relUnit.ROI2.width << " " << relUnit.ROI2.height << std::endl;

			for (size_t j = 0; j < relUnit.vPoint1.size(); j++)
			{
				cv::Point &pt = relUnit.vPoint1[j];
				fs << pt.x << " " << pt.y << " ";
			}
			fs << std::endl;

			for (size_t j = 0; j < relUnit.vPoint2.size(); j++)
			{
				cv::Point &pt = relUnit.vPoint2[j];
				fs << pt.x << " " << pt.y << " ";
			}
			fs << std::endl;
		}
	}

	void loadRelationUnit(std::ifstream &fs, char &typeC, RelationUnit &relUnit)
	{
		fs >> typeC;

		relUnit.relType = getRelationType(typeC);
		fs >> relUnit.ROI1.x >> relUnit.ROI1.y
			>> relUnit.ROI1.width >> relUnit.ROI1.height
			>> relUnit.ROI2.x >> relUnit.ROI2.y
			>> relUnit.ROI2.width >> relUnit.ROI2.height;
		std::string line;
		std::stringstream ioStr;
		std::getline(fs, line);

		//The string of points1
		std::getline(fs, line);
		relUnit.vPoint1.clear();
		relUnit.vPoint2.clear();

		ioStr << line;
		int x, y;
		while (ioStr >> x && ioStr >> y)
		{
			relUnit.vPoint1.push_back(cv::Point(x, y));
		}

		//The string of points2
		std::getline(fs, line);
		ioStr.clear();
		ioStr << line;
		while (ioStr >> x && ioStr >> y)
		{
			relUnit.vPoint2.push_back(cv::Point(x, y));
		}
	}

	char getTypeC(RelationType &relType)
	{
		char result;
		switch (relType)
		{
		case H:
			result = 'H';
			break;
		case SUB:
			result = 'B';
			break;
		case SUP:
			result = 'P';
			break;
		case V:
			result = 'V';
			break;
		case INS:
			result = 'I';
			break;
		case MROOT:
			result = 'M';
			break;
		default:
			break;
		}

		return result;
	}

	RelationType getRelationType(char &typeC)
	{
		RelationType type;
		switch (typeC)
		{
		case 'H':
			type = RelationType::H;
			break;
		case 'B':
			type = RelationType::SUB;
			break;
		case 'P':
			type = RelationType::SUP;
			break;
		case 'V':
			type = RelationType::V;
			break;
		case 'I':
			type = RelationType::INS;
			break;
		case 'M':
			type = RelationType::MROOT;
			break;
		default:
			HL_CERR("Wrong type char input");
			break;
		}

		return type;
	}

	void addRelationUnit(RelationUnit &relUnit)
	{
		switch (relUnit.relType)
		{
		case H:
			vRelH.push_back(relUnit);
			break;
		case SUB:
			vRelSub.push_back(relUnit);
			break;
		case SUP:
			vRelSup.push_back(relUnit);
			break;
		case V:
			vRelV.push_back(relUnit);
			break;
		case INS:
			vRelIns.push_back(relUnit);
			break;
		case MROOT:
			vRelMroot.push_back(relUnit);
			break;
		default:
			break;
		}
	}
};
