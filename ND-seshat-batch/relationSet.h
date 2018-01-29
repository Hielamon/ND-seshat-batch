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

private:
	std::vector<RelationUnit> vRelH, vRelSub, vRelSup, vRelV, vRelIns, vRelMroot;
};
