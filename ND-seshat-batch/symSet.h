#pragma once

#include <commonMacro.h>

#include <map>
#include <vector>
#include <fstream>
#include <string>

#define STD_CENTER_Y 3.1
#define STD_TOP_Y 0.0
#define STD_BOTTOM_Y 5.0

struct StdSymInfo
{
	std::string symStr;
	double rel_y, rel_t;
	double rel_w, rel_h;
};

class SymSet
{
	std::map<std::string, int> cl2key;
	std::vector<std::string> key2cl;

	//the information of standard symbol
	std::vector<StdSymInfo> vSymStdInfo;

	//std::vector<int> type;
public:
	SymSet(){}
	~SymSet(){}

	bool load(const std::string &filename)
	{
		if (!key2cl.empty()) key2cl.clear();
		if (!vSymStdInfo.empty()) vSymStdInfo.clear();
		//if (!type.empty()) type.clear();
		if (!cl2key.empty()) cl2key.clear();

		std::fstream fs(filename, std::ios::in);
		if (!fs.is_open())
			HL_CERR("Failed to open the file " + filename);

		int claseNum;
		fs >> claseNum;

		assert(claseNum >= 0);

		//char T;
		std::string clase;
		double rel_y, rel_t, rel_w, rel_h;

		for (int i = 0; i < claseNum; i++)
		{
			fs >> rel_y >> rel_t >> rel_w >> rel_h >> clase;
			key2cl.push_back(clase);
			cl2key[clase] = i;

			StdSymInfo sym;
			sym.rel_y = rel_y;
			sym.rel_t = rel_t;
			sym.rel_w = rel_w;
			sym.rel_h = rel_h;
			sym.symStr = clase;
			vSymStdInfo.push_back(sym);
		}

		fs.close();
		return true;
	}

	std::string& strClase(int cIdx) 
	{
		return key2cl[cIdx];
	}

	int keyClase(const std::string &claseName) 
	{
		if (cl2key.find(claseName) == cl2key.end()) {
			HL_CERR("WARNING: Class " << claseName  <<" doesn't appear in symbols database ");
		}
		return cl2key[claseName];
	}

	bool checkClase(const std::string &claseName)
	{
		if (cl2key.find(claseName) == cl2key.end())
			return false;
		return true;
	}

	int getNClases() 
	{
		return key2cl.size();
	}

	StdSymInfo &stdInfoClase(int cIdx)
	{
		return vSymStdInfo[cIdx];
	}

};