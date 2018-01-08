#pragma once

#include <commonMacro.h>

#include <map>
#include <vector>
#include <fstream>
#include <string>

class SymSet
{
	std::map<std::string, int> cl2key;
	std::vector<std::string> key2cl;
	std::vector<int> type;
public:
	SymSet(){}
	~SymSet(){}

	bool load(const std::string &filename)
	{
		if(!key2cl.empty()) key2cl.clear();
		if (!type.empty()) type.clear();
		if (!cl2key.empty()) cl2key.clear();

		std::fstream fs(filename, std::ios::in);
		if (!fs.is_open())
			HL_CERR("Failed to open the file " + filename);

		int claseNum;
		fs >> claseNum;

		assert(claseNum >= 0);

		char T;
		std::string clase;

		for (int i = 0; i < claseNum; i++)
		{
			fs >> clase >> T;
			key2cl.push_back(clase);
			cl2key[clase] = i;

			if (T == 'n')		type.push_back(0); //Centroid
			else if (T == 'a')  type.push_back(1); //Ascender
			else if (T == 'd')  type.push_back(2); //Descender
			else if (T == 'm')  type.push_back(3); //Middle
			else
			{
				fprintf(stderr, "SymSet: Error reading symbol types\n");
				exit(-1);
			}
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

	int symType(int k) 
	{
		return type[k];
	}
};