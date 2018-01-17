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
#ifndef _GRAMMAR_
#define _GRAMMAR_

#include <cstdio>
#include <string>
#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <memory>
#include "production.h"
#include "symSet.h"

inline bool isFillChar(char c)
{
	switch (c)
	{
	case ' ':
	case '\t':
	case '\n':
	case '\r':
		return true;
	default:
		return false;
	}
}

inline int split(const std::string &str, std::vector<std::string> &res)
{
	char tokensaux[2 * 1024];
	int n = 0, i = 0, j = 0;

	while (isFillChar(str[i]))  i++;

	while (str[i]) {
		if (str[i] == '\"') {
			i++;
			while (str[i] && str[i] != '\"') {
				tokensaux[j] = str[i];
				i++; j++;
			}
			i++;
		}
		else {
			while (str[i] && !isFillChar(str[i])) {
				tokensaux[j] = str[i];
				i++; j++;
			}
		}
		tokensaux[j++] = 0;
		n++;
		while (str[i] && isFillChar(str[i]))  i++;
	}

	if (!res.empty())res.clear();
	res.resize(n);
	for (i = 0, j = 0; i < n; i++) {
		int tlen = strlen(&tokensaux[j]) + 1;
		res[i] = std::string(tokensaux + j);
		j += tlen;
	}

	return n;
}

inline bool nextLine(std::ifstream &fs, std::string &line)
{
	do {
		if (!std::getline(fs, line))
		{
			//std::cout << "End of file" << std::endl;
			return false;
		}
		//std::cout << line << std::endl;
			
	} while (line.empty() || line[0] == '#');

	return true;
}

class Grammar
{
public:
	std::map<std::string, int> noTerminales;
	std::vector<int> initsyms;
	std::vector<bool> esInit;
	std::shared_ptr<SymSet> pSymSet;

	enum PBTYPE
	{
		H, SUP, SUB, V, VE, INS, MRT, SSE
	};

	std::string strType(PBTYPE pType)
	{
		std::string result = "None-defined";
		switch (pType)
		{
		case Grammar::H:
			result = "H";
			break;
		case Grammar::SUP:
			result = "SUP";
			break;
		case Grammar::SUB:
			result = "SUB";
			break;
		case Grammar::V:
			result = "V";
			break;
		case Grammar::VE:
			result = "VE";
			break;
		case Grammar::INS:
			result = "INS";
			break;
		case Grammar::MRT:
			result = "MRT";
			break;
		case Grammar::SSE:
			result = "SSE";
			break;
		default:
			break;
		}
		return result;
	}

	//TODO: use the glogal prods variable
	std::vector<std::shared_ptr<ProductionB>> prodsH, prodsSup, prodsSub;
	std::vector<std::shared_ptr<ProductionB>> prodsV, prodsVe, prodsIns, prodsMrt, prodsSSE;
	std::vector<std::shared_ptr<ProductionT>> prodTerms;
	Grammar(const std::string &gramFName, const std::shared_ptr<SymSet> &pSymSet_)
	{
		reSetup(gramFName, pSymSet_);
	}
	Grammar() {}
	~Grammar() {}

	const std::string &key2str(int k)
	{
		for (auto iter = noTerminales.begin(); iter != noTerminales.end(); iter++)
		{
			if (iter->second == k)
			{
				return iter->first;
			}
		}

		return "NULL";
	}

	bool isValid()
	{
		//TODO : more things need be checked
		return pSymSet.use_count();
	}

	bool reSetup(const std::string &gramFName, const std::shared_ptr<SymSet> &pSymSet_)
	{
		clear();

		pSymSet = pSymSet_;

		std::ifstream gramfs(gramFName, std::ios::in);
		
		if (!gramfs.is_open())
			HL_CERR("Error loading grammar " + gramFName);

		std::stringstream ioStr;
		std::string linea, tok1, tok2, aux;

		//Read nonterminal symbols
		while (nextLine(gramfs, linea) && linea != "START")
		{
			ioStr.clear();
			ioStr << linea;
			ioStr >> tok1;
			addNoTerminal(tok1);
		}

		//Read start symbol(s) of the grammar
		while (nextLine(gramfs, linea) && linea != "PTERM")
		{
			ioStr.clear();
			ioStr << linea;
			ioStr >> tok1;
			addInitSym(tok1);
		}

		//Read terminal productions
		while (nextLine(gramfs, linea) && linea != "PBIN")
		{
			float pr; 
			ioStr.clear();
			ioStr << linea;
			ioStr >> pr >> tok1 >> tok2 >> aux;
			addTerminal(pr, tok1, tok2, aux);
		}

		//Read binary productions
		while (nextLine(gramfs, linea))
		{
			std::vector<std::string> vTokens;
			int ntoks = split(linea, vTokens);

			if (ntoks != 7)
				HL_CERR("Error: Grammar not valid (PBIN)");

			addBinaryRule(atof(vTokens[0].c_str()), vTokens[1], vTokens[2], vTokens[3], vTokens[4], vTokens[5], vTokens[6]);
		}

		gramfs.close();

		esInit.resize(noTerminales.size(), false);

		for (size_t i = 0; i < initsyms.size(); i++)
			esInit[initsyms[i]] = true;

		return true;
	}

private:

	void clear()
	{
		if (!noTerminales.empty()) noTerminales.clear();
		if (!initsyms.empty()) initsyms.clear();
		if (!esInit.empty()) esInit.clear();
		if (!prodTerms.empty()) prodTerms.clear();
		if (!prodsH.empty()) prodsH.clear();
		if (!prodsSup.empty()) prodsSup.clear();
		if (!prodsSub.empty()) prodsSub.clear();
		if (!prodsV.empty()) prodsV.clear();
		if (!prodsVe.empty()) prodsVe.clear();
		if (!prodsIns.empty()) prodsIns.clear();
		if (!prodsMrt.empty()) prodsMrt.clear();
		if (!prodsSSE.empty()) prodsSSE.clear();
		if (pSymSet.use_count())pSymSet.reset();
	}

	void addInitSym(const std::string& initSymbol)
	{
		if (noTerminales.find(initSymbol) == noTerminales.end())
			HL_CERR("addInitSym: Non-terminal " << initSymbol << " not defined.");

		initsyms.push_back(noTerminales[initSymbol]);
	}

	void addNoTerminal(const std::string& ntSymbol)
	{
		int key = noTerminales.size();
		noTerminales[ntSymbol] = key;
	}

	void addTerminal(float pr, const std::string &S, const std::string &T, const std::string &tex)
	{
		if (noTerminales.find(S) == noTerminales.end())
			HL_CERR("addInitSym: Non-terminal " << S << " not defined.");

		bool create = true;
		int id = pSymSet->keyClase(T);
		if (id < 0)
			HL_CERR("ERROR: " << S << " -> " << T << " (id < 0)");

		for (size_t i = 0; i < prodTerms.size(); i++)
		{
			if (prodTerms[i]->getNoTerm() == noTerminales[S]) {
				prodTerms[i]->setClase(id, pr, tex, 'i');
				create = false;
				break;
			}
		}

		if (create) {
			std::shared_ptr<ProductionT> pt = std::make_shared<ProductionT>(noTerminales[S], pSymSet->getNClases(), S);
			pt->setClase(id, pr, tex, 'i');
			prodTerms.push_back(pt);
		}
	}

	bool addBinaryRule(float pr, const std::string &ruleType, const std::string &S, const std::string &A, const std::string &B,
					   const std::string &out, const std::string &merge)
	{
		if (noTerminales.find(S) == noTerminales.end())
			HL_CERR("Rule: Non-terminal " << S << " not defined.");
		if (noTerminales.find(A) == noTerminales.end())
			HL_CERR("Rule: Non-terminal " << A << " not defined.");
		if (noTerminales.find(B) == noTerminales.end())
			HL_CERR("Rule: Non-terminal " << B << " not defined.");

		std::shared_ptr<ProductionB> pd;

		if (ruleType == "H")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionH>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsH.push_back(pd);
		}
		else if (ruleType == "V")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionV>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsV.push_back(pd);
		}
		else if (ruleType == "Ve")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionVe>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsVe.push_back(pd);
		}
		else if (ruleType == "Sup")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionSup>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsSup.push_back(pd);
		}
		else if (ruleType == "Sub")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionSub>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsSub.push_back(pd);
		}
		else if (ruleType == "SSE")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionSSE>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsSSE.push_back(pd);
		}
		else if (ruleType == "Ins")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionIns>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsIns.push_back(pd);
		}
		else if (ruleType == "Mrt")
		{
			pd = std::static_pointer_cast<ProductionB>(std::make_shared<ProductionMrt>(
				noTerminales[S], noTerminales[A], noTerminales[B], S, A, B, pr, out));
			pd->setMerges(merge[0]);
			prodsMrt.push_back(pd);
		}
		else
			HL_CERR("Error: Binary rule type " << ruleType << " nor valid\n");


	}
};

#endif
