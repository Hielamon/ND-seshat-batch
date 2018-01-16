#define MAIN_FILE
#include <commonMacro.h>
#include <TraverFolder.h>
#include <iostream>
#include <direct.h>
#include <io.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "grammar.h"
#include "symSet.h"
#include "sample.h"
#include "meParser.h"

std::string batchFileList = "D:/Funny-Works/Academic-Codes/HandWritten/Datasets/TidyDatasets/TidyDatasets/UniformTestSet/filename.txt";
std::string VOC2007CharMapFName = "D:/Funny-Works/Academic-Codes/HandWritten/Datasets/VOC2007/charmap_.txt";
bool IsShowSample = !false;
bool saveResult = !true;
bool withGT = true;
bool symErrStop = false;

inline bool getSymbolMap(const std::string &filename, std::map<std::string, int> &symbolMap)
{
	if (!symbolMap.empty()) symbolMap.clear();
	std::fstream fs(filename, std::ios::in);
	if (!fs.is_open())
		HL_CERR_RETURN_FALSE("Failed to open the file " + filename);
	int charMapSize;
	fs >> charMapSize;
	for (size_t i = 0; i < charMapSize; i++)
	{
		std::string symbolName;
		int mapID;
		fs >> symbolName >> mapID;
		if (symbolMap.find(symbolName) == symbolMap.end())
		{
			symbolMap[symbolName] = mapID;
		}
		else
			HL_CERR_RETURN_FALSE("Symbol is duplicated in charmap file " + filename);
	}
	fs.close();
}

bool checkLatex(std::string &latex1, std::string &latex2)
{
	int i = 0, j = 0;
	while (i < latex1.size() && j < latex2.size())
	{
		if (latex1[i] != ' ' && latex2[j] != ' ')
		{
			if (latex1[i] != latex2[j]) break;
			i++;
			j++;
		}
		else
		{
			if (latex1[i] == ' ') i++;
			if (latex2[j] == ' ') j++;
		}
	}

	while (i < latex1.size() && latex1[i] == ' ')i++;
	while (j < latex2.size() && latex2[j] == ' ')j++;

	return i == latex1.size() && j == latex2.size();
}

bool checkLatexNew(std::string &latex1, std::string &latex2)
{
	std::string s1, s2;
	for (size_t i = 0; i < latex1.size(); i++)
		if (latex1[i] != ' ')s1.push_back(latex1[i]);

	for (size_t i = 0; i < latex2.size(); i++)
		if (latex2[i] != ' ')s2.push_back(latex2[i]);

	std::vector<std::vector<std::string>> replaceMap = {
		{ "{\\int}", "int" },
		{ "{\\log}", "log" },
		{ "\\int", "int" },
		{ "\\log", "log" },
		{ "\\\\", "\\" },
		{ "\\{", "{" },
		{ "\\}", "}" },
		//{ "\\lg", "lg" }
	};

	size_t pos = std::string::npos;
	for (size_t i = 0; i < replaceMap.size(); i++)
	{
		std::string &replaceStr = replaceMap[i][0], &aimStr = replaceMap[i][1];
		while ((pos = s1.find(replaceStr)) != std::string::npos)
			s1.replace(pos, replaceStr.size(), aimStr);

		while ((pos = s2.find(replaceStr)) != std::string::npos)
			s2.replace(pos, replaceStr.size(), aimStr);
	}

	std::cout << "s1     = " << s1 << std::endl;
	std::cout << "s2(GT) = " << s2 << std::endl;
	return s1 == s2;
}

void printUsage()
{

}

int parseCmdArgs(int argc, char** argv)
{

	if (argc == 1)
	{
		printUsage();
		return 0;
	}

	for (int i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-help" || std::string(argv[i]) == "/?")
		{
			printUsage();
			return 0;
		}
		else if (std::string(argv[i]) == "-batch_file")
		{
			batchFileList = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-withGT")
		{
			if (std::string(argv[i + 1]) == "yes")
				withGT = true;
			else if (std::string(argv[i + 1]) == "no")
				withGT = false;
			i++;
		}
		else if (std::string(argv[i]) == "-showSample")
		{
			if (std::string(argv[i + 1]) == "yes")
				IsShowSample = true;
			else if (std::string(argv[i + 1]) == "no")
				IsShowSample = false;
			i++;
		}
		else if (std::string(argv[i]) == "-save")
		{
			if (std::string(argv[i + 1]) == "yes")
				saveResult = true;
			else if (std::string(argv[i + 1]) == "no")
				saveResult = false;
			i++;
		}
		else if (std::string(argv[i]) == "-symErrStop")
		{
			if (std::string(argv[i + 1]) == "yes")
				symErrStop = true;
			else if (std::string(argv[i + 1]) == "no")
				symErrStop = false;
			i++;
		}
		else
		{
			return 0;
		}
	}

	return 1;
}

int main(int argc, char *argv[])
{
	parseCmdArgs(argc, argv);

	std::shared_ptr<SymSet> pSymSet = std::make_shared<SymSet>();
	pSymSet->load("Grammar/symbol_nd.types");

	std::shared_ptr<Grammar> pG = std::make_shared<Grammar>();
	pG->reSetup("Grammar/mathexp.gram", pSymSet);

	char gmmfile[1024] = "Grammar/sparels.gmm";
	std::shared_ptr<GMM> pGMM = std::make_shared<GMM>(gmmfile);

	MeParser meparser(pSymSet, pG, pGMM);

	std::shared_ptr<Sample> sample = std::make_shared<Sample>(pSymSet);

	std::map<std::string, int> symbolMap;
	getSymbolMap(VOC2007CharMapFName, symbolMap);

	std::fstream fs(batchFileList, std::ios::in);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << batchFileList);

	std::string sampleFileName, line, gtLatex, imgPath;
	std::vector<std::string> vSymErrorNames;
	std::vector<std::string> vParseErrorNames, vParseErrorLatexs, vParseErrorGTLatexs;
	std::vector<std::string> vTrueNames, vTrueLatexs, vTrueGTLatexs;

	std::string sepLine(120, '-');

	std::stringstream ioStr;
	while (std::getline(fs, line) && !line.empty())
	{
		std::cout << "\n\n" << sepLine << std::endl;
		ioStr.clear();
		ioStr.str("");
		ioStr << line;
		ioStr >> sampleFileName;
		std::cout << "Excute file : " << sampleFileName << std::endl;
		if (sample->LoadFromUnifromFile(sampleFileName, gtLatex, imgPath, symbolMap, withGT))
		{
			std::string latexResult = meparser.parse(sample);

			if (gtLatex.empty())
			{
				vTrueNames.push_back(sampleFileName);
				vTrueLatexs.push_back(latexResult);
				vTrueGTLatexs.push_back("");
				//std::cout << "EMPTY" << std::endl;
			}
			else
			{
				if (checkLatexNew(latexResult, gtLatex))
				{
					vTrueNames.push_back(sampleFileName);
					vTrueLatexs.push_back(latexResult);
					vTrueGTLatexs.push_back(gtLatex);
					std::cout << "True Sample" << std::endl;
				}
				else
				{
					//Parse to wrong Result
					vParseErrorNames.push_back(sampleFileName);
					vParseErrorLatexs.push_back(latexResult);
					vParseErrorGTLatexs.push_back(gtLatex);
					std::cout << "False Sample" << std::endl;
				}
				std::cout << "GT Latex : " << gtLatex << std::endl;
			}
			std::cout << "Image Path : " << imgPath << std::endl;
			size_t pos = sampleFileName.rfind("\\");
			std::string pureFileName;
			if (pos == std::string::npos)
				pureFileName.assign(sampleFileName.begin(), sampleFileName.end() - 4);
			else
				pureFileName.assign(sampleFileName.begin() + pos, sampleFileName.end() - 4);
			if (IsShowSample)
			{
				sample->ShowSample(pureFileName);
				cv::destroyWindow(pureFileName);
			}


		}
		else
		{
			//invalid samples , need to fix the symbol set or grammar rules
			vSymErrorNames.push_back(sampleFileName);
			if (symErrStop)
				system("PAUSE");
		}
	}

	fs.close();
	int nsErr = vSymErrorNames.size(), npErr = vParseErrorNames.size(), ntrue = vTrueNames.size();
	std::cout << "Total samples : " << nsErr + npErr + ntrue << std::endl;

	std::cout << "Symbol Error samples : " << nsErr << std::endl;
	std::cout << "Accurancy : " << ntrue / float(npErr + ntrue) << "( " << ntrue << "/" << npErr + ntrue << " )" << std::endl;

	if (!saveResult)return 0;
	std::cout << "\n\n" << sepLine << std::endl;
	std::cout << "Start to save the result" << std::endl;

	std::string resultDir = "Result";
	//Check and Create the folder
	if (_access(resultDir.c_str(), 0) == -1)
		_mkdir(resultDir.c_str());

	std::string symErrorFName = resultDir + "/" + "symErrorSampleList.txt";
	std::string parseErrorFName = resultDir + "/" + "parseErrorSampleList.txt";
	std::string trueFName = resultDir + "/" + "trueSampleList.txt";

	//Save the symbol set error file list
	fs.open(symErrorFName, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << symErrorFName);

	for (size_t i = 0; i < vSymErrorNames.size(); i++)
		fs << vSymErrorNames[i] << std::endl;

	fs.close();

	//Save the latex ground truth set error file list
	fs.open(parseErrorFName, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << parseErrorFName);

	for (size_t i = 0; i < vParseErrorNames.size(); i++)
		fs << vParseErrorNames[i] << " ------>" << vParseErrorLatexs[i] << "------GT>" << vParseErrorGTLatexs[i] << std::endl;

	fs.close();

	//Save the true sample file list
	fs.open(trueFName, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << trueFName);

	for (size_t i = 0; i < vTrueNames.size(); i++)
		fs << vTrueNames[i] << " ------>" << vTrueLatexs[i] << "------GT>" << vTrueGTLatexs[i] << std::endl;

	fs.close();

	return 0;
}