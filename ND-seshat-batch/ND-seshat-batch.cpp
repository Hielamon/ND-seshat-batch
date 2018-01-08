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

std::string testXML = "VOC2007/Annotations/hwf_0000007.xml";
std::string imgPath = "VOC2007/JPEGImages/";
std::string batchFileList = "testfilename.txt";
//std::string batchFileList = "VOC2007/Annotations/filename.txt";
std::string latexListMapFName = "VOC2007/latexListMap.txt";
bool IsShowSample = false;
bool saveResult = true;

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
		{"\\log", "log"},
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

bool getLatexListMap(const std::string &filename, std::map<std::string, std::string> &latexListMap)
{
	std::fstream fs(filename, std::ios::in);
	if (!fs.is_open())
		HL_CERR_RETURN_FALSE("Failed to open the file " << filename);

	std::string line;
	while (std::getline(fs, line) && !line.empty())
	{
		size_t spos = line.find(" ");
		std::string latexID, latex;
		latexID.assign(line.begin(), line.begin() + spos);
		latex.assign(line.begin() + spos + 1, line.end());
		latexListMap[latexID] = latex;
	}

	return true;
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
		else if (std::string(argv[i]) == "-xml")
		{
			testXML = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-imgdir")
		{
			imgPath = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-batch_file")
		{
			batchFileList = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-latexGT_file")
		{
			latexListMapFName = std::string(argv[i + 1]);
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

	std::map<std::string, std::string> latexListMap;
	getLatexListMap(latexListMapFName, latexListMap);

	std::fstream fs(batchFileList, std::ios::in);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << batchFileList);

	std::string sampleFileName, line, gtLatexID;
	std::map<std::string, std::string>::iterator it = latexListMap.end();
	std::vector<std::string> vSymErrorNames, vLatexErrorNames;
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
		if (sample->LoadFromVOC2007XML(sampleFileName, imgPath, gtLatexID))
		{
			std::string latexResult = meparser.parse(sample);

			if (latexListMap.empty())
			{
				vTrueNames.push_back(sampleFileName);
				vTrueLatexs.push_back(latexResult);
				vTrueGTLatexs.push_back("");
				//std::cout << "EMPTY" << std::endl;
			}
			else
			{
				if ((it = latexListMap.find(gtLatexID)) != latexListMap.end())
				{
					if (checkLatexNew(latexResult, it->second))
					{
						vTrueNames.push_back(sampleFileName);
						vTrueLatexs.push_back(latexResult);
						vTrueGTLatexs.push_back(it->second);
						std::cout << "True Sample" << std::endl;
					}
					else
					{
						//Parse to wrong Result
						vParseErrorNames.push_back(sampleFileName);
						vParseErrorLatexs.push_back(latexResult);
						vParseErrorGTLatexs.push_back(it->second);
						std::cout << "False Sample" << std::endl;
					}
					std::cout << "GT Latex : " << it->second << std::endl;
				}
				else
				{
					if (it == latexListMap.end())
					{
						//invalid samples, does not have a ground truth latex
						vLatexErrorNames.push_back(sampleFileName);
					}
				}
			}

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
		}
	}

	fs.close();
	int nsErr = vSymErrorNames.size(), nlatErr = vLatexErrorNames.size(),
		npErr = vParseErrorNames.size(), ntrue = vTrueNames.size();
	std::cout << "Total samples : " << nsErr + nlatErr + npErr + ntrue << std::endl;

	std::cout << "Symbol Error samples : " << nsErr << std::endl;
	std::cout << "No Ground Truth samples : " << nlatErr << std::endl;
	std::cout << "Accurancy : " << ntrue / float(npErr + ntrue) << "( " << ntrue << "/" << npErr + ntrue << " )" << std::endl;

	if (!saveResult)return 0;
	std::cout << "\n\n" << sepLine << std::endl;
	std::cout << "Start to save the result" << std::endl;

	std::string resultDir = "Result";
	//Check and Create the folder
	if (_access(resultDir.c_str(), 0) == -1)
		_mkdir(resultDir.c_str());

	std::string symErrorFName = resultDir + "/" + "symErrorSampleList.txt";
	std::string latexErrorFName = resultDir + "/" + "latexErrorSampleList.txt";
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
	fs.open(latexErrorFName, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " << latexErrorFName);

	for (size_t i = 0; i < vLatexErrorNames.size(); i++)
		fs << vLatexErrorNames[i] << std::endl;

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