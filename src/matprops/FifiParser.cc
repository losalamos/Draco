//----------------------------------*-C++-*----------------------------------//
// FifiParser.cc
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/FifiParser.hh"

#include "ds++/Assert.hh"

#include <iostream>
using std::istream;
using std::cin;
using std::cerr;
using std::cout;
using std::endl;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <map>

#include <algorithm>

#include <strstream>
using std::istrstream;
using std::ostrstream;
using std::ends;

using namespace XTM;

FifiParser::DataTypeMap FifiParser::dataTypeMap;

const int FifiParser::lineSize = 80;

FifiParser::FifiParser( std::istream &is_)
    : is(is_)
{
    is.seekg(0);
    curLinePos.lineNo = 1;
    curLinePos.filePosition = is.tellg();

    prevLinePos = curLinePos;
	
    parse();
}

typedef FifiParser::MaterialInfo MaterialInfo;

const MaterialInfo &FifiParser::getMaterialInfo(MaterialId materialId) const
{
    MatInfoMap::const_iterator matiter = materialInfoMap.find(materialId);

    if (matiter == materialInfoMap.end())
    {
	ostrstream os;
	os << "FifiParser::getMaterialInfo: "
	   << "material: " << materialId << " not found." << ends;
	
	throw ParseError(os.str());
    }

    return (*matiter).second;    
}

bool FifiParser::parseMaterial()
{
    string card;
    
    if (!getCard(card))
    {
	if (!is.eof())
	{
	    throw ParseError(string("FifiParser::parseMaterial: ")
			     + "stream error.");
	}
	else
	{
	    throw ParseError(string("FifiParser::parseMaterial: ")
			     + "material id not found at eof.");
	}
    }

    cerr << "parseMaterial: material id card: <" << card << ">" << endl << endl;

    istrstream iss(card.c_str());
    
    MaterialId matid;
    iss >> matid;
	    
    if (!iss)
	throw ParseError(string("FifiParser::parseMaterial: ")
			 + "material id not found.");

    cerr << "found material: " << matid << endl;
    
    // Ensure that this a new material id.
	    
    if (materialInfoMap.find(matid) != materialInfoMap.end())
    {
	ostrstream os;
	os << "FifiParser::parseMaterial: "
	   << "Duplicate material: " << matid << ends;
	throw ParseError(os.str());
    }

    KeywordPosMap keywordPosMap;
    bool foundNextMaterial = false;
    
    while (getCard(card))
    {
	cerr << "parseMaterial: keyword card: <" << card << ">" << endl << endl;

	istrstream iss(card.c_str());
	string keyword;
	iss >> keyword;

	// Once we see the next material card, it is time to buggy.
	
	if (keyword == "material")
	{
	    // rewind to beginning of the material line.

	    setPosition(getPrevPosition());

	    foundNextMaterial = true;
	    break;
	}

	// If we successfully parsed a keyword clause, then
	// store off the position of the beginning of the data block,
	// just at the beginning of the keyword line.
	
	Position begOfKeywordBlock = getPrevPosition();

	if (parseKeywordBlock(keyword))
	{
	    // Ensure that this keyword is new to this material.
	    
	    if (keywordPosMap.find(keyword) != keywordPosMap.end())
	    {
		ostrstream os;
		os << "FifiParser::parseMaterial: "
		   << "Duplicate keyword: " << keyword
		   << " for material: " << matid << ends;
		throw ParseError(os.str());
	    }

	    typedef KeywordPosMap::value_type value_type;
	    
	    std::pair<KeywordPosMap::iterator, bool> retval =
		keywordPosMap.insert(value_type(keyword, begOfKeywordBlock));
	    
	    // Make sure the insertion succeeded.
	
	    Assert(retval.second);	
	}
	else
	{
	    throw ParseError(string("FifiParser::parseMaterial: ")
			     + "parseKeywordBlock Unsuccessful.");
	}
    }

    if (foundNextMaterial || is.eof())
    {
	cerr << "Inserting keyword map for material: " << matid << endl;

	// Cache the position map for this material

	typedef MatInfoMap::value_type value_type;
	    
	std::pair<MatInfoMap::iterator, bool> retval =
	    materialInfoMap.insert(value_type(matid,
					      MaterialInfo(matid,
							   keywordPosMap)));
	    
	// Make sure the insertion succeeded.
	
	Assert(retval.second);
	    
	return true;
    }

    throw ParseError(string("FifiParser::parseMaterial: ")
		     + "stream error.");
    return false;
}

class NonWhiteSpace
{
  public:
    bool operator()(char c) { return !isspace(c); }
};
    
istream &FifiParser::getCard(string &str) const
{
    using std::find;
    using std::find_if;
    
    const string comments = "*#";

    string buf;
    buf.reserve(2*lineSize);

    bool foundData = false;

    while (!foundData)
    {
	cerr << "Getting line: " << curLinePos.lineNo << endl;
	
	std::getline(is, buf, '\n');
	
	if (!is)
	    return is;

	prevLinePos = curLinePos;
	curLinePos.filePosition = is.tellg();
	curLinePos.lineNo++;

	// No sense processing an empty line.
	
	if (buf.length() == 0)
	{
	    cerr << "Empty line" << endl;
	}
	else
	{
	    cerr << "About to check for comments." << endl;
	
	    // Is this a comment line?  If so, keep on scanning in lines.
	
	    bool isComment =
		(comments.end() != find(comments.begin(), comments.end(),
					buf[0]));

	    if (!isComment)
	    {
		cerr << "About to check for non-whitespace data." << endl;
		
		// Is this a blank line?  If so keep on scanning in lines.

		string::iterator bufFound = find_if(buf.begin(), buf.end(),
						    NonWhiteSpace());

		// A blank line is considered blank if the first lineSize
		// characters are whitespace.
	    
		foundData = (bufFound - buf.begin() < lineSize) &&
		    (bufFound != buf.end());
	    }
	}
	
    }

    // Set the return string to the buffer.
    
    str.assign(buf, 0, lineSize);

    return is;
}

void FifiParser::parse()
{
    string card;
    
    while (getCard(card))
    {
	istrstream iss(card.c_str());
	string keyword;
	iss >> keyword;

	if (!iss)
	    throw ParseError(string("FifiParser::parse: ")
			     + "keyword not found.");
	
	cerr << "parse found keyword: <" << keyword
	     << "> on line: " << curLinePos.lineNo - 1 << endl << endl;
	
	if (keyword == "material")
	{
	    parseMaterial();
	}
	else
	{
	    throw ParseError(string("FifiParser::parse: ")
			     + "Unexpected keyword: "
			     + keyword + ".");
	}
    }

    if (!is.eof())
	throw ParseError(string("FifiParser::parse: ")
			 + "stream error");
}

bool FifiParser::getData(vector<double> &data_) const
{
    // The file is assumed positioned to the start of the data.
    
    double val;
    vector<double> data;
    
    bool done = false;

    while (!done)
    {
	string card;

	if (!getCard(card))
	{
	    if (is.eof())
	    {
		// Set the answer.
    
		data_ = data;

		return data.size() > 0;
	    }
	    else
		throw ParseError(string("FifiParser::getData: ")
				 + "stream error");
	}

	cerr << "getData: data card: <" << card << ">" << endl << endl;

	int nvalsReadThisLine = 0;
	
	istrstream iss(card.c_str());
	while (iss >> val)
	{
	    cerr << "read data value: " << val << endl;
	    data.push_back(val);
	    nvalsReadThisLine++;
	}

	if (nvalsReadThisLine == 0)
	    done = true;
    }

    // We just found the next keyword.
    // We must reset the position to the beginning of the keyword line.

    setPosition(getPrevPosition());

    // Set the answer.
    
    data_ = data;

    return data.size() > 0;
}

bool FifiParser::getData(vector<string> &data_) const
{
    // The file is assumed positioned to the start of the data.
    
    string val;
    vector<string> data;
    
    // There is only one CHARACTER data card.
    
    string card;

    if (!getCard(card))
	throw ParseError(string("FifiParser::getData: ")
			 + "stream error");

    istrstream iss(card.c_str());
    while (iss >> val)
    {
	data.push_back(val);
    }

    // Set the answer.
    
    data_ = data;

    return data.size() > 0;
}

bool FifiParser::parseKeywordBlock(const string &keyword)
{

    DataType dataType = getDataType(keyword);

    bool success;

    vector<double> ddata;
    vector<string> cdata;
    
    switch (dataType)
    {
    case REAL:
	success = getData(ddata);
	break;
    case CHARACTER:
	success = getData(cdata);
	break;
    default:
	throw ParseError(string("FifiParser::parseKeywordBlock: ")
			 + "Unknown keyword: " + keyword + ".");
    }

    return success;
}

FifiParser::DataType FifiParser::getDataType(const std::string &keyword)
{
    typedef DataTypeMap::value_type value_type;

    // If this is the first time through, load up the table.
    
    static bool first = true;

    if (first)
    {
	first = false;

	// rsmg is handled separately

	dataTypeMap.insert(value_type("general", CHARACTER));
	dataTypeMap.insert(value_type("tgrid",   REAL));
	dataTypeMap.insert(value_type("rgrid",   REAL));
	dataTypeMap.insert(value_type("hnugrid", REAL));
	dataTypeMap.insert(value_type("ramg",    REAL));
	dataTypeMap.insert(value_type("rtmg",    REAL));
	dataTypeMap.insert(value_type("pmg",     REAL));
	dataTypeMap.insert(value_type("rgray",   REAL));
	dataTypeMap.insert(value_type("pgray",   REAL));
	dataTypeMap.insert(value_type("p",       REAL));
	dataTypeMap.insert(value_type("e",       REAL));
	dataTypeMap.insert(value_type("tfree",   REAL));
	dataTypeMap.insert(value_type("pelect",  REAL));
	dataTypeMap.insert(value_type("eelect",  REAL));
	dataTypeMap.insert(value_type("pnuc",    REAL));
	dataTypeMap.insert(value_type("enuc",    REAL));

	cerr << "dataTypeMap created" << endl;
    }

    // keyword rsmg can have a positive integer after it

    cerr << "Is keyword: <" << keyword << "> part of the rsmg family?" << endl;
    
    if (keyword == "rsmg")
    {
	return REAL;
    }
    else if (keyword.length() > 4 &&
	     keyword.compare(0, 4, "rsmg", 4) == 0)
    {
	cerr << "We have rsmg#" << endl;
	
	string numStr = keyword.substr(4, keyword.length()-4);

	cerr << "numStr: <" << numStr << ">" << endl;
	
	istrstream iss(numStr.c_str());

	// check to see if it is a valid integer after "rsmg"
	
	int num;
	iss >> num;
	if (!iss || num < 0)
	    return UNKNOWN;

	return REAL;
    }

    cerr << "Searching dataTypeMap for keyword: <" << keyword
	 << ">" << endl;
    
    DataTypeMap::iterator iter = dataTypeMap.find(keyword);

    // Did we find the keyword?
    
    if (iter == dataTypeMap.end())
	return UNKNOWN;

    cerr << "Good Keyword: " << keyword << endl;
    
    return (*iter).second;
}

void FifiParser::setPosition(const Position &position) const
{
    is.clear();
    is.seekg(position.filePosition);
    Position tmpLinePos = curLinePos;
    curLinePos = position;
    cerr << "Setting position beggining of line: "
	 << position.lineNo << endl;
    prevLinePos = tmpLinePos;
}

const FifiParser::Position &FifiParser::getPosition(MaterialId materialId,
					  const string &keyword) const
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    return getPosition(matInfo, keyword);
}

const FifiParser::Position &FifiParser::getPosition(const MaterialInfo &matInfo,
					  const string &keyword) const
{
    const KeywordPosMap &keywordPosMap = matInfo.keywordPosMap;
    KeywordPosMap::const_iterator keyiter = keywordPosMap.find(keyword);
    
    if (keyiter == keywordPosMap.end())
    {
	ostrstream os;
	os << "FifiParser::getPosition: "
	   << "keyword: " << keyword << " not found for "
	   << "material: " << matInfo.matid << "." << ends;
	
	throw ParseError(os.str());
    }

    return (*keyiter).second;
}

bool FifiParser::getData(MaterialId matid, const string &keyword,
			 vector<double> &data)
{
    setPosition(getPosition(matid, keyword));

	// Get the actual keyword card

    string str;
    if (!getCard(str))
    {
	ostrstream os;
	os << "FifiParser::getData:"
	   << " stream error for material: " << matid
	   << " and keyword: <"
	   << keyword
	   << ">" << endl;
	throw ParseError(os.str());
    }

    cerr << "Found keyword: <"
	 << keyword
	 << "> for material: " << matid << endl;

    getData(data);

    return true;
}

bool FifiParser::hasKeyword(MaterialId materialId, const string &keyword) const
{
    return getMaterialInfo(materialId).hasKeyword(keyword);
}

bool FifiParser::MaterialInfo::hasKeyword(const string &keyword) const
{
    return keywordPosMap.find(keyword) != keywordPosMap.end();
}

