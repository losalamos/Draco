//----------------------------------*-C++-*----------------------------------//
// FifiParser.cc
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "FifiParser.hh"

#include "ds++/Assert.hh"

#include <iostream>
using std::istream;
using std::cin;
#ifdef DEBUG_FIFIPARSER
using std::cerr;
#endif
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

using namespace rtt_matprops;

// C++ demands that we define our static members in the implementation file
// instead of the header file.

// The dataTypeMap is a map of keywords to their datatypes,
// e.g. keyword "ramg" is of type REAL.

FifiParser::DataTypeMap FifiParser::dataTypeMap;

// C++ demands that we define our static members in the implementation file
// instead of the header file, even the const ones.

// Fifi files are by definition only 80 characters long, period!!!

const int FifiParser::lineSize = 80;

//---------------------------------------------------------------------------//
// FifiParser c-tor:
//   Create a FifiParser from an open and valid input stream.
//---------------------------------------------------------------------------//

FifiParser::FifiParser( std::istream &is_)
    : is(is_)
{
    // Set the current position to the beginning of the file.
    
    is.seekg(0);
    curLinePos.lineNo = 1;
    curLinePos.filePosition = is.tellg();

    // There is no previous line position, so also set it to the beginning
    // of the file.  (Other methods rely on this behaviour.)
    
    prevLinePos = curLinePos;

    // Start the parsing process.
    // This will store the file positions associated with materials and
    // their respective keywords.
    
    parse();
}

//------------------------------------------------------------------------//
// getMaterialInfo:
//   Return a reference to the MaterialInfo referred to by materialId.
//------------------------------------------------------------------------//

const FifiParser::MaterialInfo &FifiParser::getMaterialInfo(
    MaterialId materialId) const
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

//------------------------------------------------------------------------//
// parseMaterial:
//   Parse the block of text beginning with the line after the
//   "material" keyword, and ending at the next "material" keyword.
//   Add all of the intervening keywords to the material's keyword-position
//   map.
//------------------------------------------------------------------------//

bool FifiParser::parseMaterial()
{
    string card;

    // This "card" (line) should only have the material id (an integer) on it.
    // It better be there.
    
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

#ifdef DEBUG_FIFIPARSER
    cerr << "parseMaterial: material id card: <" << card << ">" << endl << endl;
#endif

    // Convert the text of the "card" into an integer material id.
    
    istrstream iss(card.c_str());
    
    MaterialId matid;
    iss >> matid;
	    
    // Make sure it is a valid integer.

    if (!iss)
	throw ParseError(string("FifiParser::parseMaterial: ")
			 + "material id not found.");

#ifdef DEBUG_FIFIPARSER
    cerr << "found material: " << matid << endl;
#endif
    
    // Ensure that this a new material id.
	    
    if (materialInfoMap.find(matid) != materialInfoMap.end())
    {
	ostrstream os;
	os << "FifiParser::parseMaterial: "
	   << "Duplicate material: " << matid << ends;
	throw ParseError(os.str());
    }

    // Start loading up a keyword/position map for this new material.
    
    KeywordPosMap keywordPosMap;
    bool foundNextMaterial = false;
    
    while (getCard(card))
    {

#ifdef DEBUG_FIFIPARSER
	cerr << "parseMaterial: keyword card: <" << card << ">" << endl << endl;
#endif

	// This "card" (line) will contain some kind of keyword.
	// If it is a "material" keyword, then we have seen the whole
	// material block.  Otherwise we must store the keyword's position
	// into this material's keyword/position map.
	
	istrstream iss(card.c_str());
	string keyword;
	iss >> keyword;

	// Once we see the next material card, it is time to boogie.
	
	if (keyword == "material")
	{
	    // Rewind to beginning of the material line,
	    // so that the next read will re-read the material keyword.

	    setPosition(getPrevPosition());

	    foundNextMaterial = true;
	    break;
	}

	// If we successfully parsed a keyword clause, then
	// store off the position of the beginning of the data block,
	// just at the beginning of the keyword line.
	
	Position begOfKeywordBlock = getPrevPosition();

	// Parse the entire keyword block (usually data "cards").
	
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

	    // Insert the keyword's position into the material's map.
	    
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

    // We must clean up if we've found the next material or hit the
    // end of file.

    if (foundNextMaterial || is.eof())
    {

#ifdef DEBUG_FIFIPARSER
	cerr << "Inserting keyword map for material: " << matid << endl;
#endif

	// Cache the position map for this material into the FifiParser's
	// material position map.

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
}

// NonWhiteSpace: A functor to detect non-whitespace characters.

class NonWhiteSpace
{
  public:
    bool operator()(char c) { return !isspace(c); }
};

//---------------------------------------------------------------------------//
// getCard:
//    Get the next "card" (line) from the input file, saving
//    the position of the beginning of the next line into curLinePos.
//    This method returns the stream (to check its stream state).
//---------------------------------------------------------------------------//

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

#ifdef DEBUG_FIFIPARSER
	cerr << "Getting line: " << curLinePos.lineNo << endl;
#endif

	// Get the line into the buffer, up to but not including the '\n'
	
	std::getline(is, buf, '\n');

	// Is everything OK?
	
	if (!is)
	    return is;

	// Update the current line position to point to the beginning
	// of the **next line**.
	
	prevLinePos = curLinePos;
	curLinePos.filePosition = is.tellg();
	curLinePos.lineNo++;

	// No sense processing an empty line.
	
	if (buf.length() == 0)
	{
#ifdef DEBUG_FIFIPARSER
	    cerr << "Empty line" << endl;
#endif
	}
	else
	{
	    
#ifdef DEBUG_FIFIPARSER
	    cerr << "About to check for comments." << endl;
#endif
	
	    // Is this a comment line?  If so, keep on scanning in lines.
	
	    bool isComment =
		(comments.end() != find(comments.begin(), comments.end(),
					buf[0]));

	    if (!isComment)
	    {
		
#ifdef DEBUG_FIFIPARSER
		cerr << "About to check for non-whitespace data." << endl;
#endif
		
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

    // Set the return string to the buffer (up to the line size).
    
    str.assign(buf, 0, lineSize);

    return is;
}

//---------------------------------------------------------------------------//
// parse:
//   Perform the overall parsing, material text block by material
//   text block.
//---------------------------------------------------------------------------//

void FifiParser::parse()
{
    string card;

    // Go through the material "cards" (lines).
    // We should **only** see material keywords, since the rest of the
    // keywords and data will be swallowed up by parseMaterial().
    
    while (getCard(card))
    {
	istrstream iss(card.c_str());
	string keyword;
	iss >> keyword;

	if (!iss)
	    throw ParseError(string("FifiParser::parse: ")
			     + "keyword not found.");
	
#ifdef DEBUG_FIFIPARSER
	cerr << "parse found keyword: <" << keyword
	     << "> on line: " << curLinePos.lineNo - 1 << endl << endl;
#endif
	
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

    // If getCard() returned a non-zero stream state, then it had
    // better be due to an end of file!
    
    if (!is.eof())
	throw ParseError(string("FifiParser::parse: ")
			 + "stream error");
}

//---------------------------------------------------------------------------//
// getData:
//    Get the data from a data block, following a data keyword.
//---------------------------------------------------------------------------//

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

#ifdef DEBUG_FIFIPARSER
	cerr << "getData: data card: <" << card << ">" << endl << endl;
#endif

	int nvalsReadThisLine = 0;
	
	istrstream iss(card.c_str());
	while (iss >> val)
	{
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
    }

    // keyword rsmg can have a positive integer after it

#ifdef DEBUG_FIFIPARSER
    cerr << "Is keyword: <" << keyword << "> part of the rsmg family?" << endl;
#endif
    
    if (keyword == "rsmg")
    {
	return REAL;
    }
    else if (keyword.length() > 4 &&
	     keyword.compare(0, 4, "rsmg", 4) == 0)
    {
#ifdef DEBUG_FIFIPARSER
	cerr << "We have rsmg#" << endl;
#endif
	
	string numStr = keyword.substr(4, keyword.length()-4);

#ifdef DEBUG_FIFIPARSER
	cerr << "numStr: <" << numStr << ">" << endl;
#endif
	
	istrstream iss(numStr.c_str());

	// check to see if it is a valid integer after "rsmg"
	
	int num;
	iss >> num;
	if (!iss || num < 0)
	    return UNKNOWN;

	return REAL;
    }

#ifdef DEBUG_FIFIPARSER
    cerr << "Searching dataTypeMap for keyword: <" << keyword
	 << ">" << endl;
#endif
    
    DataTypeMap::iterator iter = dataTypeMap.find(keyword);

    // Did we find the keyword?
    
    if (iter == dataTypeMap.end())
	return UNKNOWN;

#ifdef DEBUG_FIFIPARSER
    cerr << "Good Keyword: " << keyword << endl;
#endif
    
    return (*iter).second;
}

void FifiParser::setPosition(const Position &position) const
{
    is.clear();
    is.seekg(position.filePosition);
    Position tmpLinePos = curLinePos;
    curLinePos = position;
#ifdef DEBUG_FIFIPARSER
    cerr << "Setting position beggining of line: "
	 << position.lineNo << endl;
#endif
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

#ifdef DEBUG_FIFIPARSER
    cerr << "Found keyword: <"
	 << keyword
	 << "> for material: " << matid << endl;
#endif

    getData(data);

    return true;
}

bool FifiParser::hasKeyword(MaterialId materialId, const string &keyword) const
{
    return getMaterialInfo(materialId).hasKeyword(keyword);
}

bool FifiParser::hasMaterial(MaterialId materialId) const
{
    return materialInfoMap.find(materialId) != materialInfoMap.end();
}

bool FifiParser::MaterialInfo::hasKeyword(const string &keyword) const
{
    return keywordPosMap.find(keyword) != keywordPosMap.end();
}
