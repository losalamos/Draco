//----------------------------------*-C++-*----------------------------------//
// FifiParser.hh
// Randy M. Roberts
// Mon May  4 11:12:47 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_FifiParser_hh__
#define __matprops_FifiParser_hh__

#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace rtt_matprops {
    
//===========================================================================//
// class FifiParser - 
//
// Date created :
// Purpose      : This class is responsible for reading a Fifi file,
//                storing the file positions associated with materials and
//                their keywords, and retrieving the same.
//                This class is used by the FifiMatPropsReader class.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class FifiParser
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    //========================================================================//
    // class ParseError:
    //   A handy class in which to throw an exception.
    //========================================================================//
    
    class ParseError : public std::runtime_error
    {
      public:
	ParseError(const std::string &str) : std::runtime_error(str) { }
    };

  private:
    
    // A shorthand name for a long type.
    
    typedef std::istream::pos_type pos_type;

    // A nested struct to store the linenumber and position of the
    // beginning of a line.
    
    struct Position
    {
	int lineNo;
	pos_type filePosition;
    };
    
    // A mapping of a keyword, within a material, to the beginning
    // of the line containing the keyword.
    
    typedef std::map<std::string, Position> KeywordPosMap;

    // How to reference a material.
    
    typedef int MaterialId;
    
    //========================================================================//
    // class MaterialInfo:
    //    A nested struct that contains the material id and a mapping
    //    of all the line positions for the keywords within the material.
    //========================================================================//
    
    struct MaterialInfo
    {
	MaterialId          matid;
	KeywordPosMap       keywordPosMap;
	
	MaterialInfo()
	{
	    //*empty*
	}
	MaterialInfo(MaterialId matid_, const KeywordPosMap &keywordPosMap_)
	    : matid(matid_), keywordPosMap(keywordPosMap_)
	{
	    //*empty*
	}
	bool hasKeyword(const std::string &keyword) const;
    };

    // A mapping of a material id to a MaterialInfo object.
    
    typedef std::map<MaterialId, MaterialInfo> MatInfoMap;

    // Keywords pertain to two types of data, REAL, and CHARACTER.
    // If we've never seen this keyword then it is considered of UNKNOWN type.
    
    enum DataType {REAL, CHARACTER, UNKNOWN};

    // A mapping of known keywords to their respective types.
    // This mapping is set up once at the beginning of parsing.
    
    typedef std::map<std::string, DataType> DataTypeMap;
    
    // DATA
    
  private:

    // A mapping of known keywords to their respective types.
    // This mapping is set up once at the beginning of parsing.
    
    static DataTypeMap dataTypeMap;

    // Fifi files are only 80 characters long.
    
    static const int lineSize;

    // The curLinePos points to the beginning of the line
    // *** about to be read ***.
    
    mutable Position curLinePos;

    // The prevLinePos points to the previous value of
    // curLinePos.
    
    mutable Position prevLinePos;

    // The input stream that we will parse.
    
    mutable std::istream &is;

    // The mapping of a material id to a MaterialInfo object
    // used to find the file position of the keywords within the material
    // block.
    
    MatInfoMap materialInfoMap;

    // DISALLOWED DEFAULT METHODS
    
  private:

    // We do not allow copies of the FifiParser.
    
    FifiParser(const FifiParser &rhs);
    FifiParser& operator=(const FifiParser &rhs);

  public:

    // CREATORS
    
    FifiParser(std::istream &is_);

    // MANIPULATORS
    
    bool getData(MaterialId matid, const std::string &keyword,
		 std::vector<double> &data);
    
    // ACCESSORS

    bool hasKeyword(MaterialId materialId, const std::string &keyword) const;
    bool hasMaterial(MaterialId materialId) const;
    
  private:
    
    // IMPLEMENTATION

  private:

    const MaterialInfo &getMaterialInfo(MaterialId materialId) const;

    void setPosition(const Position &position) const;
    
    const Position &getPrevPosition() const { return prevLinePos; }

    const Position &getPosition() const { return curLinePos; }

    const Position &getPosition(MaterialId materialId,
				const std::string &keyword) const;
	
    const Position &getPosition(const MaterialInfo &materialInfo,
				const std::string &keyword) const;
	
    void parse();

    std::istream &getCard(std::string &str) const;

    bool parseMaterial();

    bool parseKeywordBlock(const std::string &keyword);

    static DataType getDataType(const std::string &keyword);

    bool getData(std::vector<double> &data_) const;

    bool getData(std::vector<std::string> &data_) const;

};

} // end of rtt_matprops namespace

#endif                          // __matprops_FifiParser_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/FifiParser.hh
//---------------------------------------------------------------------------//
