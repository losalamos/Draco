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

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class FifiParser - 
//
// Date created :
// Purpose      :
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

    class ParseError : public std::runtime_error
    {
      public:
	ParseError(const std::string &str) : std::runtime_error(str) { }
    };
    
    typedef std::istream::pos_type pos_type;
    
    struct Position
    {
	int lineNo;
	pos_type filePosition;
    };
    
    typedef std::map<std::string, Position> KeywordPosMap;

    typedef int MaterialId;
    
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
    
    typedef std::map<MaterialId, MaterialInfo> MatInfoMap;

    enum DataType {REAL, CHARACTER, UNKNOWN};

    typedef std::map<std::string, DataType> DataTypeMap;
    
    // DATA
    
  private:

    static DataTypeMap dataTypeMap;

    static const int lineSize;

    // The curLinePos points to the beginning of the line
    // *** about to be read ***.
    
    mutable Position curLinePos;

    // The prevLinePos points to the previous value of
    // curLinePos.
    
    mutable Position prevLinePos;
    
    mutable std::istream &is;

    MatInfoMap materialInfoMap;

    // DISALLOWED DEFAULT METHODS
    
  private:
    
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

END_NS_XTM  // namespace XTM

#endif                          // __matprops_FifiParser_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/FifiParser.hh
//---------------------------------------------------------------------------//
