//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Token.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 13:15:29 MST 2003
 * \brief Tokens for use with simple parsers
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * This header file defines a class representing a lexical token, for use
 * in simple parsing systems for analysis codes.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Token_HH
#define CCS4_Token_HH

#include <string>
#include "ds++/Assert.hh"

namespace rtt_parser 
{

//-------------------------------------------------------------------------//
/*! 
 * \brief Token types recognized by a Token_Stream.
 */

enum Token_Type 
{
    END,    
    //!< The identifier <CODE>end</CODE>, denoting that the Parse_Table
    //!< should return control to its client.  Can be used to implement 
    //!< nested parse tables.
    
    EXIT,    
    //!< Denotes that the end of the Token_Stream has been reached.
    //!< The Token_Stream will continue to return EXIT indefinitely
    //!< once its end has been reached.
    
    KEYWORD, 
    //!< A sequence of one or more C++ identifiers separated by whitespace.
    
    REAL,
    //!< A valid C++ floating-point constant.
    
    INTEGER,
    //!< A valid C++ integer constant.
    
    STRING,
    //!< A valid C++ string constant. 
    
    ERROR,
    //!< The error token, indicating something wrong with the token stream.
    //!< For example, a file-based token stream would return this token if
    //!< the file failed to open.
    
    OTHER
    /*! A single character that does not belong to one of the regular token
     * types described above.
     */
};

//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Description of a token.
 *
 * Contains the type, value, and location of a parsed token.
 */

class Token
{
  public:
    inline Token(Token_Type t, std::string const &loc);
    inline Token(char c,       std::string const &loc);

    //! Construct a Token with specified type, text, and location.
    inline Token(Token_Type ty, std::string const &tx, std::string const &loc);
    
    //! Return the token type.
    Token_Type Type() const { return type_; }
    
    //! Return the token text.
    std::string Text() const { return text_; }
    
    //! Return the location information.
    std::string Location() const { return location_; }
    
  private:
    Token_Type type_;      //!< Type of this token
    std::string text_;     //!< Text of this token
    std::string location_; //!< Location information (such as file and line)

    bool check_class_invariants() const;
};

// For checking of assertions
bool Is_Integer_Text(char const *string);
bool Is_Keyword_Text(char const *string);
bool Is_Real_Text   (char const *string);
bool Is_String_Text (char const *string);

bool operator==(Token const &, Token const &);

//-------------------------------------------------------------------------//
/*! 
 * \param ty
 * Type of the Token.  Must equal KEYWORD, REAL, INTEGER, or STRING.
 *
 * \param tx
 * Text of the Token.
 *
 * \param loc
 * The token location.
 */

Token::Token(Token_Type const ty, 
	     std::string const &tx, 
	     std::string const &loc)
    : 
    type_(ty), text_(tx), location_(loc)
{
    Require(ty==KEYWORD || ty==REAL || ty==INTEGER || ty==STRING || ty==OTHER);
    Require(ty!=KEYWORD || Is_Keyword_Text(tx.c_str()));
    Require(ty!=REAL || Is_Real_Text(tx.c_str()));
    Require(ty!=INTEGER || Is_Integer_Text(tx.c_str()));
    Require(ty!=STRING || Is_String_Text(tx.c_str()));

    Ensure(check_class_invariants());
    Ensure(Type()==ty);
    Ensure(Text()==tx);
    Ensure(Location()==loc);
}

//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Construct a token of type OTHER.
 *
 * \param c
 * The token text (a single character)
 *
 * \param loc
 * The token location.
 *
 * \post <code>Type()==OTHER</code>
 * \post <code>Text().size()==1 && Text()[0]==c</code>
 * \post <code>Location()==loc</code>
 */

Token::Token(char c, std::string const &loc)
    : type_(OTHER), text_(1, c), location_(loc)
{
    Ensure(check_class_invariants());
    Ensure(Type()==OTHER);
    Ensure(Text().size()==1 && Text()[0]==c);
    Ensure(Location()==loc);
}

//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Construct a token of type END, EXIT, or ERROR.
 *
 * \param t
 * Token type to create; must be one of END, EXIT, or ERROR.
 *
 * \param loc
 * The token location
 *
 * These token types have no associated text.
 *
 * \pre <code>t==END || t==EXIT || t==ERROR</code>
 * \post <code>Type()==t</code>
 * \post <code>Text()==""</code>
 * \post <code>Location()==loc</code>
 */

Token::Token(Token_Type t, std::string const &loc)
    : type_(t), text_(), location_(loc)
{
    Require(t==END || t==EXIT || t==ERROR);

    Ensure(check_class_invariants());
    Ensure(Type()==t);
    Ensure(Text()=="");
    Ensure(Location()==loc);
}

}  // namespace rtt_parser

#endif  // CCS4_Token_HH
//--------------------------------------------------------------------//
//                      end of Token.hh
//--------------------------------------------------------------------//
