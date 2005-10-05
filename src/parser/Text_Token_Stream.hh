//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Text_Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 13:15:29 MST 2003
 * \brief Definition of the Text_Token_Stream class.
 * \note Copyright @ 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/12/03): 
 *    Fix indentation and comments.
 * 2) kgbudge (03/18/03): 
 *    Fix bugs and add capability to scan comments.
 * 3) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Text_Token_Stream_HH
#define CCS4_Text_Token_Stream_HH

#include <set>
#include "Token_Stream.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Abstract text-based token stream for simple parsers.
 *
 * A Text_Token_Stream obtains its stream of tokens by scanning a stream
 * of text characters, supplied through the getc member of its child
 * classes.
 *
 * C and C++ comments are treated as breaking whitespace.
 *
 * Null characters are not permitted in the character stream.  They are used
 * internally to indicate the end of file or an error condition.
 *
 * \invariant The whitespace list is a superset of the nonbreaking whitespace
 * list.
 */

class Text_Token_Stream : public Token_Stream
{
  protected:

    //! Construct a Text_Token_Stream.
    Text_Token_Stream();

    //! Construct a Text_Token_Stream.
    Text_Token_Stream(std::set<char> const &);
    
  public:

    unsigned Line() const { Ensure(line>0); return line; }
    //!< Return the current line in the text stream.

    std::set<char> const &Whitespace() const { return whitespace; }

    //! Does the Token_Stream consider \c c to be whitespace?
    bool is_whitespace(char c) const;

    //! Does the Token_Stream consider <i>c</i> to be nonbreaking
    bool is_nb_whitespace(char c) const;
       
    //! Check the class invariants.
    bool check_class_invariants() const;

    //! The default whitespace definition
    static std::set<char> const default_whitespace;

  protected:

    // IMPLEMENTATION
    
    //! Scan the next token.
    virtual Token fill();
    
    //! Push a character onto the back of the character queue.
    void character_push_back(char c);

    virtual void fill_character_buffer() = 0;
    //!< Move one or more characters from the text stream into the character
    //!< buffer. 
    
    virtual bool error() const = 0;
    //!< Has an I/O error occurred while reading the text stream?
    
    virtual bool end() const = 0;
    //!< Has the end of the text stream been reached?
   
    //! Returns a location string.    
    virtual std::string location() const = 0;

    //! Rewind the file token stream.
    virtual void Rewind() = 0;

    //! Pop a character off the internal buffer. 
    char pop_char();
    //! Peek ahead at the internal buffer. 
    char peek(unsigned pos = 0);

    //! Skip initial whitespace, if any.
    void eat_whitespace();

    // The following scan_ functions are for numeric scanning.  The names
    // reflect the context-free grammar given by Stroustrup in appendix A 
    // of _The C++ Programming Language_.  However, we do not presently
    // recognize type suffixes on either integers or floats.
    unsigned scan_floating_literal();
    unsigned scan_digit_sequence(unsigned &);
    unsigned scan_exponent_part(unsigned &);
    unsigned scan_fractional_constant(unsigned &);
    
    unsigned scan_integer_literal();
    unsigned scan_decimal_literal(unsigned &);
    unsigned scan_hexadecimal_literal(unsigned &);
    unsigned scan_octal_literal(unsigned &);

  private:    

    // IMPLEMENTATION

    // DATA

    std::deque<char> buffer;
    //!< Character buffer. Refilled as needed using fill_character_buffer()

    std::set<char> whitespace; //!< the whitespace character list
         
    unsigned line;  //!< Current line in input file.
};

}  // namespace rtt_parser

#endif  // CCS4_Text_Token_Stream_HH
//--------------------------------------------------------------------//
//                      end of Text_Token_Stream.hh
//--------------------------------------------------------------------//
