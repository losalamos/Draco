//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/DracoStrings.hh
 * \author Kelly G. Thompson <kgt@lanl.gov
 * \date   Wednesday, Aug 23, 2017, 12:48 pm
 * \brief  Enscapulates common string manipulations.
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoStrings_hh
#define rtt_dsxx_DracoStrings_hh

#include "ds++/config.h"
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace rtt_dsxx {

//----------------------------------------------------------------------------//
/*!
 * \brief Convert a number into a string
 *
 * \param[in] num The integral or float-point type number that will be converted
 *                to a string.
 * \param[in] precision How many digits will be preserved? (default: 23)
 * \return A string representation of the numeric value
 *
 * \tparam T type of number that will be converted to a string. Typically a
 *           double or int.
 *
 * \sa http://public.research.att.com/~bs/bs_faq2.html
 */
template <typename T>
std::string to_string(T const num, unsigned int const precision = 23) {
  std::ostringstream s;
  s.precision(precision);
  s << num;
  return s.str();
}

//----------------------------------------------------------------------------//
/*!
 * \brief trim whitespace (or other characters) from before and after main text.
 *
 * \param[in] str The string that will be processed
 * \param[in] whitespace A set of characters that will be removed.
 *                 (default: " \t")
 * \return A new, probably shortened, string without unwanted leading/training
 *         characters.
 */
DLL_PUBLIC_dsxx std::string trim(std::string const &str,
                                 std::string const &whitespace = " \t");

//----------------------------------------------------------------------------//
/*!
 * \brief Removes all specified characters from a string.
 *
 * \param[in] str The string that will be processed
 * \param[in] chars_to_remove A set of characters (as a std::string) that will
 *         be removed.
 * \return A new, possibly shortened, string that does not contain the unwanted
 *         characters.
 */
DLL_PUBLIC_dsxx std::string prune(std::string const &orig_str,
                                  std::string const &chars_to_remove);

//----------------------------------------------------------------------------//
/*!
 * \brief Split a string into a vector<string> using a specified delimiter.
 *
 * \param[in] str The string that will be split
 * \param[in] delimiters Set of character that separate each 'word'.
 *                 (default: " ")
 * \param[in] keepEmptyStrings Should empty values be saved in the result?
 *                 (default: false)
 * \return a vector of strings. The delimiter is always removed.
 */
DLL_PUBLIC_dsxx std::vector<std::string>
tokenize(std::string const &str, std::string const &delimiters = " ",
         bool keepEmptyStrings = false);

//----------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a floating-point type or an integral type
 *        without error checking.
 *
 * \param[in] str The string that contains a number.
 * \return A numeric value.
 *
 * \tparam The template type will be deduced based on the return type ("-> T").
 */
template <typename T> auto parse_number_impl(std::string const &str) -> T;

// specializations for these types are devined in DracoStrings.cc
template <>
DLL_PUBLIC_dsxx auto parse_number_impl<int>(std::string const &str) -> int;
template <>
DLL_PUBLIC_dsxx auto parse_number_impl<long>(std::string const &str) -> long;
template <>
DLL_PUBLIC_dsxx auto parse_number_impl<unsigned long>(std::string const &str)
    -> unsigned long;
template <>
DLL_PUBLIC_dsxx auto parse_number_impl<float>(std::string const &str) -> float;
template <>
DLL_PUBLIC_dsxx auto parse_number_impl<double>(std::string const &str)
    -> double;

//----------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a floating-point type or an integral type with
 *        error checking.
 *
 * \param[in] str The string that contains a number.
 * \param[in] verbose Should the function print conversion message warnings
 *                 (default: true).
 * \return A numeric value.
 *
 * \tparam The template type will be deduced based on the return type ("-> T").
 *
 * These functions will throw exceptions if the conversion is invalid or if the
 * converted value is out-of-range. This generic function calls
 * unchecked parse_number_impl specialization that is appropriate.
 *
 * Consider catching thrown conversion errors using code similar to this:

 * \code
 try {
    parse_number<T>(str);
  } catch (std::invalid_argument &e) {
    if no conversion could be performed
    std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
              << "No valid conversion from string to a numeric value could be "
              << "found.\n"
              << "\tstring = \"" << str << "\"\n" << std::endl;
    throw e;
  } catch (std::out_of_range &e) {
    if the converted value would fall out of the range of the result type
    or if the underlying function (std::strtol or std::strtoull) sets
    errno to ERANGE.
    std::cerr << "\n==ERROR==\nrtt_dsxx::parse_number:: "
              << "Type conversion from string to a numeric value resulted in "
              << "a value that is out of range.\n"
              << "\tstring = \"" << str << "\"\n" << std::endl;
    throw;
  } catch (...) {
    everything else
    Insist(false, "Unknown failures in call to std::sto[x] from parse_number.");
  }
 * \endcode
 */
template <typename T>
auto parse_number(std::string const &str, bool verbose = true) -> T;

//----------------------------------------------------------------------------//
/*!
 * \brief Convert a string into a vector of floating-point or an integral
 *        values.
 *
 * \param[in] str The string that contains a number.
 * \param[in] range_symbols Parenthesis or brances that mark the begining or end
 *                 of the value range. (default: "{}")
 * \param[in] delimiter A character that separates each numeric entry.
 *                 (default: ",")
 * \return A vector of numeric values.
 *
 * \tparam plain-old-data type for storing numeric values, such as double,
 *                 unsigned, float or int.
 *
 * std::string foo = "{ 1.0, 2.0, 3.0 }" will be converted to
 * std::vector<double> bar = { 1.0, 2.0, 3.0 }
 */
template <typename T>
std::vector<T> string_to_numvec(std::string const &str,
                                std::string const &range_symbols = "{}",
                                std::string const &delimiters = ",");

//----------------------------------------------------------------------------//
/*!
 * \brief Parse an ostream and create a map containing [word]:[num_occurances]
 *
 * \param[in] data An ostringstream that will be parsed
 * \param[in] verbose Echo the data map to stdout when this function is run
 *                 (default: false).
 * \return a map<string,uint> that contains [word]:[num_occurances]
 */
DLL_PUBLIC_dsxx std::map<std::string, unsigned>
get_word_count(std::ostringstream const &data, bool verbose = false);

//----------------------------------------------------------------------------//
/*!
 * \brief Parse a file and create a map containing [word]:[num_occurances]
 *
 * \param[in] filename The name of the file (full path recommended) that will be
 *                 parsed.
 * \param[in] verbose Echo the data map to stdout when this function is run
 *                 (default: false).
 * \return a map<string,uint> that contains [word]:[num_occurances]
 */
DLL_PUBLIC_dsxx std::map<std::string, unsigned>
get_word_count(std::string const &filename, bool verbose = false);

} // namespace rtt_dsxx

// Template functions
#include "DracoStrings.i.hh"

#endif // rtt_dsxx_DracoStrings_hh

//---------------------------------------------------------------------------//
// end of DracoStrings.hh
//---------------------------------------------------------------------------//
