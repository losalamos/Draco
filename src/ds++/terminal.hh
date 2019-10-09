//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/terminal.hh
 * \author Ondrej Certik
 * \date   Sat Oct 05 2019
 * \brief  Terminal class that provides colored output.
 * \note   https://github.com/certik/terminal/blob/master/terminal.h
 *
 * This file is all platform independent, it contains the logic to build the
 * features that users need in a terminal application.
 *
 * The ANSI escape sequences used here are supported by all terminals (Linux,
 * macOS, Windows). All the functionality here must work on all platforms. The
 * Windows terminal is probably the most limiting, and so we restrict to the
 * capabilities that it supports, as documented at:
 *
 * https://docs.microsoft.com/en-us/windows/console/console-virtual-terminal-sequences
 *
 * \todo Consider an enum class for colors that derives from $LS_COLORS on
 *       Linux.  This would allow color selection based on users's terminal
 *       colors (e.g.: light vs dark scheme).
 */
//---------------------------------------------------------------------------//

#ifndef TERMINAL_H
#define TERMINAL_H

#include "terminal_base.hh"
#include "utf8_decode.hh"
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#define CTRL_KEY(k) ((k)&0x1f)
#define ALT_KEY(k) (k + 128)

namespace Term {

enum class style {
  reset = 0,
  bold = 1,
  dim = 2,
  italic = 3,
  underline = 4,
  blink = 5,
  rblink = 6,
  reversed = 7,
  conceal = 8,
  crossed = 9
};

enum class fg {
  black = 30,
  red = 31,
  green = 32,
  yellow = 33,
  blue = 34,
  magenta = 35,
  cyan = 36,
  gray = 37,
  reset = 39
};

enum class bg {
  black = 40,
  red = 41,
  green = 42,
  yellow = 43,
  blue = 44,
  magenta = 45,
  cyan = 46,
  gray = 47,
  reset = 49
};

enum class fgB {
  black = 90,
  red = 91,
  green = 92,
  yellow = 93,
  blue = 94,
  magenta = 95,
  cyan = 96,
  gray = 97
};

enum class bgB {
  black = 100,
  red = 101,
  green = 102,
  yellow = 103,
  blue = 104,
  magenta = 105,
  cyan = 106,
  gray = 107
};

inline std::string cursor_off() { return "\x1b[?25l"; }
inline std::string cursor_on() { return "\x1b[?25h"; }

// If an attempt is made to move the cursor out of the window, the result is
// undefined.
inline std::string move_cursor(size_t row, size_t col) {
  return "\x1b[" + std::to_string(row) + ";" + std::to_string(col) + "H";
}

// If an attempt is made to move the cursor to the right of the right margin,
// the cursor stops at the right margin.
inline std::string move_cursor_right(int col) {
  return "\x1b[" + std::to_string(col) + "C";
}

// If an attempt is made to move the cursor below the bottom margin, the cursor
// stops at the bottom margin.
inline std::string move_cursor_down(int row) {
  return "\x1b[" + std::to_string(row) + "B";
}

inline std::string cursor_position_report() { return "\x1b[6n"; }
inline std::string erase_to_eol() { return "\x1b[K"; }

enum Key {
  BACKSPACE = 1000,
  ENTER,
  ALT_ENTER,
  TAB,
  ARROW_LEFT,
  ARROW_RIGHT,
  ARROW_UP,
  ARROW_DOWN,
  CTRL_UP,
  CTRL_DOWN,
  CTRL_RIGHT,
  CTRL_LEFT,
  NUMERIC_5,
  DEL,
  HOME,
  INSERT,
  END,
  PAGE_UP,
  PAGE_DOWN,
  ESC,
  F1,
  F2,
  F3,
  F4,
  F5,
  F6,
  F7,
  F8,
  F9,
  F10,
  F11,
  F12,
};

class Terminal : public BaseTerminal {
  bool restore_screen_;

public:
  DLL_PUBLIC_dsxx static int term_initialized;

  Terminal(bool enable_keyboard = false, bool disable_ctrl_c = true)
      : BaseTerminal(enable_keyboard, disable_ctrl_c), restore_screen_{false} {
    Term::Terminal::term_initialized = 1;
  }
  virtual ~Terminal() { restore_screen(); }
  void restore_screen();
  void save_screen();
  inline void write(const std::string &s) const {
    std::cout << s << std::flush;
  }

  // Waits for a key press, translates escape codes
  int read_key() const;

  // If there was a key press, returns the translated key from escape codes,
  // otherwise returns 0. If the escape code is not supported, returns a
  // negative number.
  int read_key0() const;
  void get_cursor_position(int &rows, int &cols) const;

  // This function takes about 23ms, so it should only be used as a fallback
  void get_term_size_slow(int &rows, int &cols) const;
};

/*----------------------------------------------------------------------------*/

template <typename T> std::string color(T const value) {
  if (Term::Terminal::term_initialized > 0)
    return "\033[" + std::to_string(static_cast<int>(value)) + "m";
  else
    return std::string();
}

inline std::string remove_color(std::string const &colored_string) {
  std::regex color_regex("\033[" "[[:digit:]]+[m]");
  return std::regex_replace(colored_string, color_regex, "");
}

/*----------------------------------------------------------------------------*/

void codepoint_to_utf8(std::string &s, char32_t c);

// Converts an UTF8 string to UTF32.
std::u32string utf8_to_utf32(const std::string &s);

// Converts an UTF32 string to UTF8.
std::string utf32_to_utf8(const std::u32string &s);

/* Represents a rectangular window, as a 2D array of characters and their
 * attributes. The render method can convert this internal representation to a
 * string that when printed will show the Window on the screen.
 *
 * Note: the characters are represented by char32_t, representing their UTF-32
 * code point. The natural way to represent a character in a terminal would be
 * a "Unicode grapheme cluster", but due to a lack of a good library for C++
 * that could handle those, we simply use a Unicode code point as a character.
 */
class Window {
private:
  size_t x0, y0;               // top-left corner of the window on the screen
  size_t w, h;                 // width and height of the window
  std::vector<char32_t> chars; // the characters in row first order
  std::vector<fg> m_fg;
  std::vector<bg> m_bg;
  std::vector<style> m_style;

public:
  Window(size_t x0, size_t y0, size_t w, size_t h)
      : x0{x0}, y0{y0}, w{w}, h{h}, chars(w * h, ' '), m_fg(w * h, fg::reset),
        m_bg(w * h, bg::reset), m_style(w * h, style::reset) {}

  inline char32_t get_char(size_t x, size_t y) {
    return chars[(y - 1) * w + (x - 1)];
  }
  inline void set_char(size_t x, size_t y, char32_t c) {
    chars[(y - 1) * w + (x - 1)] = c;
  }
  inline fg get_fg(size_t x, size_t y) { return m_fg[(y - 1) * w + (x - 1)]; }
  inline void set_fg(size_t x, size_t y, fg c) {
    m_fg[(y - 1) * w + (x - 1)] = c;
  }
  inline bg get_bg(size_t x, size_t y) { return m_bg[(y - 1) * w + (x - 1)]; }
  inline void set_bg(size_t x, size_t y, bg c) {
    m_bg[(y - 1) * w + (x - 1)] = c;
  }
  inline style get_style(size_t x, size_t y) {
    return m_style[(y - 1) * w + (x - 1)];
  }
  inline void set_style(size_t x, size_t y, style c) {
    m_style[(y - 1) * w + (x - 1)] = c;
  }
  void print_str(int x, int y, const std::string &s);
  void fill_fg(int x1, int y1, int x2, int y2, fg color);
  void fill_bg(int x1, int y1, int x2, int y2, bg color);
  void fill_style(int x1, int y1, int x2, int y2, style color);
  inline void print_border(bool unicode = true) {
    print_rect(1, 1, w, h, unicode);
  }
  void print_rect(size_t x1, size_t y1, size_t x2, size_t y2,
                  bool unicode = true);
  void clear();
  std::string render();
};

} // namespace Term

#endif // TERMINAL_H

//---------------------------------------------------------------------------//
// end of ds++/terminal.hh
//---------------------------------------------------------------------------//
