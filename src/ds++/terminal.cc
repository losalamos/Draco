//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/terminal.cc
 * \author Ondrej Certik
 * \date   Sat Oct 05 2019
 * \brief  Terminal class that provides colored output.
 * \note   https://github.com/certik/terminal/blob/master/terminal.h */
//---------------------------------------------------------------------------//

#include "terminal.hh"

// Default initialization of global data
int Term::Terminal::term_initialized = 0;

namespace Term {

void Terminal::restore_screen() {
  if (restore_screen_) {
    write("\033[?1049l"); // restore screen
    write("\033"
          "8"); // restore current cursor position
    restore_screen_ = false;
  }
}

void Terminal::save_screen() {
  restore_screen_ = true;
  write("\033"
        "7");           // save current cursor position
  write("\033[?1049h"); // save screen
}

// Waits for a key press, translates escape codes
int Terminal::read_key() const {
  int key;
  while ((key = read_key0()) == 0) {
  }
  return key;
}

// If there was a key press, returns the translated key from escape codes,
// otherwise returns 0. If the escape code is not supported, returns a negative
// number.
int Terminal::read_key0() const {
  char c;
  if (!read_raw(&c))
    return 0;

  if (c == '\x1b') {
    char seq[4];

    if (!read_raw(&seq[0]))
      return Key::ESC;
    if (!read_raw(&seq[1])) {
      if (seq[0] >= 'a' && seq[0] <= 'z') {
        // gnome-term, Windows Console
        return ALT_KEY(seq[0]);
      }
      if (seq[0] == 13) {
        // gnome-term
        return Key::ALT_ENTER;
      }
      return -1;
    }

    if (seq[0] == '[') {
      if (seq[1] >= '0' && seq[1] <= '9') {
        if (!read_raw(&seq[2])) {
          return -2;
        }
        if (seq[2] == '~') {
          switch (seq[1]) {
          case '1':
            return Key::HOME;
          case '2':
            return Key::INSERT;
          case '3':
            return Key::DEL;
          case '4':
            return Key::END;
          case '5':
            return Key::PAGE_UP;
          case '6':
            return Key::PAGE_DOWN;
          case '7':
            return Key::HOME;
          case '8':
            return Key::END;
          }
        } else if (seq[2] == ';') {
          if (seq[1] == '1') {
            if (!read_raw(&seq[2])) {
              return -10;
            }
            if (!read_raw(&seq[3])) {
              return -11;
            }
            if (seq[2] == '5') {
              switch (seq[3]) {
              case 'A':
                return Key::CTRL_UP;
              case 'B':
                return Key::CTRL_DOWN;
              case 'C':
                return Key::CTRL_RIGHT;
              case 'D':
                return Key::CTRL_LEFT;
              }
            }
            return -12;
          }
        } else {
          if (seq[2] >= '0' && seq[2] <= '9') {
            if (!read_raw(&seq[3])) {
              return -3;
            }
            if (seq[3] == '~') {
              if (seq[1] == '1') {
                switch (seq[2]) {
                case '5':
                  return Key::F5;
                case '7':
                  return Key::F6;
                case '8':
                  return Key::F7;
                case '9':
                  return Key::F8;
                }
              } else if (seq[1] == '2') {
                switch (seq[2]) {
                case '0':
                  return Key::F9;
                case '1':
                  return Key::F10;
                case '3':
                  return Key::F11;
                case '4':
                  return Key::F12;
                }
              }
            }
          }
        }
      } else {
        switch (seq[1]) {
        case 'A':
          return Key::ARROW_UP;
        case 'B':
          return Key::ARROW_DOWN;
        case 'C':
          return Key::ARROW_RIGHT;
        case 'D':
          return Key::ARROW_LEFT;
        case 'E':
          return Key::NUMERIC_5;
        case 'H':
          return Key::HOME;
        case 'F':
          return Key::END;
        }
      }
    } else if (seq[0] == 'O') {
      switch (seq[1]) {
      case 'F':
        return Key::END;
      case 'H':
        return Key::HOME;
      case 'P':
        return Key::F1;
      case 'Q':
        return Key::F2;
      case 'R':
        return Key::F3;
      case 'S':
        return Key::F4;
      }
    }
    return -4;
  } else {
    switch (c) {
    case 9:
      return Key::TAB;
    case 13:
      return Key::ENTER;
    case 127:
      return Key::BACKSPACE;
    }
    if (c == -61) {
      if (!read_raw(&c)) {
        return -8;
      } else {
        if (c >= -95 && c <= -70) {
          // xterm
          return ALT_KEY(c + 'a' - (-95));
        }
        return -9;
      }
    } else if (c == -62) {
      if (!read_raw(&c)) {
        return -10;
      } else {
        if (c == -115) {
          // xterm
          return Key::ALT_ENTER;
        }
        return -11;
      }
    }
    return c;
  }
}

void Terminal::get_cursor_position(int &rows, int &cols) const {
  char buf[32];
  unsigned int i = 0;
  write(cursor_position_report());
  while (i < sizeof(buf) - 1) {
    while (!read_raw(&buf[i])) {
    };
    if (buf[i] == 'R')
      break;
    i++;
  }
  buf[i] = '\0';
  if (i < 5) {
    throw std::runtime_error("get_cursor_position(): too short response");
  }
  if (buf[0] != '\x1b' || buf[1] != '[') {
    throw std::runtime_error("get_cursor_position(): Invalid response");
  }
  if (sscanf(&buf[2], "%d;%d", &rows, &cols) != 2) {
    throw std::runtime_error("get_cursor_position(): Invalid response");
  }
}

// This function takes about 23ms, so it should only be used as a fall back
void Terminal::get_term_size_slow(int &rows, int &cols) const {
  struct CursorOff {
    const Terminal &term;
    CursorOff(const Terminal &term) : term{term} { term.write(cursor_off()); }
    ~CursorOff() { term.write(cursor_on()); }
  };
  CursorOff cursor_off(*this);
  int old_row, old_col;
  get_cursor_position(old_row, old_col);
  write(move_cursor_right(999) + move_cursor_down(999));
  get_cursor_position(rows, cols);
  write(move_cursor(old_row, old_col));
}

/*----------------------------------------------------------------------------*/

void codepoint_to_utf8(std::string &s, char32_t c) {
  int nbytes;
  if (c < 0x80) {
    nbytes = 1;
  } else if (c < 0x800) {
    nbytes = 2;
  } else if (c < 0x10000) {
    nbytes = 3;
  } else if (c <= 0x0010FFFF) {
    nbytes = 4;
  } else {
    throw std::runtime_error("Invalid UTF32 codepoint.");
  }
  char32_t u1('x'), u2('x'), u3('x'), u4('x');
  static const unsigned char mask[4] = {0x00, 0xC0, 0xE0, 0xF0};
  switch (nbytes) {
  case 4:
    u4 = ((c | 0x80) & 0xBF);
    c >>= 6; /* fall through */
  case 3:
    u3 = ((c | 0x80) & 0xBF);
    c >>= 6; /* fall through */
  case 2:
    u2 = ((c | 0x80) & 0xBF);
    c >>= 6; /* fall through */
  case 1:
    u1 = (c | mask[nbytes - 1]);
  }
  switch (nbytes) {
  case 1:
    s.push_back(static_cast<char>(u1));
    break;
  case 2:
    s.push_back(static_cast<char>(u1));
    s.push_back(static_cast<char>(u2));
    break;
  case 3:
    s.push_back(static_cast<char>(u1));
    s.push_back(static_cast<char>(u2));
    s.push_back(static_cast<char>(u3));
    break;
  case 4:
    s.push_back(static_cast<char>(u1));
    s.push_back(static_cast<char>(u2));
    s.push_back(static_cast<char>(u3));
    s.push_back(static_cast<char>(u4));
    break;
  }
}

/*----------------------------------------------------------------------------*/

// Converts an UTF8 string to UTF32.
std::u32string utf8_to_utf32(const std::string &s) {
  uint32_t codepoint;
  uint8_t state = UTF8_ACCEPT;
  std::u32string r;
  for (size_t i = 0; i < s.size(); i++) {
    state = utf8_decode_step(state, s[i], &codepoint);
    if (state == UTF8_ACCEPT) {
      r.push_back(codepoint);
    }
    if (state == UTF8_REJECT) {
      throw std::runtime_error("Invalid byte in UTF8 encoded string");
    }
  }
  if (state != UTF8_ACCEPT) {
    throw std::runtime_error("Expected more bytes in UTF8 encoded string");
  }
  return r;
}

// Converts an UTF32 string to UTF8.
std::string utf32_to_utf8(const std::u32string &s) {
  std::string r;
  for (size_t i = 0; i < s.size(); i++) {
    codepoint_to_utf8(r, s[i]);
  }
  return r;
}

void Window::print_str(int x, int y, const std::string &s) {
  std::u32string s2 = utf8_to_utf32(s);
  for (size_t i = 0; i < s2.size(); i++) {
    size_t xpos = x + i;
    if (xpos < w) {
      set_char(xpos, y, s2[i]);
    } else {
      // String is out of the window
      return;
    }
  }
}

void Window::fill_fg(int x1, int y1, int x2, int y2, fg color) {
  for (int j = y1; j <= y2; j++) {
    for (int i = x1; i <= x2; i++) {
      set_fg(i, j, color);
    }
  }
}

void Window::fill_bg(int x1, int y1, int x2, int y2, bg color) {
  for (int j = y1; j <= y2; j++) {
    for (int i = x1; i <= x2; i++) {
      set_bg(i, j, color);
    }
  }
}

void Window::fill_style(int x1, int y1, int x2, int y2, style color) {
  for (int j = y1; j <= y2; j++) {
    for (int i = x1; i <= x2; i++) {
      set_style(i, j, color);
    }
  }
}

void Window::print_rect(size_t x1, size_t y1, size_t x2, size_t y2,
                        bool unicode) {
  std::u32string border = utf8_to_utf32("│─┌┐└┘");
  if (unicode) {
    for (size_t j = y1 + 1; j <= y2 - 1; j++) {
      set_char(x1, j, border[0]);
      set_char(x2, j, border[0]);
    }
    for (size_t i = x1 + 1; i <= x2 - 1; i++) {
      set_char(i, y1, border[1]);
      set_char(i, y2, border[1]);
    }
    set_char(x1, y1, border[2]);
    set_char(x2, y1, border[3]);
    set_char(x1, y2, border[4]);
    set_char(x2, y2, border[5]);
  } else {
    for (size_t j = y1 + 1; j <= y2 - 1; j++) {
      set_char(x1, j, '|');
      set_char(x2, j, '|');
    }
    for (size_t i = x1 + 1; i <= x2 - 1; i++) {
      set_char(i, y1, '-');
      set_char(i, y2, '-');
    }
    set_char(x1, y1, '+');
    set_char(x2, y1, '+');
    set_char(x1, y2, '+');
    set_char(x2, y2, '+');
  }
}

void Window::clear() {
  for (size_t j = 1; j <= h; j++) {
    for (size_t i = 1; i <= w; i++) {
      set_char(i, j, ' ');
      set_fg(i, j, fg::reset);
      set_bg(i, j, bg::reset);
      set_style(i, j, style::reset);
    }
  }
}

std::string Window::render() {
  std::string out;
  out.append(cursor_off());
  fg current_fg = fg::reset;
  bg current_bg = bg::reset;
  style current_style = style::reset;
  for (size_t j = 1; j <= h; j++) {
    out.append(move_cursor(y0 + j - 1, x0));
    for (size_t i = 1; i <= w; i++) {
      bool update_fg = false;
      bool update_bg = false;
      bool update_style = false;
      if (current_fg != get_fg(i, j)) {
        current_fg = get_fg(i, j);
        update_fg = true;
      }
      if (current_bg != get_bg(i, j)) {
        current_bg = get_bg(i, j);
        update_bg = true;
      }
      if (current_style != get_style(i, j)) {
        current_style = get_style(i, j);
        update_style = true;
        if (current_style == style::reset) {
          // style::reset resets fg and bg colors too, we have to set them again
          // if they are non-default, but if fg or bg colors are reset, we do
          // not update them, as style::reset already did that.
          update_fg = (current_fg != fg::reset);
          update_bg = (current_bg != bg::reset);
        }
      }
      // Set style first, as style::reset will reset colors too
      if (update_style)
        out.append(color(get_style(i, j)));
      if (update_fg)
        out.append(color(get_fg(i, j)));
      if (update_bg)
        out.append(color(get_bg(i, j)));
      codepoint_to_utf8(out, get_char(i, j));
    }
  }
  if (current_fg != fg::reset)
    out.append(color(fg::reset));
  if (current_bg != bg::reset)
    out.append(color(bg::reset));
  if (current_style != style::reset)
    out.append(color(style::reset));
  out.append(cursor_on());
  return out;
}

} // namespace Term

//---------------------------------------------------------------------------//
// end of ds++/terminal.cc
//---------------------------------------------------------------------------//
