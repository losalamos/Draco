" vim settings to comply with Draco/Jayenne style guide
" To include the settings and commands in this file, copy it
" to your home directory or add this command to your ~/.vimrc file:
" "source <path_to_draco>/environment/vim/.vimrc"

" Use spaces instead of tabs when indenting with << or >>
set expandtab

" Be smart when using tabs
set smarttab

" 1 tab == 4 spaces
set shiftwidth=4
set tabstop=4

" function that kills trailing whitespaces
fun! RemoveTrailing()
    " save the cursor position
    let l:save_cursor = getpos('.')
    " find replace all trailing spaces to end of line
    %s/\s\+$//e
    " reset cursor position
    call setpos('.', l:save_cursor)
endfun

" Remove trailing spaces when you write a file
autocmd BufWritePre * :call RemoveTrailing()

" ----------------------------------------------------------------------
" sets a gray warning line at 81 columns
if exists('+colorcolumn')
  highlight ColorColumn ctermbg=gray
  set colorcolumn=81
else
  au BufWinEnter * let w:m2=matchadd('ErrorMsg', '\%>80v.\+', -1)
endif

" --------------------------------------------------------------------
" fprettify
" https://github.com/pseewald/fprettify
" Use: 'gq'
" autocmd Filetype fortran setlocal formatprg=fprettify\ --silent

" ----------------------------------------------------------------------
" clang-format integrations
" https://clang.llvm.org/docs/ClangFormat.html#vim-integration
" Requires:
"  - Python 3 or later.
"  - clang-format.py, which is located at $LLVM_ROOT/share/clang/clang-format.py
"  - export PYTHONPATH=${PYTHONPATH:+${PYTHONPATH}:}$LLVM_ROOT/share/clang
" Use: 'C-I' to format current line

map <C-I> :pyf ${LLVM_ROOT}/share/clang/clang-format.py<cr>
imap <C-I> <c-o>:pyf ${LLVM_ROOT}/share/clang/clang-format.py<cr>

fun! CFonsave()
  let l:lines="all"
  pyf ${LLVM_ROOT}/share/clang/clang-format.py
endfun

" Run clang-format when you write a file.
autocmd BufWritePre *.h,*.hh,*.cc,*.cpp :call CFonsave()
