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

" fix trailing spaces when you write a file
autocmd BufWritePre * :call RemoveTrailing()

" sets a gray warning line at 81 columns
if exists('+colorcolumn')
  highlight ColorColumn ctermbg=gray
  set colorcolumn=81
else
  au BufWinEnter * let w:m2=matchadd('ErrorMsg', '\%>80v.\+', -1)
endif
