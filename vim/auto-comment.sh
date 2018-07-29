To turn off auto-insert of comments, you can add these lines to your vimrc:

if you don't have vimrc, just touch vimrc and then input the following:

augroup auto_comment
    au!
    au FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o
augroup END
