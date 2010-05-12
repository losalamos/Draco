#!/bin/bash2
###############################################################################
## .bash_logout
##
## $Id$
###############################################################################

## Kill any ssh agent

ssh-agent -k

## Remove files at end of session:

rm -rf ~/.DCOPserver* ~/.kxml* ~/.MCOP* ~/.mcop* ~/.xsession-errors
# rm -rf ~/.kde ~/Desktop
# rm -rf ~/.gnome-desktop ~/.*~ ~/.gnome_private ~/.viminfo
rm -rf ~/.autosave ~/.saves-* ~/nsmail ~/.netscape ~/.flexlmrc
# cd ${HOME}; cleanemacs

