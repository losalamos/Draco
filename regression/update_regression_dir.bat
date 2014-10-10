@echo off
set SVN_SSH="c:\\Program Files\\TortoiseSVN\\bin\\TortoisePlink.exe"
set SVN="c:\\Program Files\\TortoiseSVN\\bin\\svn.exe"
set LOGDIR=e:\cdash\logs
cd /d e:\cdash\draco\regression
%SVN% update > %LOGDIR%\update_regression_dir.log

rem 2&>1


