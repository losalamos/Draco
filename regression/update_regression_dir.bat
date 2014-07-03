@echo off
set SVN_SSH="c:\\Program Files\\TortoiseSVN\\bin\\TortoisePlink.exe"
set SVN=c:\\Program Files\\TortoiseSVN\\bin\\svn.exe
set LOGDIR=d:\cdash\logs
cd /d d:\cdash\draco\regression
%SVN% update > %LOGDIR%\update_regression_dir.log

rem 2&>1


