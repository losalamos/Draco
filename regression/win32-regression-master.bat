@echo off
rem ---------------------------------------------------------------------------
rem File  : regression/win32-regression-master.bat
rem Date  : Tuesday, May 31, 2016, 14:48 pm
rem Author: Kelly Thompson
rem Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
rem         All rights are reserved.
rem ---------------------------------------------------------------------------

rem Resource for batch file syntax and commands: https://ss64.com/nt/

rem Use Task Scheduler to create a task that runs this script every night at 
rem midnight.

rem Consider using date-based log file names.
rem c:\myscript.%date:~-4%%date:~4,2%%date:~7,2%.%time::=%.log 2>&1

set logdir=c:\regress\logs

rem Try to avoid LNK1109 errors (cannot remove file).
set TEMP=c:\work\temp
set TMP=c:\work\temp

rem Ensure Git/bin/sh.exe is not in the PATH.  It interferes with mingw operations.
rem set PATH=%PATH:c:\Program Files\Git\bin;=%
set PATH=%PATH:C:\Users\107638\AppData\Local\Programs\Git\bin\;=%

rem Fix LANL proxy issues
C:\Users\107638\AppData\Local\Programs\Git\bin\git.exe config --global http.proxy http://proxyout.lanl.gov:8080
rem Don't set these because they prevent posting to rtt.lanl.gov/cdash.
rem set HTTP_PROXY=http://proxyout.lanl.gov:8080
rem set HTTPS_PROXY=http://proxyout.lanl.gov:8080
set HTTP_PROXY=
set HTTPS_PROXY=

rem Change '/c' to '/k' to keep the command window open after these commands 
rem finish.

%comspec% /c ""c:\regress\draco\regression\update_regression_dir.bat"" > %logdir%\update_regression_dir.log 2>&1
%comspec% /c ""c:\regress\draco\regression\win32-regress.bat"" > %logdir%\regression-master.log 2>&1

rem ---------------------------------------------------------------------------
rem end file regression/win32-regression-master.bat
rem ---------------------------------------------------------------------------
