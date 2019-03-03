@echo off
rem ---------------------------------------------------------------------------
rem File  : regression/update_regrssion_dir.bat
rem Date  : Tuesday, May 31, 2016, 14:48 pm
rem Author: Kelly Thompson
rem Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
rem         All rights are reserved.
rem ---------------------------------------------------------------------------

set LOGDIR=c:\regress\logs
set GIT="C:\Users\107638\AppData\Local\Programs\Git\bin\git.exe"
set PATH="C:\Users\107638\AppData\Local\Programs\Git\bin";%PATH%

echo Updating regression system...
echo.
echo LOGDIR = %LOGDIR%
echo GIT    = %GIT%
echo.

echo cd c:\regress\draco
cd /d c:\regress\draco
echo %GIT% pull
%GIT% pull
echo cd c:\regress\jayenne
cd /d c:\regress\jayenne
echo %GIT% pull
%GIT% pull

rem ---------------------------------------------------------------------------
rem End regression/update_regression_dir.bat
rem ---------------------------------------------------------------------------
