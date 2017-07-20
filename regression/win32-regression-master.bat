@echo off
rem ---------------------------------------------------------------------------
rem File  : regression/win32-regression-master.bat
rem Date  : Tuesday, May 31, 2016, 14:48 pm
rem Author: Kelly Thompson
rem Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
rem         All rights are reserved.
rem ---------------------------------------------------------------------------

rem Use Task Scheduler to create a task that runs this script every night at 
rem midnight.

rem Consider using date-based log file names.
rem c:\myscript.%date:~-4%%date:~4,2%%date:~7,2%.%time::=%.log 2>&1

set logdir=e:\regress\logs

rem Change '/c' to '/k' to keep the command window open after these commands 
rem finish.

%comspec% /c ""e:\regress\draco\regression\update_regression_dir.bat"" > %logdir%\update_regression_dir.log 2>&1
%comspec% /c ""e:\regress\draco\regression\win32-regress.bat"" > %logdir%\regression-master.log 2>&1

rem ---------------------------------------------------------------------------
rem end file regression/win32-regression-master.bat
rem ---------------------------------------------------------------------------