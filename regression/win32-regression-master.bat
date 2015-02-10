@echo off
%comspec% /c ""e:\cdash\draco\regression\win32-regress.bat"" x86 > e:\cdash\logs\regression-master.log

rem c:\myscript.%date:~-4%%date:~4,2%%date:~7,2%.%time::=%.log 2>&1
