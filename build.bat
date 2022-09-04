@echo Compile programm...
@echo off
C:/TDM-GCC-32/bin/g++.exe -g .\main.cpp -o .\main.exe -lbgi -lgdi32 -lcomdlg32 -luuid -loleaut32 -lole32  
if %errorlevel% == 0 (@echo Run application... && @echo off && .\main.exe)