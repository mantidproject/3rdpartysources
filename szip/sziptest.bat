@ECHO OFF

REM This batch file is used to test whether SZIP library is built correctly.
REM By Xuan Bai
REM Created:       9/02/2004
REM Last Modfied:  9/02/2004


echo *****************************************************************************
echo                         SZIP Library Test -- %1 version
echo *****************************************************************************

echo.
echo =============================================================================
echo                         gentest -- static
echo =============================================================================

cd all\gentest\%1

gentest

echo.
echo =============================================================================
echo                         gentest -- dll
echo =============================================================================

cd ..\..\gentestdll\%1

copy ..\..\libdll\%1\szlibdll.dll %SystemRoot%\system >temp.txt

del temp.txt

gentestdll

echo.
echo =============================================================================
echo                         burst -- static
echo =============================================================================

cd ..\..\burst_szip\%1

burst_szip 8 8 10000 16 image.in

copy image.in ..\..\example\%1\image.8.in > temp.txt

del temp.txt

echo.
echo =============================================================================
echo                         example -- static
echo =============================================================================

cd ..\..\example\%1

example

cd ..\..\..






