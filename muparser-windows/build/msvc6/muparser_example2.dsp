# Microsoft Developer Studio Project File - Name="muparser_example2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=example2 - Win32 Release Static
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "muparser_example2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "muparser_example2.mak" CFG="example2 - Win32 Release Static"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "example2 - Win32 Debug DLL" (based on "Win32 (x86) Console Application")
!MESSAGE "example2 - Win32 Debug Static" (based on "Win32 (x86) Console Application")
!MESSAGE "example2 - Win32 Release DLL" (based on "Win32 (x86) Console Application")
!MESSAGE "example2 - Win32 Release Static" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "example2 - Win32 Debug DLL"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "..\samples\example2"
# PROP BASE Intermediate_Dir "vc_shared_dbg\example2"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\samples\example2"
# PROP Intermediate_Dir "vc_shared_dbg\example2"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /FD /MDd /GR /EHsc /Od /W4 /I "..\include" /Zi /Gm /GZ /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "_DEBUG" /c
# ADD CPP /nologo /FD /MDd /GR /EHsc /Od /W4 /I "..\include" /Zi /Gm /GZ /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "_DEBUG" /c
# ADD BASE RSC /l 0x409 /d "_CONSOLE" /d "_DEBUG" /i ..\include
# ADD RSC /l 0x409 /d "_CONSOLE" /d "_DEBUG" /i ..\include
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib" /debug
# ADD LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib" /debug

!ELSEIF  "$(CFG)" == "example2 - Win32 Debug Static"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "..\samples\example2"
# PROP BASE Intermediate_Dir "vc_static_dbg\example2"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\samples\example2"
# PROP Intermediate_Dir "vc_static_dbg\example2"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /FD /MDd /GR /EHsc /Od /W4 /I "..\include" /Zi /Gm /GZ /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "_DEBUG" /c
# ADD CPP /nologo /FD /MDd /GR /EHsc /Od /W4 /I "..\include" /Zi /Gm /GZ /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "_DEBUG" /c
# ADD BASE RSC /l 0x409 /d "_CONSOLE" /d "_DEBUG" /i ..\include
# ADD RSC /l 0x409 /d "_CONSOLE" /d "_DEBUG" /i ..\include
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib" /debug
# ADD LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib" /debug

!ELSEIF  "$(CFG)" == "example2 - Win32 Release DLL"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "..\samples\example2"
# PROP BASE Intermediate_Dir "vc_shared_rel\example2"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\samples\example2"
# PROP Intermediate_Dir "vc_shared_rel\example2"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /FD /MD /GR /EHsc /O2 /w /I "..\include" /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "NDEBUG" /c
# ADD CPP /nologo /FD /MD /GR /EHsc /O2 /w /I "..\include" /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "NDEBUG" /c
# ADD BASE RSC /l 0x409 /d "_CONSOLE" /d "NDEBUG" /i ..\include
# ADD RSC /l 0x409 /d "_CONSOLE" /d "NDEBUG" /i ..\include
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib"
# ADD LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib"

!ELSEIF  "$(CFG)" == "example2 - Win32 Release Static"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "..\samples\example2"
# PROP BASE Intermediate_Dir "vc_static_rel\example2"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\samples\example2"
# PROP Intermediate_Dir "vc_static_rel\example2"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /FD /MD /GR /EHsc /O2 /w /I "..\include" /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "NDEBUG" /c
# ADD CPP /nologo /FD /MD /GR /EHsc /O2 /w /I "..\include" /Fd..\samples\example2\example2.pdb /D "WIN32" /D "_CONSOLE" /D "NDEBUG" /c
# ADD BASE RSC /l 0x409 /d "_CONSOLE" /d "NDEBUG" /i ..\include
# ADD RSC /l 0x409 /d "_CONSOLE" /d "NDEBUG" /i ..\include
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib"
# ADD LINK32 ..\lib\muparser.lib /nologo /machine:i386 /out:"..\samples\example2\example2.exe" /subsystem:console /libpath:"..\lib"

!ENDIF

# Begin Target

# Name "example2 - Win32 Debug DLL"
# Name "example2 - Win32 Debug Static"
# Name "example2 - Win32 Release DLL"
# Name "example2 - Win32 Release Static"
# Begin Group "Source Files"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\samples\example2\AssemblyInfo.cpp
# End Source File
# Begin Source File

SOURCE=.\samples\example2\WndError.cpp
# End Source File
# Begin Source File

SOURCE=.\samples\example2\WndMain.cpp
# End Source File
# Begin Source File

SOURCE=.\samples\example2\stdafx.cpp
# End Source File
# End Group
# End Target
# End Project

