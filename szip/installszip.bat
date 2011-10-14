mkdir szipbin
cd szipbin
mkdir lib
mkdir dll
mkdir include
cd ..
copy all\lib\release\szlib.lib szipbin\lib
copy all\libdll\release\szlibdll.dll szipbin\dll
copy all\libdll\release\szlibdll.lib szipbin\dll
copy src\szlib.h szipbin\include
copy src\szip_adpt.h szipbin\include
copy src\ricehdf.h szipbin\include
copy src\SZconfig.h szipbin\include