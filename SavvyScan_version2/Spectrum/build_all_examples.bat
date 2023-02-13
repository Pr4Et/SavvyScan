@echo off

rem This script checks all subdirectories for .vcproj files and tries to build
rem each project as 32bit and 64bit application, both debug and release.

setlocal EnableDelayedExpansion

set LIST=

rem Konvertierung vcproj -> vcxproj:
rem devenv /Upgrade xyz.vcproj
rem danach gibt es Linkerfehler, weil die Libs falsch eingebunden sind.
rem dagegen hilft sed:
rem sed -i "s/CustomBuild/Library/" %%P
rem danach sind noch die Ausgabepfade Murks
rem wieder sed:
rem sed -i "/<Link>/{n;s/.*/      <OutputFile>$(TargetPath)<\/OutputFile>/ }" %%P

rem ältere msbuild Versionen geben komische Fehler
set MSBUILD="C:\Program Files (x86)\MSBuild\14.0\Bin\MSBuild.exe"

rem loop over all directories
FOR /F %%D in ('dir /A:D /B') DO (
    pushd %%D

    rem loop over all .vcproj files in directory
    FOR /F %%P in ('dir /B *.vcxproj 2^>nul') DO (
        rem vcbuild /r /platform:Win32 %%P Debug
        rem 
        %MSBUILD% /target:Rebuild /property:Configuration=Debug;Platform=x86 %%P
        if ERRORLEVEL 1 (
            echo.
            echo ERROR: %%D/%%P Win32 Debug FAILED
            set LIST=!LIST! "%%D/%%P Win32 Debug" 
        )
        rem vcbuild /r /platform:Win32 %%P Release
        rem sed -i "s/CustomBuild/Library/" %%P
        rem sed -i "/<Link>/{n;s/.*/      <OutputFile>$(TargetPath)<\/OutputFile>/ }" %%P
        %MSBUILD% /target:Rebuild /property:Configuration=Release;Platform=x86 %%P
        if ERRORLEVEL 1 (
            echo.
            echo ERROR: %%D/%%P Win32 Release FAILED
            set LIST=!LIST! "%%D/%%P Win32 Release" 
        )
        rem vcbuild /r /platform:x64   %%P Debug
        rem sed -i "s/CustomBuild/Library/" %%P
        rem sed -i "/<Link>/{n;s/.*/      <OutputFile>$(TargetPath)<\/OutputFile>/ }" %%P
        %MSBUILD% /target:Rebuild /property:Configuration=Debug;Platform=x64 %%P
        if ERRORLEVEL 1 (
            echo.
            echo ERROR: %%D/%%P x64 Debug FAILED
            set LIST=!LIST! "%%D/%%P x64 Debug" 
        )
        rem vcbuild /r /platform:x64   %%P Release
        rem sed -i "s/CustomBuild/Library/" %%P
        rem sed -i "/<Link>/{n;s/.*/      <OutputFile>$(TargetPath)<\/OutputFile>/ }" %%P
        %MSBUILD% /target:Rebuild /property:Configuration=Release;Platform=x64 %%P
        if ERRORLEVEL 1 (
            echo.
            echo ERROR: %%D/%%P x64 Release FAILED
            set LIST=!LIST! "%%D/%%P x64 Release" 
        )
    )

    popd
)

echo.
echo *****************************************
IF DEFINED LIST (
    echo The following examples failed to compile:
    echo.
    FOR %%a in (%LIST%) DO echo %%a
    echo.
    exit /B 1
) ELSE (
    echo All examples compiled successfully
    echo.
    exit /B 0
)

