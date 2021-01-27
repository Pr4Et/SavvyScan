;How to compile:
;Examples:
;Win32 with no CUDA:
;"\Program Files\Inno Setup 5\ISCC.exe" /dImodVersion=4.7.10 imod_install.iss
;Win64 with CUDA:
;"\Program Files (x86)\Inno Setup 5\ISCC.exe" /dImodVersion=4.8.47 /dWin64="" /dCuda=2.1 imod_install.iss

[Setup]
AppPublisher=BL3DEMC, University of Colorado
AppPublisherURL=http://bio3d.colorado.edu/imod/
AppSupportURL=http://bio3d.colorado.edu/imod/
UsePreviousAppDir=no
;Turns off user's ability to change install directory
DisableDirPage=yes
AppName=IMOD
#ifdef Win64
#ifdef Cuda
AppVerName=IMOD version {#ImodVersion}{#NoLibs} for win64 with CUDA{#Cuda}
#else
AppVerName=IMOD version {#ImodVersion}{#NoLibs} for win64
#endif
#else
AppVerName=IMOD version {#ImodVersion}{#NoLibs}
#endif
DefaultDirName={code:getAppDir}
ChangesEnvironment=yes
AlwaysShowDirOnReadyPage=yes
CreateUninstallRegKey=no
Uninstallable=no
UpdateUninstallLogAppName=no
PrivilegesRequired=none
;Function calls do not work with OutputBaseFilename because it is used for packaging the installable.
#ifdef Win64
#ifdef Cuda
OutputBaseFilename=imod_{#ImodVersion}{#NoLibs}_win64_CUDA{#Cuda}
#else
OutputBaseFilename=imod_{#ImodVersion}{#NoLibs}_win64
#endif
#else
OutputBaseFilename=imod_{#ImodVersion}{#NoLibs}_win
#endif


[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"


[Files]
;Must be hard coded.
Source: "fixCygPython.sh"; DestDir: "{app}"; Flags: ignoreversion
Source: installIMOD; DestDir: "{app}"; Flags: deleteafterinstall ignoreversion
#ifdef Win64
#ifdef Cuda
Source: "imod_{#ImodVersion}{#NoLibs}_win64_CUDA{#Cuda}.tar.gz"; DestDir: "{app}"; Flags: deleteafterinstall ignoreversion; Afterinstall: installIMOD
#else
Source: "imod_{#ImodVersion}{#NoLibs}_win64.tar.gz"; DestDir: "{app}"; Flags: deleteafterinstall ignoreversion; Afterinstall: installIMOD
#endif
#else
Source: "imod_{#ImodVersion}{#NoLibs}_win.tar.gz"; DestDir: "{app}"; Flags: deleteafterinstall ignoreversion; Afterinstall: installIMOD
#endif


[Registry]
;All installations
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "IMOD_DIR"; ValueData: "{app}\IMOD"; Check: (not isFailed) and IsAdminLoggedOn
Root: HKCU; Subkey: "Environment"; ValueType: string; ValueName: "HOME"; ValueData: "{code:getHome}"; Check: not isFailed



;Cygwin installations
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: none; ValueName: "IMOD_PLUGIN_DIR"; Flags: deletevalue; Check: (not isFailed) and isCygwin and IsAdminLoggedOn
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: none; ValueName: "IMOD_CALIB_DIR"; Flags: deletevalue; Check: (not isFailed) and isCygwin and IsAdminLoggedOn
Root: HKCU; Subkey: "Environment"; ValueType: string; ValueName: "IMOD_DIR"; ValueData: "{app}\IMOD"; Check: (not isFailed) and isCygwin and (not IsAdminLoggedOn)
;  CYGWIN env variable
Root: HKCU; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "CYGWIN"; ValueData: "nodosfilewarning"; Flags: createvalueifdoesntexist; Check: (not isFailed) and isCygwin and IsAdminLoggedOn
Root: HKCU; Subkey: "Environment"; ValueType: string; ValueName: "CYGWIN"; ValueData: "nodosfilewarning"; Flags: createvalueifdoesntexist; Check: (not isFailed) and isCygwin and (not IsAdminLoggedOn)
;Windows-only installations - admin privledges only
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "IMOD_PLUGIN_DIR"; ValueData: "{app}\IMOD\lib\imodplug"; Check: (not isFailed) and (not isCygwin) and IsAdminLoggedOn
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "IMOD_CALIB_DIR"; ValueData: "C:\ProgramData\IMOD"; Check: (not isFailed) and (not isCygwin) and IsAdminLoggedOn and DirExists('C:\ProgramData')
;Cygwin installations - the combination of single user and Windows-only does not work.
;Duplicate line - see if this can be removed
Root: HKCU; Subkey: "Environment"; ValueType: string; ValueName: "IMOD_DIR"; ValueData: "{app}\IMOD"; Check: (not isFailed) and isCygwin and (not IsAdminLoggedOn)


[Run]
Filename: "{cmd}"; Parameters: "/C move {code:getQuotedAppDir}\installIMOD.log {code:getQuotedAppDir}\IMOD"; WorkingDir: "{app}"; StatusMsg: "Moving installIMOD.log to IMOD..."; Check: (not isFailed)
Filename: "{cmd}"; Parameters: "/C move {code:getQuotedAppDir}\fixCygPython.sh {code:getQuotedAppDir}\IMOD"; WorkingDir: "{app}"; StatusMsg: "Moving fixCygPython.sh to IMOD..."; Check: (not isFailed)
Filename: "{cmd}"; Parameters: "/C move {code:getQuotedAppDir}\fixCygPython.log {code:getQuotedAppDir}\IMOD"; WorkingDir: "{app}"; StatusMsg: "Moving fixCygPython.log to IMOD..."; Check: (not isFailed)


[Code]
var
  Failed: Boolean;
  Cygwin: Boolean;
  Psutil: Boolean;
  BadCygwin: Boolean;
  AppDir: String;
  CygwinDir: String;
  PythonDir: String;
  PythonDirInSysPath: Boolean;
  
function setupCygwin(const rootKey: Integer; const cygwinKeyName: String; const pathName: String): Boolean;
//If a installation of Cygwin that contains Python is found, set global variables and return true.
//Runs fixCygPython if the "python" link exists
var
  path: String;
  bin: String;
begin
  Result := False;
  if RegQueryStringValue(rootKey, cygwinKeyName, pathName, path) and (path <> '') then begin
    path := trim(path);
    if DirExists(path) then begin
      bin := AddBackslash(path) + 'bin';
      if DirExists(bin) then begin
        if fileExists(bin + '\python.exe') or fileExists(bin + '\python') then begin
          AppDir := AddBackslash(path) + 'usr\local';
          CygwinDir := path;
          PythonDir := bin;
          Cygwin := true;
          Result := True;
          BadCygwin := False;
        end else begin
          BadCygwin := True;
        end;
      end;
    end;
  end;
end;

function pythonDirRegistrySearch(const rootKey: Integer; const subKeyName: String): Boolean;
//If a compatible version of Python is found, set the global variables and return true.
//Checks for psutil
var
  versionArray: TArrayOfString;
  version: String;
  index: Integer;
  path: String;
begin
  Result := False;
  if RegGetSubkeyNames(rootKey, subKeyName, versionArray) then begin
    for index := 0 to GetArrayLength(versionArray)-1 do begin
      version := versionArray[index];
      //Take either the first valid python with psutil or, if none of them have psutil,
      //will take the last valid one.
      if RegQueryStringValue(rootKey, AddBackslash(subKeyName) + version + '\InstallPath',
                       '', path) then begin
        if DirExists(path) and fileExists(AddBackslash(path) + 'python.exe') then begin
          AppDir := 'C:\Program Files'
          CygwinDir := '';
          PythonDir := path;
          Cygwin := False;
          Result := True;
          if DirExists(AddBackslash(path) + 'Lib\site-packages\psutil') then begin
            Psutil := True;
            break;
          end;
        end;
      end;
    end;
  end;
end;

function pythonUserQuery(const fallbackPythonDir: String): Boolean;
//Looks for python in a directory.  Returns true if the user chooses a python installation.
var
  path: String;
  folderChosen: Boolean;
begin
  Result := False;
  folderChosen := False;
  path := ExpandConstant('{sd}');
  if length(fallbackPythonDir) = 0 then begin
    folderChosen := BrowseForFolder('A Python installation has not been found.  Choose the Python installation you would like to use with IMOD', path, False);
  end else begin
    folderChosen := BrowseForFolder('A Python installation with psutil has not been found.  Choose the Python installation you would like to use with IMOD, or press Cancel to use ' + fallbackPythonDir, path, False);
  end;
  if folderChosen and DirExists(path) and fileExists(AddBackslash(path) + 'python.exe') then begin
    AppDir := 'C:\Program Files'
    CygwinDir := '';
    PythonDir := path;
    Cygwin := False;
    Result := True;
    if DirExists(AddBackslash(path) + 'Lib\site-packages\psutil') then begin
      Psutil := True;
    end;
  end;
end;

function pythonDirectorySearch(const dir: String; var fallbackPythonDir: String): Boolean;
//Looks for python in a directory.  Sets fallbackPythonDir to the last Python found, if
//it wasn't already set.  Returns true if a python with psutil is found
var
  findRec: TFindRec;
  path: String;
  setFallback: Boolean;
begin
  Result := False;
  dir := AddBackslash(dir);
  if FindFirst((dir + '*'), findRec) then begin
    try
      //The last python is likely to be the most recent version.
      setFallback := length(fallbackPythonDir) = 0;
      repeat
        if findRec.Attributes and FILE_ATTRIBUTE_DIRECTORY = 1 then begin
          path := dir + findRec.Name;
          if fileExists(AddBackslash(path) + 'python.exe') then begin
            //The last python is the fallback
            if setFallback then begin
              fallbackPythonDir := PythonDir;
            end;
            AppDir := 'C:\Program Files'
            CygwinDir := '';
            PythonDir := path;
            Cygwin := False;
            if DirExists(AddBackslash(path) + 'Lib\site-packages\psutil') then begin
              Psutil := True;
              Result := True;
              break;
            end;
          end;
        end;
      until not FindNext(findRec);
    finally
      FindClose(findRec);
    end;
  end;
end;

function getSysPath(var hkey: Integer; var subkey: String): String;
var
  sysPath: String;
begin
  if not IsAdminLoggedOn() then begin
    hkey := HKEY_CURRENT_USER;
    subkey := 'Environment';
  end else begin
    hkey := HKEY_LOCAL_MACHINE;
    subkey := 'SYSTEM\CurrentControlSet\Control\Session Manager\Environment';
  end;
  RegQueryStringValue(hkey, subkey, 'Path', sysPath);
  Result := Trim(sysPath);
end;

function pythonSysPathSearch(var fallbackPythonDir: String): Boolean;
//Looks for python PATH.  Sets fallbackPythonDir to the first Python found, if it wasn't
//already set.  Returns true if a python with psutil is found
var
  sysPath: String;
  path: String;
  position: Integer;
  index: Integer;
  hkey: Integer;
  subkey: String;
begin
  Result := False;
  sysPath := getSysPath(hkey, subkey);
  index := 1;
  if (not (Length(sysPath) = 0)) and (not (Pos('Python', sysPath) = 0)) then begin
    repeat
      position := Pos(';', sysPath);
      if not (position = 0) then begin
        path := copy(sysPath, index, position - index)
        if (not (length(path) = 0)) and (not (Pos('Python', sysPath) = 0)) and dirExists(path) and fileExists(AddBackslash(path) + 'python.exe') then begin
          if length(fallbackPythonDir) = 0 then begin
            fallbackPythonDir := PythonDir;
          end;
          PythonDirInSysPath := True;
          AppDir := 'C:\Program Files'
          CygwinDir := '';
          PythonDir := path;
          Cygwin := False;
          if DirExists(AddBackslash(path) + 'Lib\site-packages\psutil') then begin
            Psutil := True;
            Result := True;
          end;          
        end;
        index := position + 1;
        sysPath := copy(sysPath, index, length(sysPath));
      end;
    until position = 0;
  end;
end;

function pythonRegistrySearch(const rootKey: Integer; var fallbackPythonDir: String): Boolean;
//Looks for all possible places under this rootKey. Sets fallbackPythonDir to the first
//Python found, if it wasn't already set.  Returns true if a python with psutil is found
var
  path: String;
begin
  Result := False;
  Psutil := False;
  //App Path entry refers to the most recently installed python.  It is not always there
  //because it is deleted when any Python is uninstalled.
  if RegQueryStringValue(rootKey, 'Software\Microsoft\Windows\CurrentVersion\App Paths\Python.exe', '', path) and fileExists(path) then begin
    if length(fallbackPythonDir) = 0 then begin
      fallbackPythonDir := PythonDir;
    end;
    path := ExtractFileDir(path);
    AppDir := 'C:\Program Files'
    CygwinDir := '';
    PythonDir := path;
    Cygwin := False;
    if DirExists(AddBackslash(path) + 'Lib\site-packages\psutil') then begin
      Result := True;
      Psutil := True;
    end;
  end;
  if (not Psutil) and pythonDirRegistrySearch(rootKey, 'Software\Python\PythonCore') then begin
    if length(fallbackPythonDir) = 0 then begin
      fallbackPythonDir := PythonDir;
    end;
    if Psutil then begin
      Result := True;
    end;
  end;
  if (not Psutil) and pythonDirRegistrySearch(rootKey, 'Software\Wow6432Node\Python\PythonCore') then begin
    //On Windows 7, this is where 32 bit software is registered on a 64 bit OS.
    if length(fallbackPythonDir) = 0 then begin
      fallbackPythonDir := PythonDir;
    end;
    if Psutil then begin
      Result := True;
    end;
  end;
  //FallbackPythonDir should contain what is hopefully the best version of python
  //available.  All the other variables set when python is found contain constants.
  if Result and (not Psutil) then begin
    PythonDir := fallbackPythonDir;
  end;
end;

function setupWindows(const win32Comp: Boolean): Boolean;
//Looks for a python installation with psutil.  Fallback is a python without psutil
//Returns true if a python is found
var
  fallbackPythonDir: String;
begin
  Result := True;
  fallbackPythonDir := '';
  if not pythonRegistrySearch(HKEY_LOCAL_MACHINE_32, fallbackPythonDir) then begin
    if win32Comp or (not pythonRegistrySearch(HKEY_LOCAL_MACHINE_64, fallbackPythonDir)) then begin
      if not pythonRegistrySearch(HKEY_CURRENT_USER_32, fallbackPythonDir) then begin
        if win32Comp or (not pythonRegistrySearch(HKEY_CURRENT_USER_64, fallbackPythonDir)) then begin
          if not pythonSysPathSearch(fallbackPythonDir) then begin
            if not pythonDirectorySearch(ExpandConstant('{sd}'), fallbackPythonDir) then begin
              if not pythonDirectorySearch(ExpandConstant('{pf}'), fallbackPythonDir) then begin
                if win32Comp or not pythonDirectorySearch(ExpandConstant('{pf32}'), fallbackPythonDir) then begin
                  if not pythonUserQuery(fallbackPythonDir) then begin
                    //Either python was not found, or psutil is missing from all pythons
                    if length(fallbackPythonDir) = 0 then begin
                      //No python was found
                      MsgBox('A compatible Python version has not been installed on this computer.  Unable to install.', mbCriticalError, MB_OK);
                      Result := False;
                    end else begin
                      //Python was found, but not psutil
                      PythonDir := fallbackPythonDir;
                    end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
  end;
end;

function InitializeSetup(): Boolean;
//Initializes the global variables.
//Check 64-bit status, checks whether the right packages are installed, and sets the global variables.
//Return false if 64-bit status is wrong, or the right packages are not installed.
//If false is returned, the installer terminates.
var
  Win32Comp : Boolean;
begin
  Failed := False;
  BadCygwin := False;
  Result := True;
  win32Comp := (not isWin64());
#ifdef Win64
  //Make sure this computer is a 64-bit computer.
  if (not isWin64()) then begin
    MsgBox('Not a 64-bit computer.  Unable to install', mbCriticalError, MB_OK);
    Result := False;
  end else begin
#endif
    //find out if Cygwin is installed
    if win32Comp or (not setupCygwin(HKEY_LOCAL_MACHINE_64, 'SOFTWARE\Cygwin\setup', 'rootdir')) then begin
      if not setupCygwin(HKEY_LOCAL_MACHINE_32, 'SOFTWARE\Cygwin\setup', 'rootdir') then begin
        if win32Comp or (not setupCygwin(HKEY_CURRENT_USER_64, 'SOFTWARE\Cygwin\setup', 'rootdir')) then begin  
          if not setupCygwin(HKEY_CURRENT_USER_32, 'SOFTWARE\Cygwin\setup', 'rootdir') then begin
            if not setupCygwin(HKEY_LOCAL_MACHINE_32, 'SOFTWARE\Cygnus Solutions\Cygwin\mounts v2\/', 'native') then begin                 
              if win32Comp or not (setupCygwin(HKEY_LOCAL_MACHINE_64, 'SOFTWARE\Cygnus Solutions\Cygwin\mounts v2\/', 'native')) then begin
                if not setupCygwin(HKEY_CURRENT_USER_32, 'SOFTWARE\Cygnus Solutions\Cygwin\mounts v2\/', 'native') then begin
                  if win32Comp or (not setupCygwin(HKEY_CURRENT_USER_64, 'SOFTWARE\Cygnus Solutions\Cygwin\mounts v2\/', 'native')) then begin
                    //find out if Windows Python is installed.  A single-user Windows-only install will not work because a non-admin
                    //account cannot write to Program Files.
                    if not setupWindows(win32Comp) then begin
                      Result := False;
                    end;                               
                  end;              
                end;                    
              end; 
            end;
          end;
        end;
      end;
    end;
#ifdef Win64
  end;
#endif
  //If Cygwin is installed but does not not have Python installed, and there is a valid Windows Python installed,
  //ask the user before going with the Windows-only installation.
  if BadCygwin and (AppDir <> '') and (MsgBox('There is no Python installed in Cygwin.  ' +
                    'If you proceed with this'#13#10'installation, ' +
                    'some IMOD programs will not be runnable from Cygwin terminals.'#13#10#13#10
                    'To use IMOD with a Cygwin terminal, ' +
                    'exit this installer and use the Cygwin installer to install Python.'#13#10#13#10
                    'Exit?', mbConfirmation, MB_YESNO) = IDYES) then begin
    MsgBox('IMOD had not been installed.  Exiting installer.', mbInformation, MB_OK);
    Result := False;
  end;
end;

function getAppDir(const dummy: String): String;
begin
  Result := AppDir;
end;

function getQuotedAppDir(const dummy: String): String;
begin
  Result := '"' + AppDir + '"';
end;

function isFailed(): Boolean;
begin
  Result := Failed;
end;

function isCygwin(): Boolean;
begin
  Result := Cygwin;
end;

function getPythonDir(const dummy: String): String;
begin
  Result := PythonDir;
end;

function getCygwinDir(const dummy: String): String;
begin
  Result := CygwinDir;
end;

procedure installIMOD();
//Run installIMOD.
var
  command: String;
  fixCommand: String;
  outputBaseFilename: String;
  skip: String;
  returnCode: Integer;
  //Runs the installIMOD script.  Sets Failed if installIMOD failed.
begin
  if Cygwin then begin
    skip := '';
  end else begin
    skip := '-skip ';
  end;
#ifdef Win64
#ifdef Cuda
  outputBaseFilename := 'imod_{#ImodVersion}{#NoLibs}_win64_CUDA{#Cuda}'
#else
  outputBaseFilename := 'imod_{#ImodVersion}{#NoLibs}_win64'
#endif
#else
  outputBaseFilename := 'imod_{#ImodVersion}{#NoLibs}_win'
#endif
  if Cygwin and fileExists(PythonDir + '\python') then begin
    fixCommand := '/C '+ PythonDir + '\bash.exe ' + ExpandConstant('{app}')+'\fixCygPython.sh > fixCygPython.log 2>&1';
    Exec(ExpandConstant('{cmd}'), fixCommand, ExpandConstant('{app}'), SW_SHOW, ewWaitUntilTerminated, returnCode);
    //Exec(PythonDir + '\bash.exe',ExpandConstant('{app}')+'\fixCygPython.sh',PythonDir,SW_SHOW,ewWaitUntilTerminated,returnCode)
    if returnCode <> 0 then begin
      MsgBox('A script to create a python.exe in cygwin has failed, returnCode=' + IntToStr(returnCode) + '.  Error messages are in fixCygPython.log in ' + ExpandConstant('{app}') + 
      ' If IMOD install fails, try running the command "cp -L /bin/python /bin/python.exe" in a Cygwin terminal', mbCriticalError, MB_OK);
    end;
  end;
  command := '/C PATH=' + PythonDir + ';%PATH% && echo Installing IMOD.......... && python installIMOD -yes ' + skip + outputBaseFilename + '.tar.gz > installIMOD.log 2>&1';
  if Exec(ExpandConstant('{cmd}'), command, ExpandConstant('{app}'), SW_SHOW, ewWaitUntilTerminated, returnCode) then begin
    if returnCode <> 0 then begin
      Failed := True;
      MsgBox('IMOD Install has failed!.  IMOD will not be installed (ignore messages to the contrary from this ' +
             'installer).  Error messages are in installIMOD.log in ' + ExpandConstant('{app}') +
             '.  Please see http://bio3d.colorado.edu/imod for alternatives.  Working directory=' +
             ExpandConstant('{app}') + ',command=' + command + ',returnCode=' + IntToStr(returnCode), mbCriticalError, MB_OK);
    end;
  end else begin
    Failed := True;
    MsgBox('IMOD Install has failed!.  IMOD will not be installed (ignore messages to the contrary from this ' +
           'installer).  Please see http://bio3d.colorado.edu/imod for alternatives.  Working directory=' +
           ExpandConstant('{app}') + ',command=' + command + ',returnCode=' + IntToStr(returnCode), mbCriticalError, MB_OK);
  end;
end;

function getHome(const dummy: String): String;
//Return the value for the HOME variable.
var
  version: TWindowsVersion;
begin
  if Cygwin then begin
    //Cygwin
    Result := ExpandConstant(CygwinDir + '\home\{username}');
  end else begin
    GetWindowsVersionEx(version);
    if version.Major >= 6 then begin
      //Vista and on
      Result := ExpandConstant('C:\Users\{username}');
    end else begin
      //XP and earlier
      Result := ExpandConstant('C:\Documents and Settings\{username}');
    end;
  end;
end;

procedure addToSysPath(const path: String; var sysPath: String);
//Add a path to the PATH environment variable value.
var
  addPath: Boolean;
  index: Integer;
begin
  addPath := True;
  //Don't add the path if it already exists
  index := Pos(path, sysPath);
  if index <> 0 then begin
    index := index + Length(path) + 1;
    if (index <> Length(sysPath)) and (sysPath[index] <> ';') then begin
      addPath := False;
    end;
  end;
  if addPath then begin
    sysPath := path + ';' + sysPath;
  end;
end;

procedure CurStepChanged(const CurStep: TSetupStep);
var
  sysPath: String;
  hkey: Integer;
  subkey: String;
begin
  if CurStep = ssPostInstall then begin
    //Modify the PATH environment variable.
    if not Failed then begin
      sysPath := getSysPath(hkey, subkey);
      addToSysPath(AppDir + '\IMOD\bin', sysPath);
      if Cygwin then begin
        addToSysPath(CygwinDir + '\bin', sysPath);
      end else if not PythonDirInSysPath then begin
        addToSysPath(PythonDir, sysPath);
      end;
      RegWriteStringValue(hkey, subkey, 'Path', sysPath);
      if (not Cygwin) and (not Psutil) then begin
        //warn about missing psutil
        MsgBox('WARNING: the module psutil does not appear to be installed in ' + PythonDir,
               mbInformation, MB_OK);
      end;
    end;
  end;
end;
