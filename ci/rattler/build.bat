SET PYTHONUTF8=1
SET PYTHONIOENCODING=utf-8

maturin build --release
if %ERRORLEVEL% neq 0 exit %ERRORLEVEL%
%PYTHON% -m pip install --find-links=target\wheels %PKG_NAME%
if %ERRORLEVEL% neq 0 exit %ERRORLEVEL%
