@echo off
REM Windows entry point: find Rscript, prepare packages, then launch the local app.
REM Use this file's directory so relative R and data paths resolve consistently.
setlocal
cd /d "%~dp0"
set "LABSTATS_R="
REM Prefer Rscript on PATH; otherwise look for an installed R under Program Files.
for /f "delims=" %%I in ('where Rscript.exe 2^>nul') do if not defined LABSTATS_R set "LABSTATS_R=%%I"
if not defined LABSTATS_R for /d %%D in ("%ProgramFiles%\R\R-*") do if exist "%%~D\bin\Rscript.exe" set "LABSTATS_R=%%~D\bin\Rscript.exe"
if not defined LABSTATS_R (
  echo R was not found. Install R, then try again.
  pause
  exit /b 1
)
REM Stop if package setup fails rather than opening an app with missing dependencies.
"%LABSTATS_R%" setup.R
if errorlevel 1 (
  pause
  exit /b 1
)
"%LABSTATS_R%" launch.R
pause
