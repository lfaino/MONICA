@echo off
rem MONICA - double-click to start. Windows.
setlocal
if "%MONICA_IMAGE%"=="" set MONICA_IMAGE=ghcr.io/REPLACE-ME/monica:latest
if "%MONICA_PORT%"=="" set MONICA_PORT=8765
if "%MONICA_HOME%"=="" set MONICA_HOME=%USERPROFILE%\MONICA

echo MONICA - avvio in corso

where docker >nul 2>&1
if errorlevel 1 (
    echo.
    echo Docker non e' installato.
    echo Scaricalo da https://www.docker.com/products/docker-desktop/ , installalo,
    echo avvialo e poi riapri questo file.
    pause
    exit /b 1
)
docker info >nul 2>&1
if errorlevel 1 (
    echo Docker e' installato ma non e' in esecuzione: avvia Docker Desktop e riprova.
    pause
    exit /b 1
)

if not exist "%MONICA_HOME%\reads" mkdir "%MONICA_HOME%\reads"
if not exist "%MONICA_HOME%\database" mkdir "%MONICA_HOME%\database"
if not exist "%MONICA_HOME%\results" mkdir "%MONICA_HOME%\results"
echo I tuoi dati stanno in %MONICA_HOME%

docker rm -f monica >nul 2>&1
echo Scarico l'immagine (solo la prima volta, qualche minuto)...
docker pull %MONICA_IMAGE%
if errorlevel 1 ( echo Download fallito. & pause & exit /b 1 )

docker run -d --name monica -p %MONICA_PORT%:8000 ^
    -v "%MONICA_HOME%\reads:/data:ro" ^
    -v "%MONICA_HOME%\database:/db" ^
    -v "%MONICA_HOME%\results:/jobs" ^
    %MONICA_IMAGE%
if errorlevel 1 ( echo Avvio fallito. & pause & exit /b 1 )

timeout /t 8 /nobreak >nul
start "" "http://localhost:%MONICA_PORT%/"
echo.
echo MONICA resta attivo in background. Per fermarlo:  docker stop monica
pause
