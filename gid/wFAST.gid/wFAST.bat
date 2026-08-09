rem OutputFile: wfastinp.log
rem ErrorFile: wfast_gid.err

rem --------------------------------------------------
rem Clean previous results
rem --------------------------------------------------

del "%1.post.res"
del "%1.post.msh"

del "gid.dat"
del "wfastinp.log"
del "wfast_gid.err"
del "wfast_windows.err"
del "wfastinp.dat"
del "wfastinp.res"

del "wfast-po_solid.dat"
del "wfast-po_joint.dat"
del "wfastinp_point_load.dat"
del "wfastinp_surface_load.dat"
del "wfastinp_po_solid.dat"
del "wfastinp_po_joint.dat"
del "wfastinp_temperature.dat"
del "wfastinp_drho.dat"

rem ------------------------------------------------------------
rem Check input and solver
rem ------------------------------------------------------------

if not exist "%1.dat" goto no_input
if not exist "%3\exec\wFAST.exe" goto no_solver
if not exist "%3\exec\wFAST.cfg" goto no_config

rem ------------------------------------------------------------
rem Prepare wFAST input
rem ------------------------------------------------------------

rename "%1.dat" "gid.dat"

if not exist "gid.dat" goto input_error

rem ------------------------------------------------------------
rem Run wFAST
rem ------------------------------------------------------------

"%3\exec\wFAST.exe" "%3\exec\wFAST.cfg"

if errorlevel 1 goto solver_error

rem ------------------------------------------------------------
rem Check results
rem ------------------------------------------------------------

if not exist "wfastinp.post.res" goto no_result
if not exist "wfastinp.post.msh" goto no_result

rename "wfastinp.post.res" "%1.post.res"
rename "wfastinp.post.msh" "%1.post.msh"

goto end

:input_error
echo ERROR: Failed to prepare gid.dat. > "wfast_gid.err"
goto end

:no_input
echo ERROR: GiD input file was not generated. > "wfast_gid.err"
goto end

:no_solver
echo ERROR: wFAST.exe was not found. Build or install the wFAST solver first. > "wfast_gid.err"
goto end

:no_config
echo ERROR: wFAST.cfg was not found. > "wfast_gid.err"
goto end

:solver_error
echo ERROR: wFAST solver returned a non-zero exit code. >> "wfast_gid.err"
goto end

:no_result
echo ERROR: wFAST did not create the expected GiD post-processing files. >> "wfast_gid.err"
goto end

:end