# run.ps1 - Run MML_VisualizationApp
# Usage:
#   .\run.ps1                         # default demo, Release build
#   .\run.ps1 forms_3d                # typed forms single-vector scene
#   .\run.ps1 forms_field_3d          # vortex tube field scene
#   .\run.ps1 forms_em_3d             # EM dipole field scene
#   .\run.ps1 real_function           # y = f(x) example
#   .\run.ps1 field_3d                # 3D vector field
#   .\run.ps1 all                     # run all visualizations
#   .\run.ps1 forms_3d Debug          # use Debug build

param(
    [string]$Scene  = "",
    [string]$Config = "Release"
)

$root = Resolve-Path "$PSScriptRoot\..\.."
$exe  = "$root\build\src\visualization_examples\$Config\MML_VisualizationApp.exe"

if (-not (Test-Path $exe)) {
    Write-Host "Executable not found: $exe" -ForegroundColor Red
    Write-Host "Run .\build.ps1 first (or .\build.ps1 $Config for a specific config)." -ForegroundColor Yellow
    exit 1
}

if ($Scene -eq "") {
    Write-Host "Running MML_VisualizationApp (default demo, $Config) ..." -ForegroundColor Cyan
    & $exe
} else {
    Write-Host "Running MML_VisualizationApp $Scene ($Config) ..." -ForegroundColor Cyan
    & $exe $Scene
}

exit $LASTEXITCODE
