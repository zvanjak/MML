# build.ps1 - Build MML_VisualizationApp
# Usage:
#   .\build.ps1            # Release build (default)
#   .\build.ps1 Debug      # Debug build

param(
    [string]$Config = "Release"
)

$root = Resolve-Path "$PSScriptRoot\..\.."

Write-Host "Building MML_VisualizationApp ($Config) ..." -ForegroundColor Cyan

cmake --build "$root\build" `
      --config $Config `
      --target MML_VisualizationApp `
      --parallel

if ($LASTEXITCODE -ne 0) {
    Write-Host "Build FAILED (exit code $LASTEXITCODE)" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host "Build succeeded." -ForegroundColor Green
Write-Host "  Executable: $root\build\src\visualization_examples\$Config\MML_VisualizationApp.exe"
