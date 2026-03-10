# PowerShell Script to visualize Particle Visualizer 2D data files
# Using Qt MML_ParticleVisualizer2D

# Configuration
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$VisualizerPath = Join-Path $ScriptDir "MML_ParticleVisualizer2D\MML_ParticleVisualizer2D.exe"
$DataDir = Join-Path $ScriptDir "..\..\..\data\ParticleVisualizer2D\Basic"

# Color settings
$Host.UI.RawUI.ForegroundColor = "Cyan"

Write-Host "========================================" -ForegroundColor Green
Write-Host "Particle Visualizer 2D Runner (Qt)" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Check if visualizer exists
if (-not (Test-Path $VisualizerPath)) {
    Write-Host "ERROR: Visualizer not found at: $VisualizerPath" -ForegroundColor Red
    Write-Host "Please ensure the Qt visualizer is built and deployed." -ForegroundColor Yellow
    exit 1
}

# Check if data directory exists
if (-not (Test-Path $DataDir)) {
    Write-Host "ERROR: Data directory not found at: $DataDir" -ForegroundColor Red
    exit 1
}

# Get all .txt files from the data directory
$DataFiles = Get-ChildItem -Path $DataDir -Filter "*.txt" | Sort-Object Name

if ($DataFiles.Count -eq 0) {
    Write-Host "No data files found in: $DataDir" -ForegroundColor Yellow
    exit 0
}

Write-Host "Found $($DataFiles.Count) data files in:" -ForegroundColor Cyan
Write-Host "  $DataDir" -ForegroundColor Gray
Write-Host ""
Write-Host "Visualizer:" -ForegroundColor Cyan
Write-Host "  $VisualizerPath" -ForegroundColor Gray
Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Counter for processed files
$FileNumber = 0

# Process each file one by one
foreach ($File in $DataFiles) {
    $FileNumber++
    
    Write-Host "[$FileNumber/$($DataFiles.Count)] Processing: $($File.Name)" -ForegroundColor Yellow
    
    # Start the visualizer with the current data file
    $Process = Start-Process -FilePath $VisualizerPath -ArgumentList "`"$($File.FullName)`"" -PassThru -Wait
    
    if ($Process.ExitCode -ne 0) {
        Write-Host "  WARNING: Visualizer exited with code $($Process.ExitCode)" -ForegroundColor Red
    } else {
        Write-Host "  Completed successfully" -ForegroundColor Green
    }
    
    Write-Host ""
}

Write-Host "========================================" -ForegroundColor Green
Write-Host "All files processed: $FileNumber" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
