# Register MML File Association
# Run this script as Administrator to register .mml file extension
# Copyright (c) 2026 MML Visualizers

param(
    [switch]$CurrentUser,  # Register for current user only (no admin required)
    [switch]$Force         # Overwrite existing registration
)

$ErrorActionPreference = "Stop"

# Get the launcher path
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$launcherPath = Join-Path $scriptDir "MML_Launcher.exe"

# Check if launcher exists
if (-not (Test-Path $launcherPath)) {
    # Try in parent directory
    $launcherPath = Join-Path (Split-Path -Parent $scriptDir) "MML_Launcher.exe"
    if (-not (Test-Path $launcherPath)) {
        Write-Error "MML_Launcher.exe not found. Please run this script from the same directory as the launcher."
        exit 1
    }
}

$launcherPath = (Resolve-Path $launcherPath).Path
Write-Host "Launcher found at: $launcherPath" -ForegroundColor Green

# Determine registry root
if ($CurrentUser) {
    $regRoot = "HKCU:\Software\Classes"
    Write-Host "Registering for current user only..." -ForegroundColor Yellow
} else {
    $regRoot = "HKLM:\Software\Classes"
    Write-Host "Registering for all users (requires Administrator)..." -ForegroundColor Yellow
    
    # Check for admin privileges
    $isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $isAdmin) {
        Write-Host ""
        Write-Host "This script requires Administrator privileges for system-wide registration." -ForegroundColor Red
        Write-Host "Please run PowerShell as Administrator, or use -CurrentUser flag:" -ForegroundColor Yellow
        Write-Host "    .\register_mml.ps1 -CurrentUser" -ForegroundColor Cyan
        Write-Host ""
        exit 1
    }
}

# ProgID for MML files
$progId = "MML.VisualizerFile"
$progIdPath = "$regRoot\$progId"
$extPath = "$regRoot\.mml"

# Check for existing registration
if ((Test-Path $extPath) -and -not $Force) {
    $existingProgId = (Get-ItemProperty -Path $extPath -ErrorAction SilentlyContinue).'(default)'
    if ($existingProgId) {
        Write-Host ""
        Write-Host "Warning: .mml extension is already registered to: $existingProgId" -ForegroundColor Yellow
        Write-Host "Use -Force to overwrite the existing registration." -ForegroundColor Yellow
        Write-Host ""
        $response = Read-Host "Continue anyway? (y/N)"
        if ($response -ne 'y' -and $response -ne 'Y') {
            Write-Host "Registration cancelled."
            exit 0
        }
    }
}

Write-Host ""
Write-Host "Registering .mml file association..." -ForegroundColor Cyan

try {
    # Create ProgID
    Write-Host "  Creating ProgID: $progId"
    New-Item -Path $progIdPath -Force | Out-Null
    Set-ItemProperty -Path $progIdPath -Name "(default)" -Value "MML Visualizer Data File"
    
    # Set default icon (use the launcher's icon)
    $iconPath = "$progIdPath\DefaultIcon"
    Write-Host "  Setting default icon"
    New-Item -Path $iconPath -Force | Out-Null
    Set-ItemProperty -Path $iconPath -Name "(default)" -Value "`"$launcherPath`",0"
    
    # Create shell\open\command
    $shellPath = "$progIdPath\shell"
    $openPath = "$shellPath\open"
    $commandPath = "$openPath\command"
    
    Write-Host "  Creating shell open command"
    New-Item -Path $commandPath -Force | Out-Null
    Set-ItemProperty -Path $openPath -Name "(default)" -Value "Open with MML Visualizer"
    Set-ItemProperty -Path $commandPath -Name "(default)" -Value "`"$launcherPath`" `"%1`""
    
    # Register .mml extension
    Write-Host "  Registering .mml extension"
    New-Item -Path $extPath -Force | Out-Null
    Set-ItemProperty -Path $extPath -Name "(default)" -Value $progId
    Set-ItemProperty -Path $extPath -Name "Content Type" -Value "application/x-mml-visualizer"
    
    # Notify shell of changes
    Write-Host "  Notifying shell of changes..."
    
    # Use SHChangeNotify to refresh shell
    Add-Type -TypeDefinition @"
        using System;
        using System.Runtime.InteropServices;
        public class Shell32 {
            [DllImport("shell32.dll", CharSet = CharSet.Auto, SetLastError = true)]
            public static extern void SHChangeNotify(int wEventId, int uFlags, IntPtr dwItem1, IntPtr dwItem2);
        }
"@
    [Shell32]::SHChangeNotify(0x08000000, 0, [IntPtr]::Zero, [IntPtr]::Zero)  # SHCNE_ASSOCCHANGED

    Write-Host ""
    Write-Host "========================================" -ForegroundColor Green
    Write-Host "Registration complete!" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Green
    Write-Host ""
    Write-Host "The .mml file extension is now associated with MML Launcher."
    Write-Host "Double-clicking any .mml file will open the appropriate visualizer."
    Write-Host ""
    Write-Host "Supported file types:"
    Write-Host "  - REAL_FUNCTION, MULTI_REAL_FUNCTION"
    Write-Host "  - PARAMETRIC_CURVE_CARTESIAN_2D/3D"
    Write-Host "  - PARAMETRIC_SURFACE_CARTESIAN"
    Write-Host "  - PARTICLE_SIMULATION_DATA_2D/3D"
    Write-Host "  - VECTOR_FIELD_2D/3D_CARTESIAN"
    Write-Host "  - SCALAR_FUNCTION_CARTESIAN_2D/3D"
    Write-Host "  - RIGID_BODY_TRAJECTORY_3D"
    Write-Host ""
    
} catch {
    Write-Host ""
    Write-Host "Error during registration: $_" -ForegroundColor Red
    Write-Host ""
    exit 1
}
