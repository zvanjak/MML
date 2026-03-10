# Unregister MML File Association
# Run this script as Administrator to remove .mml file extension registration
# Copyright (c) 2026 MML Visualizers

param(
    [switch]$CurrentUser  # Unregister for current user only
)

$ErrorActionPreference = "Stop"

# Determine registry root
if ($CurrentUser) {
    $regRoot = "HKCU:\Software\Classes"
    Write-Host "Unregistering for current user..." -ForegroundColor Yellow
} else {
    $regRoot = "HKLM:\Software\Classes"
    Write-Host "Unregistering for all users (requires Administrator)..." -ForegroundColor Yellow
    
    # Check for admin privileges
    $isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $isAdmin) {
        Write-Host ""
        Write-Host "This script requires Administrator privileges for system-wide unregistration." -ForegroundColor Red
        Write-Host "Please run PowerShell as Administrator, or use -CurrentUser flag:" -ForegroundColor Yellow
        Write-Host "    .\unregister_mml.ps1 -CurrentUser" -ForegroundColor Cyan
        Write-Host ""
        exit 1
    }
}

$progId = "MML.VisualizerFile"
$progIdPath = "$regRoot\$progId"
$extPath = "$regRoot\.mml"

Write-Host ""
Write-Host "Removing .mml file association..." -ForegroundColor Cyan

$removed = $false

try {
    # Remove extension registration
    if (Test-Path $extPath) {
        $currentProgId = (Get-ItemProperty -Path $extPath -ErrorAction SilentlyContinue).'(default)'
        if ($currentProgId -eq $progId) {
            Write-Host "  Removing .mml extension registration"
            Remove-Item -Path $extPath -Recurse -Force
            $removed = $true
        } else {
            Write-Host "  .mml extension is registered to a different program: $currentProgId" -ForegroundColor Yellow
            Write-Host "  Skipping extension removal to avoid breaking other applications."
        }
    } else {
        Write-Host "  .mml extension not registered" -ForegroundColor Gray
    }
    
    # Remove ProgID
    if (Test-Path $progIdPath) {
        Write-Host "  Removing ProgID: $progId"
        Remove-Item -Path $progIdPath -Recurse -Force
        $removed = $true
    } else {
        Write-Host "  ProgID not found" -ForegroundColor Gray
    }
    
    if ($removed) {
        # Notify shell of changes
        Write-Host "  Notifying shell of changes..."
        
        Add-Type -TypeDefinition @"
            using System;
            using System.Runtime.InteropServices;
            public class Shell32Unreg {
                [DllImport("shell32.dll", CharSet = CharSet.Auto, SetLastError = true)]
                public static extern void SHChangeNotify(int wEventId, int uFlags, IntPtr dwItem1, IntPtr dwItem2);
            }
"@
        [Shell32Unreg]::SHChangeNotify(0x08000000, 0, [IntPtr]::Zero, [IntPtr]::Zero)
        
        Write-Host ""
        Write-Host "========================================" -ForegroundColor Green
        Write-Host "Unregistration complete!" -ForegroundColor Green
        Write-Host "========================================" -ForegroundColor Green
        Write-Host ""
    } else {
        Write-Host ""
        Write-Host "Nothing to unregister." -ForegroundColor Yellow
        Write-Host ""
    }
    
} catch {
    Write-Host ""
    Write-Host "Error during unregistration: $_" -ForegroundColor Red
    Write-Host ""
    exit 1
}
