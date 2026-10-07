# .NET And WPF Notices

MML Visualizers WPF release archives contain Windows-only WPF applications targeting .NET 8.

## Current Packaging Model

The current Windows release gate publishes WPF applications with:

```text
dotnet publish -c Release -f net8.0-windows --no-self-contained
```

That produces framework-dependent applications and does not intentionally redistribute the .NET
Desktop Runtime inside the archive. Users must install the Microsoft .NET 8 Desktop Runtime on
Windows.

## License

.NET runtime components and reference assemblies are licensed by Microsoft under their own terms,
including MIT-licensed .NET source components and Microsoft runtime distribution terms where
applicable. They are not relicensed by MML Visualizers.

## Required Release Check

If a future WPF artifact becomes self-contained, single-file, or otherwise redistributes .NET runtime
components, include the Microsoft .NET license and third-party notices from the exact SDK/runtime
used to publish that artifact.

Official .NET licensing information:

- https://github.com/dotnet/runtime/blob/main/LICENSE.TXT
- https://github.com/dotnet/runtime/blob/main/THIRD-PARTY-NOTICES.TXT
- https://dotnet.microsoft.com/platform/free
