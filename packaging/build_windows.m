function build_windows(version, outDir)
%BUILD_WINDOWS  Build the CADENCE Windows app and installer with MATLAB Compiler.
%
%   build_windows                      % version 1.0.1, output C:\CADENCE_build
%   build_windows('1.0.2')
%   build_windows('1.0.2', 'D:\builds\cadence')
%
%   Run on Windows with MATLAB R2025b and MATLAB Compiler installed
%   (check: which mcc). MATLAB Compiler cannot cross-compile, so the macOS
%   app is built separately on a Mac (Cadence2026.prj compiler task).
%
%   What it does
%     1. Puts ONLY the CADENCE source folders on the path. Build output
%        (distribution/, release/) and Claude worktrees (.claude/) hold full
%        copies of the code, some encrypted; if they were on the path the
%        compiler could bundle them instead of the source.
%     2. Builds CADENCE.exe with the helper folders listed explicitly as
%        additional files. The automatic dependency scan missed
%        GSDconverter_BVW (Data Conversion's .gsd reader) on Windows, so
%        nothing is left to the scan.
%     3. Packages CADENCE_Installer_Windows.exe (MATLAB Runtime downloaded at
%        install time).
%     4. Prints the build-log checks.
%
%   Settings mirror the macOS compiler task: installer splash Logo_v3.png,
%   installer logo Sidebar.png, and the launcher's images bundled. The icon
%   is the square cadence_icon_map_1024.png; the non-square
%   cadence_icon_illustrator.png left the default icon, and MATLAB cannot
%   read the repo's .ico (its images are PNG-compressed).
%
%   After a rebuild in the same folder, File Explorer may keep showing a
%   cached (default) icon for CADENCE.exe; copy it to a new name to check.
%
%   The installer is not code-signed, so Windows SmartScreen shows an
%   "unrecognized app" warning (More info -> Run anyway).
%   See WINDOWS_BUILD_CHECKLIST.md for signing.

    if nargin < 1 || isempty(version), version = '1.0.1'; end
    if nargin < 2 || isempty(outDir),  outDir  = 'C:\CADENCE_build'; end
    if ~ispc
        error('build_windows:platform', 'Run this on Windows; build the macOS app from the Mac compiler task.');
    end
    if isempty(which('mcc'))
        error('build_windows:compiler', 'MATLAB Compiler is not installed (which mcc is empty).');
    end

    root    = fileparts(fileparts(mfilename('fullpath')));   % repo root (this file is in packaging/)
    helpers = {'utils', 'signal_conditioning_helper', 'feature_extraction_helper', ...
               'conduction_velocity_helper', 'arrhythmia_dynamics_helper', 'validation_helper'};

    % ---- 1. clean path: repo root + helper folders only ----------------------
    restoredefaultpath;
    addpath(root);
    for k = 1:numel(helpers)
        addpath(genpath(fullfile(root, helpers{k})));
    end
    hits = which('compute_lat_50', '-all');
    if numel(hits) ~= 1
        error('build_windows:path', 'compute_lat_50 resolves to %d files; expected exactly 1.', numel(hits));
    end
    fprintf('Source   : %s\n', root);

    % ---- 2. build ---------------------------------------------------------
    buildDir = fullfile(outDir, 'build');
    pkgDir   = fullfile(outDir, 'package');
    for d = {buildDir, pkgDir}
        if isfolder(d{1}), rmdir(d{1}, 's'); end            % no stale files from an earlier build
    end
    exeIcon = fullfile(root, 'assets', 'cadence_icon_map_1024.png');
    extra   = [{fullfile(root, 'Logo_v3.png'), fullfile(root, 'Sidebar.png'), ...
                fullfile(root, 'assets', 'cadence_icon_illustrator.png')}, ...
               fullfile(root, helpers)];

    r = compiler.build.standaloneWindowsApplication(fullfile(root, 'Cadence.mlapp'), ...
        'ExecutableName', 'CADENCE', 'ExecutableVersion', version, ...
        'ExecutableIcon', exeIcon, 'AdditionalFiles', extra, ...
        'OutputDir', buildDir, 'Verbose', 'on');

    % ---- 3. installer -------------------------------------------------------
    compiler.package.installer(r, 'ApplicationName', 'CADENCE', 'Version', version, ...
        'InstallerName', 'CADENCE_Installer_Windows', 'AuthorName', 'Kedar Aras', ...
        'InstallerIcon', exeIcon, ...
        'InstallerSplash', fullfile(root, 'Logo_v3.png'), ...
        'InstallerLogo', fullfile(root, 'Sidebar.png'), ...
        'RuntimeDelivery', 'web', 'OutputDir', pkgDir);

    % ---- 4. build-log checks -----------------------------------------------
    fprintf('\nBuild checks (%s):\n', buildDir);
    report(fullfile(buildDir, 'unresolvedSymbols.txt'), 'unresolved symbols');
    report(fullfile(buildDir, 'mccExcludedFiles.log'),  'excluded files');
    fprintf('  installer: %s\n', fullfile(pkgDir, 'CADENCE_Installer_Windows.exe'));
    fprintf('Next: run CADENCE.exe and take one recording through all five modules.\n');
end


function report(file, what)
% A clean log has only its header line(s); anything more is listed.
    if ~isfile(file), fprintf('  %s: log not found\n', what); return; end
    lines = splitlines(strtrim(fileread(file)));
    body  = lines(~startsWith(lines, {'Path', 'Excluded Files'}) & strlength(strtrim(lines)) > 0);
    if isempty(body)
        fprintf('  %s: none (clean)\n', what);
    else
        fprintf('  %s: %d line(s), see %s\n', what, numel(body), file);
    end
end
