function mexPath = build_sos_nlp_mex_mingw_windows()
%BUILD_SOS_NLP_MEX_MINGW_WINDOWS Build SOS MEX with MSYS2/MINGW on Windows.
%
% This wrapper invokes the MSYS2/MINGW source build script.  The generated
% MEX has the standard SOS gateway name, sos_nlp_mex.mexw64.

root = fileparts(fileparts(mfilename("fullpath")));
bash = "C:\msys64\usr\bin\bash.exe";
script = fullfile(root, "tools", "build_sos_nlp_mex_mingw_windows.sh");

if ~isfile(bash)
    error("sos:build:missingBash", "Missing MSYS2 bash: %s", bash);
end
if ~isfile(script)
    error("sos:build:missingScript", "Missing build script: %s", script);
end

cmd = sprintf('"%s" "%s"', bash, script);
[status, out] = system(cmd);
if status ~= 0
    error("sos:build:commandFailed", "Command failed: %s\n%s", cmd, out);
end

mexPath = fullfile(root, "build", "mex", "mingw_windows", "mexw64", ...
    "sos_nlp_mex.mexw64");
fprintf("%s\n", out);
fprintf("Built %s\n", mexPath);
end
