function mexPath = build_sos_nlp_mex_nag_mac()
%BUILD_SOS_NLP_MEX_NAG_MAC Build SOS MEX with NAG Fortran on macOS.

root = fileparts(fileparts(mfilename("fullpath")));

baseMexSrc = fullfile(root, "src", "mex", "sos_nlp_mex_clean.cpp");
ioSrc = fullfile(root, "src", "fortran", "sos_mex_io.f");
sosArchive = fullfile(root, "build", "libsos_nlp_nag_mac.a");
highsArchive = fullfile(root, "third_party", "highs-install-source-mac", ...
    "lib", "libhighs.a");
nagQuickfit = "/usr/local/lib/NAG_Fortran/quickfit.o";
nagRuntime = "/usr/local/lib/NAG_Fortran/libf72rts.a";

outDir = fullfile(root, "build", "mex", "nag_mac", mexext);
objDir = fullfile(outDir, "obj");
mexSrc = fullfile(outDir, "sos_nlp_mex_nag_mac.cpp");
ioObj = fullfile(objDir, "sos_mex_io.o");

unagfor = "/usr/local/bin/unagfor";
if strlength(getenv("NAG_KUSARI_FILE")) == 0
    error("sos:build:missingNagLicense", ...
        "Set NAG_KUSARI_FILE to your NAG license file before building.");
end

mustExist(baseMexSrc, "clean MEX source");
mustExist(ioSrc, "Fortran MEX I/O helper");
mustExist(sosArchive, "NAG SOS archive");
mustExist(highsArchive, "HiGHS static archive");
mustExist(nagQuickfit, "NAG Fortran quickfit object");
mustExist(nagRuntime, "NAG Fortran static runtime archive");

if ~isfolder(objDir)
    mkdir(objDir);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

text = fileread(baseMexSrc);
text = strrep(text, "gfortran + HiGHS", ...
    "NAG Fortran macOS + source-built HiGHS");
text = strrep(text, ...
    "sos_nlp_highs_experiment/build/libsos_nlp.a + HiGHS", ...
    "build/libsos_nlp_nag_mac.a + source-built HiGHS");
fid = fopen(mexSrc, "w");
if fid < 0
    error("sos:build:openFailed", "Could not write %s", mexSrc);
end
cleanupFile = onCleanup(@() fclose(fid));
fprintf(fid, "%s", text);
clear cleanupFile

cmd = sprintf('%s -O2 -dusty -PIC -c "%s" -o "%s"', ...
    unagfor, ioSrc, ioObj);
runCommand(cmd);

oldDir = cd(outDir);
cleanup = onCleanup(@() cd(oldDir));

mex("-R2018a", "-v", "-output", "sos_nlp_mex", mexSrc, ...
    nagQuickfit, ioObj, sosArchive, highsArchive, ...
    "LDFLAGS=$LDFLAGS -Wl,-force_load," + nagRuntime + " -Wl,-dead_strip", ...
    "-lz", "-lm");

mexPath = fullfile(outDir, "sos_nlp_mex." + mexext);
fprintf("Built %s\n", mexPath);
end

function runCommand(cmd)
[status, out] = system(cmd);
if status ~= 0
    error("sos:build:commandFailed", "Command failed: %s\n%s", cmd, out);
end
end

function mustExist(pathValue, label)
if ~isfile(pathValue)
    error("sos:build:missingFile", "Missing %s: %s", label, pathValue);
end
end
