function mexPath = build_sos_nlp_mex_gfortran_mac()
%BUILD_SOS_NLP_MEX_GFORTRAN_MAC Build SOS MEX with gfortran on macOS.

root = fileparts(fileparts(mfilename("fullpath")));

baseMexSrc = fullfile(root, "src", "mex", "sos_nlp_mex_clean.cpp");
ioSrc = fullfile(root, "src", "fortran", "sos_mex_io.f");
sosArchive = fullfile(root, "build", "libsos_nlp_gfortran_mac.a");
highsArchive = fullfile(root, "third_party", "highs-install-source-mac", ...
    "lib", "libhighs.a");

outDir = fullfile(root, "build", "mex", "gfortran_mac", mexext);
objDir = fullfile(outDir, "obj");
mexSrc = fullfile(outDir, "sos_nlp_mex_gfortran_mac.cpp");
ioObj = fullfile(objDir, "sos_mex_io.o");

gfortran = "/opt/homebrew/bin/gfortran";
libgfortran = strip(string(systemCapture(gfortran + ...
    " -print-file-name=libgfortran.a")));
libquadmath = strip(string(systemCapture(gfortran + ...
    " -print-file-name=libquadmath.a")));
libgcc = strip(string(systemCapture(gfortran + ...
    " -print-file-name=libgcc.a")));

mustExist(baseMexSrc, "clean MEX source");
mustExist(ioSrc, "Fortran MEX I/O helper");
mustExist(sosArchive, "SOS archive");
mustExist(highsArchive, "HiGHS static archive");
mustExist(libgfortran, "gfortran runtime archive");
mustExist(libquadmath, "quadmath runtime archive");
mustExist(libgcc, "gcc runtime archive");

if ~isfolder(objDir)
    mkdir(objDir);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

text = fileread(baseMexSrc);
text = strrep(text, "gfortran + HiGHS", ...
    "gfortran macOS + source-built HiGHS");
text = strrep(text, ...
    "sos_nlp_highs_experiment/build/libsos_nlp.a + HiGHS", ...
    "build/libsos_nlp_gfortran_mac.a + source-built HiGHS");
fid = fopen(mexSrc, "w");
if fid < 0
    error("sos:build:openFailed", "Could not write %s", mexSrc);
end
cleanupFile = onCleanup(@() fclose(fid));
fprintf(fid, "%s", text);
clear cleanupFile

cmd = sprintf('%s -O2 -std=legacy -fallow-argument-mismatch -mmacosx-version-min=12.0 -fPIC -c "%s" -o "%s"', ...
    gfortran, ioSrc, ioObj);
runCommand(cmd);

oldDir = cd(outDir);
cleanup = onCleanup(@() cd(oldDir));

mex("-R2018a", "-v", "-output", "sos_nlp_mex", mexSrc, ...
    ioObj, sosArchive, highsArchive, ...
    "LDFLAGS=$LDFLAGS -Wl,-dead_strip", ...
    libgfortran, libquadmath, libgcc, "-lz");

mexPath = fullfile(outDir, "sos_nlp_mex." + mexext);
fprintf("Built %s\n", mexPath);
end

function out = systemCapture(cmd)
[status, out] = system(cmd);
if status ~= 0
    error("sos:build:commandFailed", "Command failed: %s\n%s", cmd, out);
end
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
