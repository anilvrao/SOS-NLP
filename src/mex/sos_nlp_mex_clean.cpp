#include "mex.h"

#include <algorithm>
#include <cctype>
#include <ctime>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

extern "C" {
void sprnlp_(int*, int*, double*, double*, double*, int*, double*, int*,
             double*, double*, int*, int*, double*, double*, int*, int*, int*,
             double*, int*, int*, int*, double*, double*, double*, int*, int*,
             int*, double*, double*, int*, int*, int*, int*, int*, double*,
             int*, int*, int*, int*, int*);

void barnlp_(int*, int*, double*, double*, double*, int*, double*, int*,
             double*, double*, int*, int*, double*, double*, int*, int*, int*,
             double*, int*, int*, int*, double*, double*, double*, int*, int*,
             int*, double*, double*, int*, int*, int*, int*, int*, double*,
             int*, int*, int*, int*, int*);

void insnlp_(char*, int);
void sosopn_(int*, char*, int);
void soscls_(int*);
}

extern "C" float etime_(float* tarray) {
    const float elapsed = static_cast<float>(std::clock()) /
                          static_cast<float>(CLOCKS_PER_SEC);
    if (tarray) {
        tarray[0] = elapsed;
        tarray[1] = 0.0f;
    }
    return elapsed;
}

namespace {

void fail(const char* id, const char* msg) {
    mexErrMsgIdAndTxt(id, "%s", msg);
}

std::string getString(const mxArray* value, const char* name, bool lowercase = true) {
    const mxArray* charValue = value;
    mxArray* converted = nullptr;
    if (mxIsClass(value, "string")) {
        if (mxGetNumberOfElements(value) != 1) {
            fail("sos_nlp_mex:type", name);
        }
        mxArray* rhs = const_cast<mxArray*>(value);
        if (mexCallMATLAB(1, &converted, 1, &rhs, "char") != 0) {
            fail("sos_nlp_mex:string", "Could not convert MATLAB string.");
        }
        charValue = converted;
    }
    if (!mxIsChar(charValue)) {
        if (converted) {
            mxDestroyArray(converted);
        }
        fail("sos_nlp_mex:type", name);
    }
    char* raw = mxArrayToString(charValue);
    if (converted) {
        mxDestroyArray(converted);
    }
    if (!raw) {
        fail("sos_nlp_mex:string", "Could not convert MATLAB string.");
    }
    std::string out(raw);
    mxFree(raw);
    if (lowercase) {
        std::transform(out.begin(), out.end(), out.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    }
    return out;
}

mxArray* makeStringCell(const std::vector<const char*>& values) {
    mxArray* cell = mxCreateCellMatrix(1, values.size());
    for (mwIndex i = 0; i < values.size(); ++i) {
        mxSetCell(cell, i, mxCreateString(values[i]));
    }
    return cell;
}

mxArray* makeInfo() {
    const char* names[] = {
        "mex", "platform", "matlabRelease", "fortranCompiler", "solvers",
        "interface", "archive", "optimalControlExcluded"
    };
    mxArray* out = mxCreateStructMatrix(1, 1, 8, names);
    mxSetField(out, 0, "mex", mxCreateString("sos_nlp_mex"));
    mxSetField(out, 0, "platform", mxCreateString("maca64"));
    mxSetField(out, 0, "matlabRelease", mxCreateString("R2025b"));
    mxSetField(out, 0, "fortranCompiler",
               mxCreateString("gfortran + HiGHS"));
    mxSetField(out, 0, "solvers",
               makeStringCell({"sprnlp", "barnlp"}));
    mxSetField(out, 0, "interface",
               mxCreateString("reverse-communication one-step gateway"));
    mxSetField(out, 0, "archive",
               mxCreateString("sos_nlp_highs_experiment/build/libsos_nlp.a + HiGHS"));
    mxSetField(out, 0, "optimalControlExcluded", mxCreateLogicalScalar(true));
    return out;
}

void callInsnlp(const std::string& command) {
    constexpr int len = 80;
    char buf[len];
    std::memset(buf, ' ', len);
    const std::size_t n = std::min(command.size(), static_cast<std::size_t>(len));
    std::memcpy(buf, command.data(), n);
    insnlp_(buf, len);
}

int getScalarIntInput(const mxArray* value, const char* message) {
    if (!mxIsDouble(value) || mxIsComplex(value) || mxGetNumberOfElements(value) != 1) {
        fail("sos_nlp_mex:type", message);
    }
    return static_cast<int>(mxGetScalar(value));
}

void openFortranOutput(int unit, const std::string& path) {
    if (path.empty()) {
        fail("sos_nlp_mex:outputPath", "Output path must not be empty.");
    }
    std::vector<char> chars(path.begin(), path.end());
    sosopn_(&unit, chars.data(), static_cast<int>(chars.size()));
}

void closeFortranOutput(int unit) {
    soscls_(&unit);
}

void resetGateway(const std::vector<int>& units) {
    for (int unit : units) {
        soscls_(&unit);
    }
}

mxArray* requireField(mxArray* state, const char* name) {
    mxArray* value = mxGetField(state, 0, name);
    if (!value) {
        mexErrMsgIdAndTxt("sos_nlp_mex:missingField",
                          "State struct is missing field '%s'.", name);
    }
    return value;
}

int scalarInt(mxArray* state, const char* name) {
    mxArray* value = requireField(state, name);
    if (!mxIsDouble(value) || mxIsComplex(value) || mxGetNumberOfElements(value) != 1) {
        mexErrMsgIdAndTxt("sos_nlp_mex:type",
                          "Field '%s' must be a real scalar double.", name);
    }
    return static_cast<int>(mxGetScalar(value));
}

double* realData(mxArray* state, const char* name, mwSize minElements) {
    mxArray* value = requireField(state, name);
    if (!mxIsDouble(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) < minElements) {
        mexErrMsgIdAndTxt("sos_nlp_mex:type",
                          "Field '%s' must be a real double array with enough elements.",
                          name);
    }
    return mxGetDoubles(value);
}

int* intData(mxArray* state, const char* name, mwSize minElements) {
    mxArray* value = requireField(state, name);
    if (!mxIsInt32(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) < minElements) {
        mexErrMsgIdAndTxt("sos_nlp_mex:type",
                          "Field '%s' must be a real int32 array with enough elements.",
                          name);
    }
    return static_cast<int*>(mxGetData(value));
}

int* scalarIntData(mxArray* state, const char* name) {
    mxArray* value = requireField(state, name);
    if (!mxIsInt32(value) || mxIsComplex(value) || mxGetNumberOfElements(value) != 1) {
        mexErrMsgIdAndTxt("sos_nlp_mex:type",
                          "Field '%s' must be a scalar int32.", name);
    }
    return static_cast<int*>(mxGetData(value));
}

mxArray* stepSolver(const std::string& solver, const mxArray* inputState) {
    if (!mxIsStruct(inputState) || mxGetNumberOfElements(inputState) != 1) {
        fail("sos_nlp_mex:type", "State must be a scalar struct.");
    }

    mxArray* state = mxDuplicateArray(inputState);

    int* irvcom = scalarIntData(state, "irvcom");
    int* irevrs = intData(state, "irevrs", 5);
    int* ndimP = scalarIntData(state, "ndim");
    int* maxresP = scalarIntData(state, "maxres");
    int* nresP = scalarIntData(state, "nres");
    int* nonzrP = scalarIntData(state, "nonzr");
    int* nonzhP = scalarIntData(state, "nonzh");
    int* maxconP = scalarIntData(state, "maxcon");
    int* mconP = scalarIntData(state, "mcon");
    int* nonzgP = scalarIntData(state, "nonzg");
    int* iferr = scalarIntData(state, "iferr");
    int* nfeval = scalarIntData(state, "nfeval");
    int* nholdP = scalarIntData(state, "nhold");
    int* niholdP = scalarIntData(state, "nihold");
    int* needed = scalarIntData(state, "needed");
    int* iernlp = scalarIntData(state, "iernlp");

    const int ndim = *ndimP;
    const int maxres = *maxresP;
    const int nonzr = *nonzrP;
    const int nonzh = *nonzhP;
    const int maxcon = *maxconP;
    const int nonzg = *nonzgP;
    const int nhold = *nholdP;
    const int nihold = *niholdP;

    if (ndim < 1 || maxres < 1 || maxcon < 1 || nonzr < 1 ||
        nonzh < 1 || nonzg < 1 || nhold < 1 || nihold < 1) {
        mxDestroyArray(state);
        fail("sos_nlp_mex:dimensions",
             "ndim, maxres, maxcon, nonzr, nonzh, nonzg, nhold, and nihold must be positive.");
    }

    double* xbar = realData(state, "xbar", ndim);
    double* xlwr = realData(state, "xlwr", ndim);
    double* xupr = realData(state, "xupr", ndim);
    int* istatv = intData(state, "istatv", ndim);
    double* vecnu = realData(state, "vecnu", ndim);
    double* fbar = realData(state, "fbar", 1);
    double* resvec = realData(state, "resvec", maxres);
    double* delf = realData(state, "delf", ndim);
    double* rmat = realData(state, "rmat", nonzr);
    int* irowr = intData(state, "irowr", nonzr);
    int* jcolr = intData(state, "jcolr", nonzr);
    double* hmat = realData(state, "hmat", nonzh);
    int* irowh = intData(state, "irowh", nonzh);
    int* jstrh = intData(state, "jstrh", ndim + 1);
    double* cbar = realData(state, "cbar", maxcon);
    double* clwr = realData(state, "clwr", maxcon);
    double* cupr = realData(state, "cupr", maxcon);
    int* istatc = intData(state, "istatc", maxcon);
    double* veclam = realData(state, "veclam", maxcon);
    double* gmat = realData(state, "gmat", nonzg);
    int* irowg = intData(state, "irowg", nonzg);
    int* jcolg = intData(state, "jcolg", nonzg);
    double* hold = realData(state, "hold", nhold);
    int* ihold = intData(state, "ihold", nihold);

    if (solver == "sprnlp") {
        sprnlp_(irvcom, irevrs, xbar, xlwr, xupr, istatv, vecnu, ndimP,
                fbar, resvec, maxresP, nresP, delf, rmat, irowr, jcolr,
                nonzrP, hmat, irowh, jstrh, nonzhP, cbar, clwr, cupr,
                istatc, maxconP, mconP, veclam, gmat, irowg, jcolg,
                nonzgP, iferr, nfeval, hold, nholdP, ihold, niholdP,
                needed, iernlp);
    } else if (solver == "barnlp") {
        barnlp_(irvcom, irevrs, xbar, xlwr, xupr, istatv, vecnu, ndimP,
                fbar, resvec, maxresP, nresP, delf, rmat, irowr, jcolr,
                nonzrP, hmat, irowh, jstrh, nonzhP, cbar, clwr, cupr,
                istatc, maxconP, mconP, veclam, gmat, irowg, jcolg,
                nonzgP, iferr, nfeval, hold, nholdP, ihold, niholdP,
                needed, iernlp);
    } else {
        mxDestroyArray(state);
        fail("sos_nlp_mex:solver", "Solver must be 'sprnlp' or 'barnlp'.");
    }

    return state;
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1) {
        plhs[0] = makeInfo();
        return;
    }

    const std::string command = getString(prhs[0], "First input must be a command string.");

    if (command == "info") {
        if (nlhs > 1) {
            fail("sos_nlp_mex:nlhs", "The 'info' command returns one output.");
        }
        plhs[0] = makeInfo();
        return;
    }

    if (command == "solvers") {
        if (nlhs > 1) {
            fail("sos_nlp_mex:nlhs", "The 'solvers' command returns one output.");
        }
        plhs[0] = makeStringCell({"sprnlp", "barnlp"});
        return;
    }

    if (command == "insnlp") {
        if (nrhs != 2 || nlhs > 0) {
            fail("sos_nlp_mex:usage", "Usage: sos_nlp_mex('insnlp', command)");
        }
        callInsnlp(getString(prhs[1], "INSNLP command must be a string.", false));
        return;
    }

    if (command == "open_output") {
        if (nrhs != 3 || nlhs > 0) {
            fail("sos_nlp_mex:usage", "Usage: sos_nlp_mex('open_output', unit, path)");
        }
        const int unit = getScalarIntInput(prhs[1], "Output unit must be a scalar double.");
        openFortranOutput(unit, getString(prhs[2], "Output path must be a string.", false));
        return;
    }

    if (command == "close_output") {
        if (nrhs != 2 || nlhs > 0) {
            fail("sos_nlp_mex:usage", "Usage: sos_nlp_mex('close_output', unit)");
        }
        const int unit = getScalarIntInput(prhs[1], "Output unit must be a scalar double.");
        closeFortranOutput(unit);
        return;
    }

    if (command == "reset") {
        if (nlhs > 0) {
            fail("sos_nlp_mex:nlhs", "The 'reset' command returns no outputs.");
        }
        std::vector<int> units;
        if (nrhs == 1) {
            units = {99};
        } else {
            for (int k = 1; k < nrhs; ++k) {
                units.push_back(getScalarIntInput(prhs[k], "Reset units must be scalar doubles."));
            }
        }
        resetGateway(units);
        return;
    }

    if (command == "step") {
        if (nrhs != 3 || nlhs > 1) {
            fail("sos_nlp_mex:usage",
                 "Usage: state = sos_nlp_mex('step', solver, state)");
        }
        const std::string solver = getString(prhs[1], "Solver must be a string.");
        plhs[0] = stepSolver(solver, prhs[2]);
        return;
    }

    fail("sos_nlp_mex:command", "Unknown command. Use 'info', 'solvers', 'reset', or 'step'.");
}
