#include "noise_simulations/param_sweep.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace lbcrypto;

Context generateContext(const HEParams& p) {
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetPlaintextModulus(p.plaintextModulus);
    parameters.SetMultiplicativeDepth(p.multDepth);
    parameters.SetStandardDeviation(static_cast<float>(p.sigma));
    // Unscaled BGV: v = c0 + c1*s = m + t*e (mod Q) holds exactly (see crypto_context.cpp).
    parameters.SetScalingTechnique(FIXEDMANUAL);
    if (p.firstModSize > 0) parameters.SetFirstModSize(p.firstModSize);
    if (p.scalingModSize > 0) parameters.SetScalingModSize(p.scalingModSize);
    if (p.numLargeDigits > 0) parameters.SetNumLargeDigits(p.numLargeDigits);

    if (p.enforce128BitSecurity) {
        parameters.SetSecurityLevel(HEStd_128_classic);
        if (p.ringDim > 0) parameters.SetRingDim(p.ringDim);  // OpenFHE throws if n is too small for q
    } else {
        parameters.SetSecurityLevel(HEStd_NotSet);
        parameters.SetRingDim(p.ringDim > 0 ? p.ringDim : 1024);
    }

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    KeyPair<DCRTPoly> keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    return Context{cc, keyPair, p.multDepth};
}

namespace {

std::string autoLabel(const HEParams& p) {
    std::ostringstream os;
    os << "n=" << (p.ringDim ? std::to_string(p.ringDim) : std::string("auto")) << " depth=" << p.multDepth
       << " sigma=" << p.sigma << " p=" << p.plaintextModulus;
    if (!p.enforce128BitSecurity) os << " (insecure)";
    return os.str();
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (quoted) {
            if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
            else if (c == '"') quoted = false;
            else cur += c;
        } else if (c == '"') {
            quoted = true;
        } else if (c == ',') {
            out.push_back(cur);
            cur.clear();
        } else if (c != '\r') {
            cur += c;
        }
    }
    out.push_back(cur);
    return out;
}

std::string csvQuote(const std::string& s) {
    std::string out = "\"";
    for (char c : s) out += (c == '"') ? std::string("\"\"") : std::string(1, c);
    return out + "\"";
}

}  // namespace

std::vector<HEParams> loadHEParamsCSV(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) throw std::runtime_error("loadHEParamsCSV: cannot open " + path);

    std::string line;
    std::vector<std::string> header;
    while (std::getline(file, line)) {  // first non-comment line is the header
        if (!line.empty() && line[0] != '#') { header = splitCsvLine(line); break; }
    }

    std::vector<HEParams> params;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto cells = splitCsvLine(line);
        HEParams p;
        for (size_t i = 0; i < header.size() && i < cells.size(); ++i) {
            const std::string& key = header[i];
            const std::string& v = cells[i];
            if (v.empty()) continue;
            if (key == "label") p.label = v;
            else if (key == "ringDim") p.ringDim = static_cast<uint32_t>(std::stoul(v));
            else if (key == "multDepth") p.multDepth = static_cast<uint32_t>(std::stoul(v));
            else if (key == "firstModSize") p.firstModSize = static_cast<uint32_t>(std::stoul(v));
            else if (key == "scalingModSize") p.scalingModSize = static_cast<uint32_t>(std::stoul(v));
            else if (key == "sigma") p.sigma = std::stod(v);
            else if (key == "plaintextModulus") p.plaintextModulus = static_cast<uint32_t>(std::stoul(v));
            else if (key == "secure") p.enforce128BitSecurity = (v == "1" || v == "true");
            else if (key == "numLargeDigits") p.numLargeDigits = static_cast<uint32_t>(std::stoul(v));
        }
        params.push_back(p);
    }
    std::cout << "Loaded " << params.size() << " parameter sets from " << path << std::endl;
    return params;
}

void runParameterSweep(const std::vector<HEParams>& paramSets, const SweepConfig& cfg) {
    namespace fs = std::filesystem;
    fs::create_directories(cfg.outDir);

    const std::string indexPath = (fs::path(cfg.outDir) / "sweep_index.csv").string();
    std::ofstream index(indexPath);
    index << "id,label,status,ringDim,log2Q,numTowers,sigma,plaintextModulus,multDepth,firstModSize,"
             "scalingModSize,secure,modReduce,numChains,maxMults,growthFile,histFile,numLargeDigits,log2P\n";

    for (size_t k = 0; k < paramSets.size(); ++k) {
        const HEParams& p = paramSets[k];
        std::ostringstream idStream;
        idStream << "p" << std::setw(2) << std::setfill('0') << k;
        const std::string id = idStream.str();
        const std::string label = p.label.empty() ? autoLabel(p) : p.label;
        std::cout << "\n===== [" << id << "] " << label << " =====" << std::endl;

        Context ctx;
        try {
            ctx = generateContext(p);
        } catch (const std::exception& e) {
            std::cerr << "Skipping " << id << ": " << e.what() << std::endl;
            std::string msg = e.what();
            const auto pos = msg.rfind("(): ");  // drop OpenFHE's "file:line:function(): " prefix
            if (pos != std::string::npos) msg = msg.substr(pos + 4);
            msg = msg.substr(0, 200);
            index << id << "," << csvQuote(label) << "," << csvQuote("error: " + msg) << ",,,,"
                  << p.sigma << "," << p.plaintextModulus << "," << p.multDepth << "," << p.firstModSize << ","
                  << p.scalingModSize << "," << (p.enforce128BitSecurity ? 1 : 0) << ","
                  << (cfg.modReduceAfterMult ? 1 : 0) << "," << cfg.numChains << "," << cfg.maxMults << ",,,"
                  << p.numLargeDigits << ",\n";
            index.flush();
            continue;
        }
        printContextDetails(ctx);

        OpChainConfig chainCfg;
        chainCfg.name = label;
        chainCfg.ops = repeatOp(HeOp::Mult, cfg.maxMults);
        chainCfg.numChains = cfg.numChains;
        chainCfg.n = cfg.n;
        chainCfg.modReduceAfterMult = cfg.modReduceAfterMult;
        chainCfg.measureStride = 1;  // every step, so the first error is found exactly
        chainCfg.stopAtFailure = false;
        chainCfg.stepsAfterFailure = cfg.stepsAfterFailure;
        chainCfg.histBinWidth = cfg.histBinWidth;

        const auto rows = runOpChainExperiment(ctx, chainCfg);

        const std::string growthFile = id + "_growth.csv";
        const std::string histFile = id + "_hist.csv";
        writeOpChainCSV((fs::path(cfg.outDir) / growthFile).string(), ctx, rows);
        writeNoiseHistogramCSV((fs::path(cfg.outDir) / histFile).string(), ctx, rows);

        const auto params = ctx.cc->GetCryptoParameters();
        const auto elementParams = params->GetElementParams();
        double sigma = p.sigma;
        if (auto rlwe = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(params))
            sigma = rlwe->GetDistributionParameter();
        std::ostringstream log2Q;  // separate stream, so the precision does not leak into sigma
        log2Q << std::fixed << std::setprecision(1) << std::log2(elementParams->GetModulus().ConvertToDouble());

        index << id << "," << csvQuote(label) << ",ok," << ctx.cc->GetRingDimension() << "," << log2Q.str() << ","
              << elementParams->GetParams().size() << "," << sigma << "," << params->GetPlaintextModulus() << ","
              << p.multDepth << "," << p.firstModSize << "," << p.scalingModSize << ","
              << (p.enforce128BitSecurity ? 1 : 0) << "," << (cfg.modReduceAfterMult ? 1 : 0) << ","
              << cfg.numChains << "," << cfg.maxMults << "," << growthFile << "," << histFile << ",";
        uint32_t digits = p.numLargeDigits;
        std::string log2P;
        if (auto rns = std::dynamic_pointer_cast<CryptoParametersRNS>(params)) {
            digits = rns->GetNumPartQ();
            if (auto pp = rns->GetParamsP()) {
                std::ostringstream os;
                os << std::fixed << std::setprecision(1) << std::log2(pp->GetModulus().ConvertToDouble());
                log2P = os.str();
            }
        }
        index << digits << "," << log2P << "\n";
        index.flush();  // a partially finished sweep is still usable
    }
    std::cout << "\nSweep index written to " << indexPath << std::endl;
}