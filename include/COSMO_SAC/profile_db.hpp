#ifndef COSMO_PROFILES
#define COSMO_PROFILES

#include "COSMO_SAC/util.hpp"
#include <Eigen/Dense>
#include <map>
#include <limits>
#include <iterator>
#include <algorithm>
#include <cctype>
#include <functional>
#include "nlohmann/json.hpp"

namespace COSMOSAC {

/**
A single sigma profile.  In this data structure (for consistency with the 2005 Virginia Tech
database of COSMO-SAC parameters), the first column is the electron density in [e/A^2], and the second column
is the probability of finding a segment with this electron density, multiplied by the segment area, in A^2
*/
class SigmaProfile {
public:
    Eigen::ArrayXd m_sigma, m_psigmaA;
    /// Default constructor
    SigmaProfile() {};
    /// Copy-by-reference constructor
    SigmaProfile(const Eigen::ArrayXd &sigma, const Eigen::ArrayXd &psigmaA) : m_sigma(sigma), m_psigmaA(psigmaA) {};
    /// Move constructor
    SigmaProfile(Eigen::ArrayXd &&sigma, Eigen::ArrayXd &&psigmaA) : m_sigma(sigma), m_psigmaA(psigmaA) {};
    const Eigen::ArrayXd &psigmaA() const { return m_psigmaA; }
    const Eigen::ArrayXd &sigma() const { return m_sigma; }
    const Eigen::ArrayXd psigma(double A_i) const { return m_psigmaA/A_i; }

};

struct SigmaProfileSet {
    SigmaProfile nhb, ///< The profile for the non-hydrogen-bonding segments
        oh, ///< The profile for the OH-bonding segments
        ot; ///< The profile for the "other" segments
};

/// The sigma profiles (and some more metadata) associated with the fluid
struct FluidProfiles {
    enum class dispersion_classes{DISP_WATER, DISP_COOH, DISP_NHB, DISP_ONLY_ACCEPTOR, DISP_DONOR_ACCEPTOR};
    SigmaProfileSet profiles;
    std::string name;
    std::size_t VTnumber;
    double A_COSMO_A2; ///< The surface area of the molecule as calculated by COSMO-SAC, in A^2
    double V_COSMO_A3; ///< The volume of the molecule as calculated by COSMO-SAC, in A^3
    dispersion_classes dispersion_flag;
    double dispersion_eoverkB;
    nlohmann::json meta; ///< Any additional metadata, stored in JSON format
};

///
class ProfileDatabase {
private:
    std::map<std::string, FluidProfiles> datadb;
public:
    FluidProfiles get_profile(const std::string & identifier) const{
        if (datadb.find(identifier) != datadb.cend()) {
            return datadb.find(identifier)->second;
        }
        else {
            throw std::invalid_argument("Unable to load the profile for the identifier [" + identifier + "]");
        }
    }
    virtual std::string normalize_identifier(const std::string &identifier) const { return identifier; };
    std::size_t length(const std::string & identifier) {
        return datadb.size();
    }
    void clear() {
        datadb.clear();
    }
    void add_to_db(const std::string & identifier, FluidProfiles &&value) {
        datadb[identifier] = value;
    }
};

struct VTFluidInfo {
    int VTIndex;
    std::string numstring;
    std::string name;
    std::string CAS;
    double V_COSMO_A3;
};
void to_json(nlohmann::json& j, const VTFluidInfo& p) {
    j = nlohmann::json{{"VTIndex", p.VTIndex}, {"numstring", p.numstring}, {"name", p.name}, {"CAS", p.CAS},  {"V_COSMO_A3", p.V_COSMO_A3}};
}
void from_json(const nlohmann::json& j, VTFluidInfo& p) {
    j.at("VTIndex").get_to(p.VTIndex);
    j.at("numstring").get_to(p.numstring);
    j.at("name").get_to(p.name);
    j.at("CAS").get_to(p.CAS);
    j.at("V_COSMO_A3").get_to(p.V_COSMO_A3);
}

/// The format of the VirginiaTech database of profiles
class VirginiaTechProfileDatabase : public ProfileDatabase{
private:
    
    std::map<std::string, VTFluidInfo> m_VTdata;
    const std::string m_dbpath;
public:    
    VirginiaTechProfileDatabase(const std::string &tab_delimited_summary, const std::string &dbpath) : m_dbpath(dbpath) {
        // First, read in the file of summary information
        auto lines = str_split(get_file_contents(tab_delimited_summary), "\n");
        if (lines.size() == 1){
            lines = str_split(get_file_contents(tab_delimited_summary), "\r");
        }
        if (lines.size() == 1){
            throw std::invalid_argument("Unable to match the line endings in file " + tab_delimited_summary);
        }
        // Iterate over the lines, storing an entry for each one
        for (auto &&line : lines) {
            VTFluidInfo fluid;
            if (line.empty()) { continue; }
            auto elements = str_split(strstrip(line), "\t");
            if (elements.size() < 6) {
                throw std::invalid_argument("Need at least 6 entries in line");
            }
            if (elements[0] == "Index No.") {
                continue;
            }
            fluid.VTIndex = mystrtoi(elements[0]);
            // Index 1 is a chemical formula, but non-standard form
            fluid.name = elements[2];
            fluid.CAS = elements[3];
            fluid.V_COSMO_A3 = mystrtod(elements[5]);
            char num[5];
            snprintf(num, sizeof(num), "%04d", static_cast<int>(fluid.VTIndex));
            fluid.numstring = num;
            if (m_VTdata.find(fluid.numstring) != m_VTdata.end()) {
                throw std::invalid_argument("Cannot have duplicate key: "+ fluid.numstring);
            }
            else{
                m_VTdata[fluid.numstring] = fluid;
            }
            
        }
    }
    std::string normalize_identifier(const std::string &identifier) const override {
        // The key is in the database, return it
        if (m_VTdata.find(identifier) != m_VTdata.end()) {
            return identifier;
        }
        for (auto &el : m_VTdata) {
            auto &fluid = el.second;
            if (fluid.name == identifier)
            { 
                // Construct the fluid identifier from the identifier that is provided
                char num[5];
                snprintf(num, sizeof(num), "%04d", static_cast<int>(fluid.VTIndex));
                return num;
            }
            if (fluid.CAS == identifier)
            { 
                // Construct the fluid identifier from the identifier that is provided
                char num[5];
                snprintf(num, sizeof(num), "%04d", static_cast<int>(fluid.VTIndex));
                return num;
            }
        }
        throw std::invalid_argument("Unable to match identifier in secondary lookup");
    }
    void add_profile(std::string identifier) {
        auto info = m_VTdata[identifier];
        // Now we load the sigma profile(s) from the file
        auto lines = str_split(get_file_contents(m_dbpath + "/VT2005-" + info.numstring + "-PROF.txt"));
        std::vector<double> sigma, psigmaA;
        FluidProfiles fluid;
        for (auto &&line : lines) {
            if (line.size() > 10) {
                auto _sigma = mystrtod(line.substr(0, 25));
                auto _psigmaA = mystrtod(line.substr(25, 25));
                sigma.push_back(_sigma);
                psigmaA.push_back(_psigmaA);
            }
        }
        fluid.name = info.name;
        fluid.VTnumber = info.VTIndex;
        fluid.V_COSMO_A3 = info.V_COSMO_A3;
        if (sigma.size() == 51 && psigmaA.size() == 51) {
            fluid.profiles.nhb = SigmaProfile(Eigen::Map<Eigen::ArrayXd>(&(sigma[0]), sigma.size()),
                Eigen::Map<Eigen::ArrayXd>(&(psigmaA[0]), psigmaA.size()));
            fluid.A_COSMO_A2 = fluid.profiles.nhb.psigmaA().sum();
        }
        else {
            throw std::invalid_argument("Don't support 2 & 3 sigma profiles yet");
        }
        add_to_db(info.numstring, std::move(fluid));
    }
    std::string to_JSON(){
        auto nspaces = 2;
        nlohmann::json j = m_VTdata;
        return j.dump(nspaces);
    }
};

struct DelawareFluidInfo {
    std::size_t DelIndex;
    std::string CAS;
    std::string formula;
    std::string name;
    std::string SMILES;
    std::string InChIString; ///< Standard InChI string
    std::string InChIKey; ///< Standard InChI key
    double V_COSMO_A3;
    nlohmann::json meta;
};
void to_json(nlohmann::json& j, const DelawareFluidInfo& p) {
    j = nlohmann::json{
        {"DelIndex", p.DelIndex}, {"CAS", p.CAS}, {"formula", p.formula},  {"name", p.name},
        {"SMILES", p.SMILES}, {"InChIString", p.InChIString}, {"InChIKey", p.InChIKey},  {"V_COSMO_A3", p.V_COSMO_A3}
};
}
void from_json(const nlohmann::json& j, DelawareFluidInfo& p) {
    j.at("DelIndex").get_to(p.DelIndex);
    j.at("CAS").get_to(p.CAS);
    j.at("formula").get_to(p.formula);
    j.at("name").get_to(p.name);
    j.at("SMILES").get_to(p.SMILES);
    j.at("InChIString").get_to(p.InChIString);
    j.at("InChIKey").get_to(p.InChIKey);
    j.at("V_COSMO_A3").get_to(p.V_COSMO_A3);
    j.at("name").get_to(p.name);
}
    
/// The format of the Delaware database of profiles
class DelawareProfileDatabase : public ProfileDatabase {
private:
    
    std::map<std::string, DelawareFluidInfo> m_Deldata;
    const std::string m_dbpath;
public:
    DelawareProfileDatabase(const std::string &space_delimited_summary, const std::string &dbpath) : m_dbpath(dbpath) {
        // First, read in the file of summary information
        auto lines = str_split(get_file_contents(space_delimited_summary), "\n");
        // Iterate over the lines, storing an entry for each one
        for (auto &&line : lines) {
            DelawareFluidInfo fluid;
            if (line.empty()) { continue; }
            auto elements = str_split(strrstrip(line), " ");
            if (elements[0] == "ID") {
                continue;
            }
            fluid.DelIndex = mystrtoi(elements[0]);
            fluid.formula = elements[1];
            fluid.CAS = elements[2];
            fluid.name = elements[3];
            fluid.SMILES = elements[4];
            fluid.InChIString = elements[5];
            fluid.InChIKey = elements[6];
            m_Deldata[fluid.InChIKey] = fluid;
        }
    }
    std::string normalize_identifier(const std::string &identifier) const override {
        // The key is in the database, return it
        if (m_Deldata.find(identifier) != m_Deldata.end()) {
            return identifier;
        }
        for (auto &el : m_Deldata) {
            auto &fluid = el.second;
            if (fluid.name == identifier)
            {
                return fluid.InChIKey;
            }
            if (fluid.CAS == identifier)
            {
                return fluid.InChIKey;
            }
        }
        throw std::invalid_argument("Unable to match identifier in secondary lookup");
    }
    void add_profile(std::string identifier) {
        auto info = m_Deldata[identifier];
        // Now we load the sigma profile(s) from the file
        auto lines = str_split(get_file_contents(m_dbpath + "/" + info.InChIKey + ".sigma"));

        std::vector<double> _sigma, _psigmaA;
        FluidProfiles fluid;
        std::string meta;
        for (auto &&line : lines) {
            if (line.substr(0,8) == "# Name: "){ std::string check_name = line.substr(8,line.size()-8); continue;}
            if (line.substr(0,8) == "# CASn: ") { std::string check_CAS = line.substr(8, line.size() - 8); continue; }
            if (line.substr(0, 8) == "# meta: ") { meta = line.substr(8, line.size() - 8); continue; }
            if (line[0] == '#'){ continue; }
            auto v = str_split(strrstrip(line), " ");
            if (v.empty() || (v.size() ==1 && v[0].empty())){ continue; }
            _sigma.push_back(mystrtod(v[0]));
            _psigmaA.push_back(mystrtod(v[1]));
        }
        Eigen::Map<Eigen::ArrayXd> sigma(&(_sigma[0]), _sigma.size());
        Eigen::Map<Eigen::ArrayXd> psigmaA(&(_psigmaA[0]), _psigmaA.size());
        fluid.name = info.name;
        
        fluid.VTnumber = info.DelIndex;
        if (sigma.size() == 51 && psigmaA.size() == 51) {
            fluid.profiles.nhb = SigmaProfile(sigma, psigmaA);
            fluid.A_COSMO_A2 = fluid.profiles.nhb.psigmaA().sum();
        }
        else if (sigma.size() == 51*3 && psigmaA.size() == 51*3){
            fluid.profiles.nhb = SigmaProfile(sigma.segment(0*51, 51), psigmaA.segment(0*51, 51));
            fluid.profiles.oh =  SigmaProfile(sigma.segment(1*51, 51), psigmaA.segment(1*51, 51));
            fluid.profiles.ot =  SigmaProfile(sigma.segment(2*51, 51), psigmaA.segment(2*51, 51));
            double check_Area2 = fluid.profiles.nhb.psigmaA().sum() + fluid.profiles.oh.psigmaA().sum() + fluid.profiles.ot.psigmaA().sum();
            fluid.A_COSMO_A2 = check_Area2;
        }
        else{
            throw std::invalid_argument("Length of sigma profile ["+std::to_string(sigma.size())+"] is neither 51 nor 51*3");
        }
        fluid.meta = nlohmann::json::parse(meta);
        fluid.V_COSMO_A3 = fluid.meta["volume [A^3]"];
        std::string flag = fluid.meta["disp. flag"];
        if (flag == "COOH") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_COOH;
        }
        else if (flag == "H2O") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_WATER;
        }
        else if (flag == "NHB") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_NHB;
        }
        else if (flag == "HB-ACCEPTOR") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_ONLY_ACCEPTOR;
        }
        else if (flag == "HB-DONOR-ACCEPTOR") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_DONOR_ACCEPTOR;
        }
        else {
            throw std::invalid_argument("Unable to match dispersion flag: \""+flag+"\"");
        }
        if (fluid.meta["disp. e/kB [K]"].is_null()){
            // The null from JSON is mapped to a NaN of C++
            fluid.dispersion_eoverkB = std::numeric_limits<double>::quiet_NaN();
        }
        else {
            fluid.dispersion_eoverkB = fluid.meta["disp. e/kB [K]"];
        }
        add_to_db(info.InChIKey, std::move(fluid));
    }
    std::string to_JSON(){
        auto nspaces = 2;
        nlohmann::json j = m_Deldata;
        return j.dump(nspaces);
    }
};

/// The DirectImport feature for importing profiles based on "filename" + "path" 
class DirectImport : public ProfileDatabase {
private:
    std::string m_dbpath;

public:
    DirectImport() {
    }

    void add_profile(const std::string& identifier, const std::string& dbpath = ".") {
        m_dbpath = dbpath;
        // Now we load the sigma profile(s) from the file
        auto lines = str_split(get_file_contents(m_dbpath + "/" + identifier + ".sigma"));

        std::vector<double> _sigma, _psigmaA;
        FluidProfiles fluid;
        std::string meta;
        for (auto &&line : lines) {
            if (line.substr(0,8) == "# Name: "){ std::string check_name = line.substr(8,line.size()-8); continue;}
            if (line.substr(0,8) == "# CASn: ") { std::string check_CAS = line.substr(8, line.size() - 8); continue; }
            if (line.substr(0, 8) == "# meta: ") { meta = line.substr(8, line.size() - 8); continue; }
            if (line[0] == '#'){ continue; }
            auto v = str_split(strrstrip(line), " ");
            if (v.empty() || (v.size() ==1 && v[0].empty())){ continue; }
            _sigma.push_back(mystrtod(v[0]));
            _psigmaA.push_back(mystrtod(v[1]));
        }
        Eigen::Map<Eigen::ArrayXd> sigma(&(_sigma[0]), _sigma.size());
        Eigen::Map<Eigen::ArrayXd> psigmaA(&(_psigmaA[0]), _psigmaA.size());
        fluid.name = identifier;

        if (sigma.size() == 51 && psigmaA.size() == 51) {
            fluid.profiles.nhb = SigmaProfile(sigma, psigmaA);
            fluid.A_COSMO_A2 = fluid.profiles.nhb.psigmaA().sum();
        }
        else if (sigma.size() == 51*3 && psigmaA.size() == 51*3){
            fluid.profiles.nhb = SigmaProfile(sigma.segment(0*51, 51), psigmaA.segment(0*51, 51));
            fluid.profiles.oh =  SigmaProfile(sigma.segment(1*51, 51), psigmaA.segment(1*51, 51));
            fluid.profiles.ot =  SigmaProfile(sigma.segment(2*51, 51), psigmaA.segment(2*51, 51));
            double check_Area2 = fluid.profiles.nhb.psigmaA().sum() + fluid.profiles.oh.psigmaA().sum() + fluid.profiles.ot.psigmaA().sum();
            fluid.A_COSMO_A2 = check_Area2;
        }
        else{
            throw std::invalid_argument("Length of sigma profile ["+std::to_string(sigma.size())+"] is neither 51 nor 51*3");
        }
        fluid.meta = nlohmann::json::parse(meta);
        fluid.V_COSMO_A3 = fluid.meta["volume [A^3]"];
        std::string flag = fluid.meta["disp. flag"];
        if (flag == "COOH") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_COOH;
        }
        else if (flag == "H2O") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_WATER;
        }
        else if (flag == "NHB") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_NHB;
        }
        else if (flag == "HB-ACCEPTOR") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_ONLY_ACCEPTOR;
        }
        else if (flag == "HB-DONOR-ACCEPTOR") {
            fluid.dispersion_flag = FluidProfiles::dispersion_classes::DISP_DONOR_ACCEPTOR;
        }
        else {
            throw std::invalid_argument("Unable to match dispersion flag: \""+flag+"\"");
        }
        if (fluid.meta["disp. e/kB [K]"].is_null()){
            // The null from JSON is mapped to a NaN of C++
            fluid.dispersion_eoverkB = std::numeric_limits<double>::quiet_NaN();
        }
        else {
            fluid.dispersion_eoverkB = fluid.meta["disp. e/kB [K]"];
        }
        add_to_db(identifier, std::move(fluid));
    }
};

} /* namespace COSMOSAC */

#endif
