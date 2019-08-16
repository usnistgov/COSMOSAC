// Includes from this library
#include "COSMO_SAC/COSMO.hpp"
using namespace COSMOSAC;
#include "HOME_PATH.h"

#include "catch2/catch.hpp"

TEST_CASE("Test infinite dilution activity coefficients", "[setup]") {

    struct el {
        std::string solute, solvent;
        double T, lngammai_inf_dil_comb, lngammai_inf_dil_res, lngammai_inf_dil_dsp;
    };

    /// From Hsieh, Fluid Phase Equil., 2014, Table 6
    std::vector<el> data = {
        { "N-PENTANE","N-HEPTANE",293.15, -0.0317, 0.0001, 0.0000 },
        { "N-PENTANE","N-TETRACOSANE", 298.15, -0.6123, 0.0012, 0.0000 },
        { "N-PENTANE","BENZENE", 298.15, -0.0075,0.8001, 0.0009 },
        { "N-PENTANE","1,2-DICHLOROETHANE", 298.15, -0.0146, 1.2032, 0.0098 },
        { "1,1-DICHLOROETHANE","1,1,1-TRICHLOROETHANE",328.2, -0.0032, 0.0787, -0.0004 },
        { "N-PENTANE","ACETONITRILE",298.15,-0.1610,3.3610, 0.1095 },
        { "N-PENTANE","ACETONE",298.15,-0.0513,1.1924, 0.3317 },
        { "BENZENE","TRIETHYLAMINE",323.5,-0.0413,0.2426, 0.0799 },
        { "1-PENTENE","ETHYL_BENZOATE",313.2,-0.0632,0.1426, 0.0514 },
        { "ISOPRENE","METHYL_ETHYL_KETONE",293.15,0.0003,-0.0206, 0.2277 },
        { "TOLUENE","ETHYL_ACETATE",298.2,0.0094,-0.0185, 0.2145 },
        { "CHLOROFORM","N-HEXANE",301,-0.0434,0.3346, 0.0224 },
        { "ACETONE","1-BUTANOL",298,-0.0255,0.5548, 0.0525 },
        { "N-PENTANE","ETHANOL",298,-0.1214,1.8569, 0.2845 },
        { "BENZENE","WATER",298,-1.0519,8.4222, 0.7955 },
        { "N-PENTANE","WATER",298.15,-1.3733,12.3241, 0.7429 },
        { "1-OCTADECANOL","WATER",298.15,-7.8329,29.9470, 0.5793 },
        { "N-HEPTANE","ACETIC_ACID",298,-0.2835,4.3649, -0.8758 }
    };
    
    for (auto &el : data) {
        SECTION("Inf. dilution of " + el.solute + " in " + el.solvent){

        // Construct the database (empty to start)
        DelawareProfileDatabase datadbDel(
            COSMO_SAC_HOME + std::string("/profiles/UD/complist.txt"),
            COSMO_SAC_HOME + std::string("/profiles/UD/sigma3/"));

        // Load the desired profiles into the database
        std::vector<std::string> names = { el.solute, el.solvent };
        CAPTURE(el.solvent);
        CAPTURE(el.solute);
        CAPTURE(el.T);
        CAPTURE(el.lngammai_inf_dil_comb);
        CAPTURE(el.lngammai_inf_dil_res);

        std::string ident_solute = datadbDel.normalize_identifier(names[0]); 
        std::string ident_solvent = datadbDel.normalize_identifier(names[1]);
        
        CAPTURE(ident_solvent); 
        CAPTURE(ident_solute);
        for (auto &i : names) {
            std::string ident;
            CHECK_NOTHROW(ident = datadbDel.normalize_identifier(i));
            CHECK_NOTHROW(datadbDel.add_profile(ident));
        }

        // Build the evaluator class
        COSMO3<double> COSMO(names, datadbDel);
        Eigen::ArrayXd z(2); z.fill(0.0); z[1] = 1 - z[0]; 
        
        // Combinatorial part
        double lngammai_comb = COSMO.get_lngamma_comb(el.T, z, 0);
        double comb_err = std::abs(lngammai_comb - el.lngammai_inf_dil_comb);
        CHECK(comb_err < 2e-4); 
        
        // Residual part
        double lngammai = COSMO.get_lngamma(el.T, z)(0);
        double res_err = std::abs(lngammai - (el.lngammai_inf_dil_comb + el.lngammai_inf_dil_res + el.lngammai_inf_dil_dsp));
        CHECK(res_err < 2e-3);

        // Dispersive part
        double lngammai_dsp = COSMO.get_lngamma_disp(z)(0);
        CAPTURE(lngammai_dsp);
        double dsp_err = std::abs(lngammai_dsp - (el.lngammai_inf_dil_dsp));
        CHECK(dsp_err < 2e-3);

        }
    }
    
}

TEST_CASE("Test slow and fast Gamma are the same for COSMO1", "[Gamma]") {

    // Construct the database (empty to start)
    VirginiaTechProfileDatabase datadb(
        COSMO_SAC_HOME + std::string("/profiles/VT2005/Sigma_Profile_Database_Index_v2.txt"),
        COSMO_SAC_HOME + std::string("/profiles/VT2005/Sigma_Profiles_v2/")
    );

    // Load the desired profiles into the database
    std::vector<std::string> names = { "N-PENTANE", "ETHANOL" };
    for (auto &i : names) {
        auto ident = datadb.normalize_identifier(i);
        datadb.add_profile(ident);
        //std::cout << datadbDel.get_profile(ident).name << std::endl;
    }

    // Build the evaluator class
    COSMO1<double> COSMO(names, datadb);

    // The parameters at which we will evaluate
    double T = 298.15;
    Eigen::ArrayXd z(2); z.fill(0.0); z[1] = 1 - z[0];

    // Carry out the calculation(s)
    std::vector<COSMO1<double>::EigenArray> outputs;
    for (auto fast : {true, false}){
        COSMO.get_mutable_COSMO_constants().fast_Gamma = fast;
        outputs.push_back(COSMO.get_lngamma(T, z));
    }
    double diff = (outputs[0] - outputs[1]).cwiseAbs().sum();
    CHECK(diff < 1e-7);
}

TEST_CASE("Test slow and fast Gamma are the same for COSMO3", "[Gamma]") {

    // Construct the database (empty to start)
    DelawareProfileDatabase datadbDel(
        COSMO_SAC_HOME + std::string("/profiles/UD/complist.txt"),
        COSMO_SAC_HOME + std::string("/profiles/UD/sigma3/"));

    // Load the desired profiles into the database
    std::vector<std::string> names = { "WATER", "ETHANOL" };
    for (auto &i : names) {
        auto ident = datadbDel.normalize_identifier(i);
        datadbDel.add_profile(ident);
        //std::cout << datadbDel.get_profile(ident).name << std::endl;
    }

    // Build the evaluator class
    COSMO3<double> COSMO(names, datadbDel);

    // The parameters at which we will evaluate
    double T = 298.15;
    Eigen::ArrayXd z(2); z.fill(0.0); z[1] = 1 - z[0];

    // Carry out the calculation(s)
    std::vector<COSMO1<double>::EigenArray> outputs;
    for (auto fast : { true, false }) {
        COSMO.get_mutable_COSMO_constants().fast_Gamma = fast;
        outputs.push_back(COSMO.get_lngamma(T, z));
    }
    double diff = (outputs[0] - outputs[1]).cwiseAbs().sum();
    CHECK(diff < 1e-8);
}

TEST_CASE("Test slow and fast DELTAW are the same", "[Gamma]") {

    // Construct the database (empty to start)
    DelawareProfileDatabase datadbDel(
        COSMO_SAC_HOME + std::string("/profiles/UD/complist.txt"),
        COSMO_SAC_HOME + std::string("/profiles/UD/sigma3/"));

    // Load the desired profiles into the database
    std::vector<std::string> names = { "N-PENTANE", "ETHANOL" };
    for (auto &i : names) {
        auto ident = datadbDel.normalize_identifier(i);
        datadbDel.add_profile(ident);
        //std::cout << datadbDel.get_profile(ident).name << std::endl;
    }

    // Build the evaluator class
    COSMO3<double> COSMO(names, datadbDel);
    double T = 298.15; 
    using profile_type = COSMOSAC::AbstractCOSMOModel<double>::profile_type;
    
    Eigen::Index ileft, w;
    std::tie(ileft, w) = COSMO.get_ileftw();
    std::vector<profile_type> types = { profile_type::NHB_PROFILE, profile_type::OH_PROFILE, profile_type::OT_PROFILE };
    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            COSMO3<double>::EigenMatrix diff = COSMO.get_DELTAW(T, types[i], types[j]).array() - COSMO.get_DELTAW_fast(T, types[i], types[j]).array();
            CHECK(diff.matrix().block(ileft, ileft, w, w).cwiseAbs().maxCoeff() < 1e-15);
        }
    }
}

//TEST_CASE("Test fast Gamma is faster", "[Gamma]") {
//
//    // Construct the database (empty to start)
//    DelawareProfileDatabase datadbDel(
//        COSMO_SAC_HOME + std::string("/profiles/UD/complist.txt"),
//        COSMO_SAC_HOME + std::string("/profiles/UD/sigma3/"));
//
//    // Load the desired profiles into the database
//    std::vector<std::string> names = { "CARBON_DIOXIDE", "WATER" };
//    for (auto &i : names) {
//        auto ident = datadbDel.normalize_identifier(i);
//        datadbDel.add_profile(ident);
//        //std::cout << datadbDel.get_profile(ident).name << std::endl;
//    }
//
//    // Build the evaluator class
//    COSMO3<double> COSMO(names, datadbDel);
//
//    // The parameters at which we will evaluate
//    double T = 298.15;
//    Eigen::ArrayXd z(2); z[0] = 0.9975; z[1] = 1 - z[0];
//
//    // Carry out the calculation(s)
//    std::vector<std::tuple<double, COSMO1<double>::EigenArray> > outputs;
//    for (auto fast : { true, false }) {
//        COSMO.get_mutable_COSMO_constants().fast_Gamma = fast;
//        auto startTime = std::chrono::high_resolution_clock::now();
//        for (auto i = 0; i < 1000; ++i) {
//            auto lngamma = COSMO.get_lngamma(T, z);
//        }
//        auto endTime = std::chrono::high_resolution_clock::now();
//        auto elapsed_s = std::chrono::duration<double>(endTime - startTime).count()/1000;
//        outputs.push_back(std::make_tuple(elapsed_s, COSMO.get_lngamma(T, z)));
//    }
//    double diff = (std::get<1>(outputs[0]) - std::get<1>(outputs[1])).cwiseAbs().sum();
//    double speedup = std::get<0>(outputs[1]) / std::get<0>(outputs[0]);
//    std::cout << std::get<0>(outputs[1]) << std::endl;
//    std::cout << std::get<0>(outputs[0]) << std::endl;
//    CAPTURE(speedup);
//    CHECK(diff < 1e-8);
//    CHECK(speedup > 1);
//}