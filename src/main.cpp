// Includes from this library
#include "COSMO_SAC/COSMO.hpp"
using namespace COSMOSAC;

#include <string>
#include <numeric>
#include <map>
#include <regex>
#include <iostream>
#include <chrono>

#include <Eigen/Dense>
#include <chrono>

void speed_test() {

    constexpr int N = 8;
    Eigen::Array<double, N, N> A51;
    Eigen::Array<double, N, 1> b51;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> AXX(N,N);
    Eigen::Array<double, Eigen::Dynamic, 1> bX(N);
    A51.setRandom();
    AXX.setRandom();
    bX.setRandom();
    b51.setRandom();
    double s1 =0, s2 = 0;
    auto tic = std::chrono::high_resolution_clock::now();
    for (auto ii = 0; ii < 1000; ++ii){
        auto m51 = A51.matrix()*b51.matrix();
        s1 += A51.matrix().eigenvalues()(0).real();
    }
    auto toc = std::chrono::high_resolution_clock::now();
    for (auto ii = 0; ii < 1000; ++ii) {
        auto mxx = AXX.matrix()*bX.matrix();
        s2 += AXX.matrix().eigenvalues()(0).real();
    }
    auto toc2 = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double>(toc - tic).count()*1e6 << std::endl;
    std::cout << std::chrono::duration<double>(toc2 - toc).count()*1e6 << std::endl;
}

// The header HOME_PATH.h will be populated with a macro COSMO_SAC_HOME set to the value of the root of the repository, 
// such that the folder "profiles" will be inside it
#include "HOME_PATH.h"

void testVT2005() {
    // Construct the database (empty to start)
    VirginiaTechProfileDatabase datadb(
        COSMO_SAC_HOME + std::string("/profiles/VT2005/Sigma_Profile_Database_Index_v2.txt"),
        COSMO_SAC_HOME + std::string("/profiles/VT2005/Sigma_Profiles_v2/")
    );

    // Load the desired profiles into the database
    std::vector<std::string> names = { "N-PENTANE", "BENZENE" };
    for (auto &i : names) {
        auto ident = datadb.normalize_identifier(i);
        datadb.add_profile(ident);
        //std::cout << datadbDel.get_profile(ident).name << std::endl;
    }

    // Build the evaluator class
    COSMO1 COSMO(names, datadb);

    // The parameters at which we will evaluate
    double T = 298.15;
    Eigen::ArrayXd z(2); z.fill(0.0); z[1] = 1 - z[0];
    COSMO1::EigenArray51 psigma_mix = COSMO.get_psigma_mix(z);
    // Carry out the calculation(s)
    COSMO1::EigenArray51 Gammafast, Gammaslow; 
    
    double r1=0,r2=0;
    auto startTime = std::chrono::high_resolution_clock::now();
    for (auto ii = 0; ii < 100; ++ii){
        Gammafast = COSMO.get_Gamma_fast(T, psigma_mix); r1 += Gammafast(25); }
    auto endTimef = std::chrono::high_resolution_clock::now();
    for (auto ii = 0; ii < 100; ++ii) {
        Gammaslow = COSMO.get_Gamma_slow(T, psigma_mix); r2 += Gammaslow(25); }
    auto endTimes = std::chrono::high_resolution_clock::now();
    std::cout << r1 << " " << std::chrono::duration<double>(endTimef - startTime).count()/100 << " s elapsed\n";
    std::cout << r2 << " " << std::chrono::duration<double>(endTimes - endTimef).count()/100 << " s elapsed\n" ;
    /*std::cout << Gammafast << std::endl;
    std::cout << "----------------------------------\n";
    std::cout << Gammaslow << std::endl;*/
}
int main() {
    testVT2005();
    //speed_test(); speed_test(); speed_test(); speed_test();
    return EXIT_SUCCESS;

    // Construct the database (empty to start)
    DelawareProfileDatabase datadbDel(COSMO_SAC_HOME+std::string("/profiles/UD/complist.txt"),
                                      COSMO_SAC_HOME+std::string("/profiles/UD/test3/"));

    // Load the desired profiles into the database
    std::vector<std::string> names = { "N-PENTANE", "BENZENE" };
    for (auto &i : names) {
        auto ident = datadbDel.normalize_identifier(i);
        datadbDel.add_profile(ident);
        //std::cout << datadbDel.get_profile(ident).name << std::endl;
    }

    // Build the evaluator class
    COSMO3 COSMO(names, datadbDel);
    
    // The parameters at which we will evaluate
    double T = 298.15;
    Eigen::ArrayXd z(2); z.fill(0.0); z[1] = 1 - z[0];
    auto startTime = std::chrono::high_resolution_clock::now();
   
        // Carry out the calculation
        double res_comb = COSMO.get_lngamma_comb(T, z, 0);
        Eigen::ArrayXd res = COSMO.get_lngamma(T, z);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double>(endTime - startTime).count() << " s elapsed\n";
    std::cout << res << std::endl;
}