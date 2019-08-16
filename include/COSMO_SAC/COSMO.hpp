#ifndef COSMO_COSMO
#define COSMO_COSMO

#include "COSMO_SAC/profile_db.hpp"

#include <chrono>

namespace COSMOSAC {
    template<typename TYPE> TYPE POW2(TYPE x){return x*x;}

    struct CombinatorialConstants {
        double q0 = 79.53, // [A^2]
            r0 = 66.69, // [A^3]
            z_coordination = 10.0;
        std::string to_string() {
            return "q0: " + std::to_string(q0) + " A^2 \nr0: " + std::to_string(r0) + " A^3\nz_coordination: " + std::to_string(z_coordination);
        }
    };

    class AbstractCOSMOModel {
    public:
        typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
        typedef Eigen::Array<double, Eigen::Dynamic, 1> EigenArray;
        enum class profile_type { NHB_PROFILE, OH_PROFILE, OT_PROFILE };
    private:
        const ProfileDatabase & datadb;
    protected:
        /// A mini convenience function to avoid verbose type casting
        static double dot(const EigenArray &v1, const EigenArray &v2) {
            return v2.matrix().dot(v1.matrix());
        }
        std::vector<FluidProfiles> build_fluids(const std::vector<std::string> &identifiers) {
            std::vector<FluidProfiles> fluids;
            for (auto & i : identifiers) {
                auto ident = datadb.normalize_identifier(i);
                // See https://stackoverflow.com/a/15531940/1360263 for explanation of why there is a "this" here
                fluids.push_back(this->get_sigma_profiles(ident));
            }
            return fluids;
        }
        std::vector<FluidProfiles> m_fluids;
        CombinatorialConstants m_comb_consts;
    public:
        AbstractCOSMOModel(const ProfileDatabase &db, const std::vector<std::string> &identifiers) : datadb(db), m_fluids(build_fluids(identifiers)) {};

        std::tuple<Eigen::Index, Eigen::Index> get_nonzero_bounds(){
            // Determine the range of entries in p(sigma) that are greater than zero, we 
            // will only calculate segment activity coefficients for these segments
            Eigen::Index min_ileft = 51, max_iright = 0;
            for (auto i = 0; i < m_fluids.size(); ++i) {
                const EigenArray psigma = m_fluids[i].profiles.nhb.psigma(m_fluids[i].A_COSMO_A2);
                Eigen::Index ileft = 0, iright = psigma.size();
                for (auto ii = 0; ii < psigma.size(); ++ii) { if (std::abs(psigma(ii) > 1e-16)) { ileft = ii; break; } }
                for (auto ii = psigma.size() - 1; ii > ileft; --ii) { if (std::abs(psigma(ii) > 1e-16)) { iright = ii; break; } }
                if (ileft < min_ileft) { min_ileft = ileft; }
                if (iright > max_iright) { max_iright = iright; }
            }
            return std::make_tuple(min_ileft, max_iright);
        }

        CombinatorialConstants &get_mutable_combinatorial_constants(){ return m_comb_consts; }
        /**
        The combinatorial part of ln(γ_i)
        */
        double get_lngamma_comb(double T, const EigenArray &x, std::size_t i) const {

            EigenArray A_A2(x.size()), ///< The surface area of each compound in A^2
                       V_A3(x.size()); ///< The volume of each compound in A^3
            for (auto i = 0; i < x.size(); ++i) {
                A_A2(i) = m_fluids[i].A_COSMO_A2;
                V_A3(i) = m_fluids[i].V_COSMO_A3;
            }

            double q0 = m_comb_consts.q0, 
                  r0 = m_comb_consts.r0, 
                  z_coordination = m_comb_consts.z_coordination;
            EigenArray q = A_A2 / q0,
                r = V_A3 / r0,
                l = z_coordination / 2 * (r - q) - (r - 1);
            double theta_i = x[i] * q[i] / dot(x, q),
                phi_i = x[i] * r[i] / dot(x, r),
                phi_i_over_x_i = r[i] / dot(x, r),
                theta_i_over_phi_i = (q[i] / dot(x, q))/(r[i] / dot(x, r));
            return (log(phi_i_over_x_i) + z_coordination / 2 * q[i] * log(theta_i_over_phi_i)
                + l[i] - phi_i_over_x_i * dot(x, l));
        }

        EigenArray get_psigma_mix(const EigenArray &z, profile_type type = profile_type::NHB_PROFILE) const {
            EigenArray psigma_mix(51); psigma_mix.fill(0);
            double xA = 0;
            for (auto i = 0; i < z.size(); ++i) {
                switch (type) {
                case profile_type::NHB_PROFILE:
                    psigma_mix += z[i] * m_fluids[i].profiles.nhb.psigmaA(); break;
                case profile_type::OH_PROFILE:
                    psigma_mix += z[i] * m_fluids[i].profiles.oh.psigmaA(); break;
                case profile_type::OT_PROFILE:
                    psigma_mix += z[i] * m_fluids[i].profiles.ot.psigmaA(); break;
                }                
                xA += z[i] * m_fluids[i].A_COSMO_A2;
            }
            psigma_mix /= xA;
            return psigma_mix;
        }

        /*
        Get the profile(s)
        If 51 elements, there is just one profile
        If 153 elements, there are 3 profiles, one for non-hydrogen bonding (nhb),
        one for hydrogen bonding (hb), and one for other (oth)
        */
        FluidProfiles get_sigma_profiles(const std::string &identifier) const {
            return datadb.get_profile(identifier);
        };

        virtual EigenArray get_lngamma(double T, const EigenArray &x) const = 0;
    };

    struct COSMO1Constants {
        double c_hb = 85580, // [kcal A^4 /(mol*e^2)]
        sigma_hb = 0.0084, // [e/A^2]
        alpha_prime = 16466.72, // [kcal A^4 /(mol*e^2)]
        AEFFPRIME = 7.5, // [A^2]
        R = 0.001987; // [kcal/(mol*K)]; Universal gas constant divided by 4184, and truncated.  Consistent with Mullins, not consistent with CODATA
        bool fast_Gamma = false; ///< True to use the accelerated version of Gamma that only evaluates the entries with p(sigma) > 0 for all components
        std::string to_string() {
            return "c_hb: " + std::to_string(c_hb) + " kcal A^4 /(mol*e^2) \nsigma_hb: " + std::to_string(sigma_hb) + " e/A^2\nalpha_prime: " + std::to_string(alpha_prime) + " kcal A^4 /(mol*e^2)\nAEFFPRIME: "+std::to_string(AEFFPRIME)+" A\nR: "+std::to_string(R) + " kcal/(mol*K)\nfast_Gamma: "+ std::to_string(fast_Gamma);
        }
    };

    class COSMO1 : public AbstractCOSMOModel {
    public:
        typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
        typedef Eigen::Array<double, 51, 51> EigenMatrix5151;
        typedef Eigen::Array<double, Eigen::Dynamic, 1> EigenArray;
        typedef Eigen::Array<double, 51, 1> EigenArray51;
    private:
        
        EigenMatrix5151 build_DELTAW() {
            double sigma_hb = m_consts.sigma_hb,
                 alpha_prime = m_consts.alpha_prime,
                 c_hb = m_consts.c_hb;
            auto delta_sigma = (0.025 - (-0.025)) / 50;
            EigenMatrix5151 DELTAW;
            for (auto m = 0; m < 51; ++m) {
                for (auto n = 0; n < 51; ++n) {
                    double sigma_m = -0.025 + delta_sigma*m,
                        sigma_n = -0.025 + delta_sigma*n,
                        sigma_acc = std::max(sigma_m, sigma_n),
                        sigma_don = std::min(sigma_m, sigma_n);
                    DELTAW(m, n) = alpha_prime / 2 * (sigma_m + sigma_n)*(sigma_m + sigma_n) + c_hb*std::max(0.0, sigma_acc - sigma_hb)*std::min(0.0, sigma_don + sigma_hb);
                }
            }
            return DELTAW;
        }
        COSMO1Constants m_consts;
        EigenMatrix5151 m_DELTAW;
        Eigen::Index ileft, w;
    public:
        COSMO1(const std::vector<std::string> &identifiers, const ProfileDatabase &datadb, const COSMO1Constants &constants = COSMO1Constants())
            : AbstractCOSMOModel(datadb, identifiers), m_consts(constants), m_DELTAW(build_DELTAW()) {
            
            Eigen::Index iL, iR;
            std::tie(iL, iR) = get_nonzero_bounds();
            this->ileft = iL; this->w = iR-iL+1;
        };

        /// Get access to the parameters that are in use
        COSMO1Constants &get_mutable_COSMO_constants(){ return m_consts; }

        EigenMatrix5151 get_DELTAW() const {
            return m_DELTAW;
        }
        EigenArray51 get_Gamma_fast(double T, const EigenArray51 &psigma) const {
            auto tic = std::chrono::high_resolution_clock::now();
            double R = m_consts.R;
            EigenArray51 Gamma, Gammanew; Gamma.setOnes();
            EigenMatrix5151 AA;
            AA.matrix().block(ileft, ileft, w, w) = Eigen::exp(-m_DELTAW.matrix().block(ileft, ileft, w, w).array() / (R*T)).rowwise()*psigma.matrix().segment(ileft, w).array().transpose();
            
            auto toc = std::chrono::high_resolution_clock::now();
            // The iterative loop for Gamma
            for (auto counter = 0; counter < 200; ++counter) {
                Gammanew.segment(ileft,w) = 1/(AA.matrix().block(ileft,ileft,w,w).array().rowwise()*Gamma.matrix().segment(ileft,w).array().transpose()).rowwise().sum();
                Gamma.segment(ileft, w) = (Gamma.segment(ileft, w) + Gammanew.segment(ileft, w)) / 2;
                double maxdiff = ((Gamma.segment(ileft, w) - Gammanew.segment(ileft, w)) / Gamma.segment(ileft, w)).cwiseAbs().real().maxCoeff();
                if (maxdiff < 1e-8) {
                    break;
                }
            }
            auto toc2 = std::chrono::high_resolution_clock::now();
            //std::cout << std::chrono::duration<double>(toc - tic).count()*1e6 << std::endl;
            //std::cout << std::chrono::duration<double>(toc2 - toc).count()*1e6 << std::endl;
            return Gamma;
        }
        EigenArray51 get_Gamma_slow(const double T, const EigenArray51 &psigma) const {
            double R = m_consts.R;
            EigenArray51 Gamma, Gammanew; Gamma.setOnes();
            EigenMatrix5151 AA = Eigen::exp(-m_DELTAW / (R*T)).rowwise()*psigma.transpose();
            for (auto counter = 0; counter < 200; ++counter) {
                Gammanew = 1 / (AA.rowwise()*Gamma.transpose()).rowwise().sum();
                Gamma = (Gamma + Gammanew) / 2;
                double maxdiff = ((Gamma - Gammanew) / Gamma).cwiseAbs().real().maxCoeff();
                if (maxdiff < 1e-8) {
                    break;
                }
            }
            return Gamma;
        }
        EigenArray51 get_Gamma(double T, const EigenArray51 &psigma) const {
            return (m_consts.fast_Gamma) ? get_Gamma_fast(T, psigma) : get_Gamma_slow(T, psigma);
        }

        /**
        The residual part of ln(γ_i)
        */
        double get_lngamma_resid(std::size_t i, double T, const EigenArray &lnGamma_mix) const
        {
            double AEFFPRIME = m_consts.AEFFPRIME;
            const EigenArray psigma = m_fluids[i].profiles.nhb.psigma(m_fluids[i].A_COSMO_A2);
            const EigenArray lnGammai = get_Gamma(T, psigma).log();
            if (m_consts.fast_Gamma) {
                return m_fluids[i].A_COSMO_A2 / AEFFPRIME*(psigma.segment(ileft,w).array()*(lnGamma_mix.segment(ileft,w).array() - lnGammai.segment(ileft, w).array())).sum();
            }
            else{
                return m_fluids[i].A_COSMO_A2 / AEFFPRIME*(psigma*(lnGamma_mix - lnGammai)).sum();
            }
        }
        EigenArray get_lngamma(double T, const EigenArray &x) const override {
            auto startTime = std::chrono::high_resolution_clock::now();
            EigenArray lngamma(x.size());
            EigenArray psigma_mix = get_psigma_mix(x);
            EigenArray lnGamma_mix = get_Gamma(T, psigma_mix).log();
            auto midTime = std::chrono::high_resolution_clock::now();
            for (Eigen::Index i = 0; i < x.size(); ++i) {
                lngamma(i) = get_lngamma_resid(i, T, lnGamma_mix) + get_lngamma_comb(T, x, i);
            }
            auto endTime = std::chrono::high_resolution_clock::now();
            //std::cout << std::chrono::duration<double>(midTime - startTime).count() << " s elapsed (mix)\n"; 
            //std::cout << std::chrono::duration<double>(endTime - midTime).count() << " s elapsed (comps)\n";
            //std::cout << std::chrono::duration<double>(endTime - startTime).count() << " s elapsed (total)\n";
            return lngamma;
        }
    };

    struct COSMO3Constants {
        double 
            AEFFPRIME = 7.25, // [A^2]
            c_OH_OH = 4013.78, // [kcal A^4/(mol e^2)]
            c_OT_OT = 932.31, // [kcal A^4 /(mol e^2)]
            c_OH_OT = 3016.43, // [kcal A^4 /(mol e^2)]
            A_ES = 6525.69, // [kcal A^4 /(mol e^2)]
            B_ES = 1.4859e8, // [kcal A^4 K^2/(mol e^2)]
            N_A = 6.022140758e23, // [mol^{-1}]
            k_B = 1.38064903e-23, // [J K^{-1}]
            R = k_B*N_A/4184, // [kcal/(mol*K)]; Universal gas constant of new redefinition of 2018, https://doi.org/10.1088/1681-7575/aa99bc
            Gamma_rel_tol = 1e-8; // relative tolerance for Gamma in iterative loop
        bool fast_Gamma = false;
        std::string to_string() {
            return "NOT IMPLEMENTED YET";
            //return "c_hb: " + std::to_string(c_hb) + " kcal A^4 /(mol*e^2) \nsigma_hb: " + std::to_string(sigma_hb) + " e/A^2\nalpha_prime: " + std::to_string(alpha_prime) + " kcal A^4 /(mol*e^2)\nAEFFPRIME: " + std::to_string(AEFFPRIME) + " A\nR: " + std::to_string(R) + " kcal/(mol*K)";
        }
    };

    class COSMO3 : public AbstractCOSMOModel {
    public:
        typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
        typedef Eigen::Array<double, 51, 51> EigenMatrix5151;
        typedef Eigen::Array<double, Eigen::Dynamic, 1> EigenArray;
        typedef Eigen::Array<double, 51, 1> EigenArray51;
        typedef Eigen::Array<double, 51*3, 1> EigenArray153;
        
    private:
        COSMO3Constants m_consts;
        Eigen::Index ileft, w;
    public:
        COSMO3(const std::vector<std::string> &identifiers, const ProfileDatabase &datadb, const COSMO3Constants &constants = COSMO3Constants())
            : AbstractCOSMOModel(datadb, identifiers), m_consts(constants) {
            Eigen::Index iL, iR;
            std::tie(iL, iR) = get_nonzero_bounds();
            this->ileft = iL; this->w = iR - iL + 1;
        };

        /// Get access to the parameters that are in use
        COSMO3Constants &get_mutable_COSMO_constants() { return m_consts; }

        std::tuple<Eigen::Index, Eigen::Index> get_ileftw() const { return std::make_tuple(ileft, w);} ;

        double get_c_hb(profile_type type1, profile_type type2) const{

            if (type1 == type2){
               if (type1 == profile_type::OH_PROFILE) {
                    return m_consts.c_OH_OH;
                }
                else if (type1 == profile_type::OT_PROFILE) {
                    return m_consts.c_OT_OT;
                }
                else {
                    return 0.0;
                }
            }
            else if ((type1 == profile_type::OH_PROFILE && type2 == profile_type::OT_PROFILE) 
                  || (type1 == profile_type::OT_PROFILE && type2 == profile_type::OH_PROFILE)) {
                return m_consts.c_OH_OT;
            }
            else {
                return 0.0;
            }
        }
        EigenMatrix5151 get_DELTAW(double T, profile_type type_t, profile_type type_s) const {
            auto delta_sigma = 2*0.025/50;
            EigenMatrix5151 DELTAW;
            double cc = get_c_hb(type_t, type_s);
            for (auto m = 0; m < 51; ++m) {
                for (auto n = 0; n < 51; ++n) {
                    double sigma_m = -0.025 + delta_sigma*m,
                         sigma_n = -0.025 + delta_sigma*n,
                         c_hb = (sigma_m*sigma_n >= 0) ? 0 : cc,
                         c_ES = m_consts.A_ES + m_consts.B_ES/(T*T);
                    DELTAW(m, n) = c_ES*POW2(sigma_m + sigma_n) - c_hb*POW2(sigma_m-sigma_n);
                }
            }
            return DELTAW;
        }
        EigenMatrix5151 get_DELTAW_fast(double T, profile_type type_t, profile_type type_s) const {
            auto delta_sigma = 2 * 0.025 / 50;
            EigenMatrix5151 DELTAW; DELTAW.setZero();
            double cc = get_c_hb(type_t, type_s);
            for (auto m = ileft; m < ileft+w+1; ++m) {
                for (auto n = ileft; n < ileft+w+1; ++n) {
                    double sigma_m = -0.025 + delta_sigma*m,
                        sigma_n = -0.025 + delta_sigma*n,
                        c_hb = (sigma_m*sigma_n >= 0) ? 0 : cc,
                        c_ES = m_consts.A_ES + m_consts.B_ES / (T*T);
                    DELTAW(m, n) = c_ES*POW2(sigma_m + sigma_n) - c_hb*POW2(sigma_m - sigma_n);
                }
            }
            return DELTAW;
        }
        EigenMatrix get_AA(double T, EigenArray153 psigmas){
            // Build the massive \Delta W matrix that is 153*153 in size
            EigenMatrix DELTAW(153, 153);
            double R = m_consts.R;
            std::vector<profile_type> types = { profile_type::NHB_PROFILE, profile_type::OH_PROFILE, profile_type::OT_PROFILE };
            for (auto i = 0; i < 3; ++i) {
                for (auto j = 0; j < 3; ++j) {
                    DELTAW.matrix().block(51 * i, 51 * j, 51, 51) = get_DELTAW(T, types[i], types[j]);
                }
            }
            return Eigen::exp(-DELTAW/(R*T)).rowwise()*psigmas.transpose();
        }
        EigenArray153 get_Gamma(double T, EigenArray153 psigmas) const { // psigmas in NHB, OH, OT order
            auto startTime = std::chrono::high_resolution_clock::now();
            double R = m_consts.R;
            EigenArray153 Gamma, Gammanew; Gamma.setOnes(); Gammanew.setOnes();
            
            if (!m_consts.fast_Gamma){

                // Build the massive \Delta W matrix that is 153*153 in size
                EigenMatrix DELTAW(153, 153);
                std::vector<profile_type> types = { profile_type::NHB_PROFILE, profile_type::OH_PROFILE, profile_type::OT_PROFILE };
                for (auto i = 0; i < 3; ++i) {
                    for (auto j = 0; j < 3; ++j) {
                        DELTAW.matrix().block(51*i, 51*j, 51, 51) = get_DELTAW(T, types[i], types[j]);
                    }
                }

                EigenMatrix AA = Eigen::exp(-DELTAW / (R*T)).rowwise()*psigmas.transpose();
                for (auto counter = 0; counter < 251; ++counter) {
                    Gammanew = 1 / (AA.rowwise()*Gamma.transpose()).rowwise().sum();
                    Gamma = (Gamma + Gammanew) / 2;
                    double maxdiff = ((Gamma - Gammanew) / Gamma).cwiseAbs().real().maxCoeff();
                    if (maxdiff < m_consts.Gamma_rel_tol) {
                        break;
                    }
                }
                return Gamma;
            }
            else {
                // Build the massive AA matrix that is 153*153 in size
                auto midTime = std::chrono::high_resolution_clock::now();
                std::vector<profile_type> types = { profile_type::NHB_PROFILE, profile_type::OH_PROFILE, profile_type::OT_PROFILE };
                std::vector<Eigen::Index> offsets = {0*51, 1*51, 2*51};
                EigenMatrix AA(153,153);
                for (auto i = 0; i < 3; ++i) {
                    for (auto j = 0; j < 3; ++j) {
                        Eigen::Index rowoffset = offsets[i], coloffset = offsets[j];
                        AA.matrix().block(rowoffset + ileft, coloffset + ileft, w, w) = Eigen::exp(-get_DELTAW_fast(T, types[i], types[j]).block(ileft, ileft, w, w).array() / (R*T)).rowwise()*psigmas.segment(coloffset+ileft,w).transpose();
                    }
                }
                auto midTime2 = std::chrono::high_resolution_clock::now();

                for (auto counter = 0; counter < 251; ++counter) {
                    for (Eigen::Index offset : {51*0, 51*1, 51*2}){
                        Gammanew.segment(offset + ileft, w) = 1 / (
                              AA.matrix().block(offset+ileft,51*0+ileft,w,w).array().rowwise()*Gamma.segment(51*0+ileft, w).transpose()
                            + AA.matrix().block(offset+ileft,51*1+ileft,w,w).array().rowwise()*Gamma.segment(51*1+ileft, w).transpose()
                            + AA.matrix().block(offset+ileft,51*2+ileft,w,w).array().rowwise()*Gamma.segment(51*2+ileft, w).transpose()
                        ).rowwise().sum();
                    }
                    for (Eigen::Index offset : {51 * 0, 51 * 1, 51 * 2}) {
                        Gamma.segment(offset + ileft, w) = (Gamma.segment(offset + ileft, w) + Gammanew.segment(offset + ileft, w)) / 2;
                    }
                    double maxdiff = ((Gamma - Gammanew) / Gamma).cwiseAbs().real().maxCoeff();
                    if (maxdiff < m_consts.Gamma_rel_tol) {
                        break;
                    }
                }
                auto endTime = std::chrono::high_resolution_clock::now();
                //std::cout << std::chrono::duration<double>(midTime - startTime).count() << " s elapsed (DELTAW)\n"; 
                //std::cout << std::chrono::duration<double>(midTime2 - midTime).count() << " s elapsed (AA)\n";
                //std::cout << std::chrono::duration<double>(endTime - midTime2).count() << " s elapsed (comps)\n";
                //std::cout << std::chrono::duration<double>(endTime - startTime).count() << " s elapsed (total)\n";
                return Gamma;
            }
        }

        /**
        The residual part of ln(γ_i)
        */
        double get_lngamma_resid(std::size_t i, double T, const EigenArray &lnGamma_mix) const
        {
            double AEFFPRIME = m_consts.AEFFPRIME;
            EigenArray psigmas(3*51);
            const auto &f = m_fluids[i];
            double A_i = f.A_COSMO_A2;
            psigmas << f.profiles.nhb.psigma(A_i), f.profiles.oh.psigma(A_i), f.profiles.ot.psigma(A_i);
            double check_sum = psigmas.sum();
            EigenArray lnGammai = get_Gamma(T, psigmas).log();
            return A_i/AEFFPRIME*(psigmas*(lnGamma_mix - lnGammai)).sum();
        }
        EigenArray get_lngamma_disp(const EigenArray &x) const {
            if (x.size() != 2){ throw std::invalid_argument("Multi-component mixtures not supported for dispersive contribution yet"); }
            double w = 0.27027; // default value
            auto cls0 = m_fluids[0].dispersion_flag, cls1 = m_fluids[1].dispersion_flag;
            using d = FluidProfiles::dispersion_classes;
            if ((cls0 == d::DISP_WATER && cls1 == d::DISP_ONLY_ACCEPTOR)
                |
                (cls1 == d::DISP_WATER && cls0 == d::DISP_ONLY_ACCEPTOR)){ 
                w = -0.27027;
            }
            else if ((cls0 == d::DISP_WATER && cls1 == d::DISP_COOH)
                |
                (cls1 == d::DISP_WATER && cls0 == d::DISP_COOH)) {
                w = -0.27027;
            }
            else if ((cls0 == d::DISP_COOH && (cls1 == d::DISP_NHB || cls1 == d::DISP_DONOR_ACCEPTOR))
                |
                (cls1 == d::DISP_COOH && (cls0 == d::DISP_NHB || cls0 == d::DISP_DONOR_ACCEPTOR))) {
                w = -0.27027;
            }

            double ekB0 = m_fluids[0].dispersion_eoverkB, ekB1 = m_fluids[1].dispersion_eoverkB;
            double A = w*(0.5*(ekB0+ekB1) - sqrt(ekB0*ekB1));
            EigenArray lngamma_dsp(2);
            lngamma_dsp(0) = A*x[1]*x[1];
            lngamma_dsp(1) = A*x[0]*x[0];
            return lngamma_dsp;
        }
        EigenArray get_lngamma(double T, const EigenArray &x) const override {
            EigenArray lngamma(x.size());
            EigenArray psigmas(3*51);
            psigmas << get_psigma_mix(x, profile_type::NHB_PROFILE), get_psigma_mix(x, profile_type::OH_PROFILE), get_psigma_mix(x, profile_type::OT_PROFILE);
            double check_sum = psigmas.sum();
            EigenArray lnGamma_mix = get_Gamma(T, psigmas).log();
            EigenArray lngamma_disp = get_lngamma_disp(x);
            for (Eigen::Index i = 0; i < x.size(); ++i) {
                lngamma(i) = get_lngamma_resid(i, T, lnGamma_mix) + get_lngamma_comb(T, x, i) + lngamma_disp(i);
            }
            return lngamma;
        }
    };

}; /* namespace COSMOSAC */

#endif
