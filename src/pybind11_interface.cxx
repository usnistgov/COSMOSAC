// Includes from this library
#include "COSMO_SAC/COSMO.hpp"
using namespace COSMOSAC; 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

void init_COSMO(py::module &m) {

    py::class_<COSMO1Constants>(m, "COSMO1Constants")
        .def_readwrite("AEFFPRIME", &COSMO1Constants::AEFFPRIME)
        .def_readwrite("alpha_prime", &COSMO1Constants::alpha_prime)
        .def_readwrite("c_hb", &COSMO1Constants::c_hb)
        .def_readwrite("R", &COSMO1Constants::R)
        .def_readwrite("sigma_hb", &COSMO1Constants::sigma_hb)
        .def_readwrite("fast_Gamma", &COSMO1Constants::fast_Gamma)
        .def("__repr__", &COSMO1Constants::to_string);

    py::class_<COSMO3Constants>(m, "COSMO3Constants")
        .def_readwrite("AEFFPRIME", &COSMO3Constants::AEFFPRIME)
        .def_readwrite("A_ES", &COSMO3Constants::A_ES)
        .def_readwrite("B_ES", &COSMO3Constants::B_ES)
        .def_readwrite("c_OH_OH", &COSMO3Constants::c_OH_OH)
        .def_readwrite("c_OH_OT", &COSMO3Constants::c_OH_OT)
        .def_readwrite("c_OT_OT", &COSMO3Constants::c_OT_OT)
        .def_readwrite("R", &COSMO3Constants::R)
        .def_readwrite("Gamma_rel_tol", &COSMO3Constants::Gamma_rel_tol)
        .def_readwrite("fast_Gamma", &COSMO3Constants::fast_Gamma)
        .def("__repr__", &COSMO3Constants::to_string)
        ;

    py::class_<CombinatorialConstants>(m, "CombinatorialConstants")
        .def_readwrite("q0", &CombinatorialConstants::q0)
        .def_readwrite("r0", &CombinatorialConstants::r0)
        .def_readwrite("z_coordination", &CombinatorialConstants::z_coordination)
        .def("__repr__", &CombinatorialConstants::to_string);

    py::class_<SigmaProfile>(m, "SigmaProfile")
        .def_property_readonly("psigmaA", [](SigmaProfile &p) { return p.psigmaA(); })
        .def_property_readonly("psigma", [](SigmaProfile &p, double A_i) { return p.psigma(A_i); })
        .def_property_readonly("sigma", [](SigmaProfile &p) { return p.sigma(); });

    py::class_<SigmaProfileSet>(m, "SigmaProfileSet")
        .def_readonly("nhb", &SigmaProfileSet::nhb)
        .def_readonly("oh", &SigmaProfileSet::oh)
        .def_readonly("ot", &SigmaProfileSet::ot)
        ;

    py::class_<FluidProfiles>(m, "FluidProfiles")
        .def_readonly("profiles", &FluidProfiles::profiles)
        .def_readonly("name", &FluidProfiles::name)
        .def_readonly("VTnumber", &FluidProfiles::VTnumber)
        .def_readonly("A_COSMO_A2", &FluidProfiles::A_COSMO_A2)
        .def_readonly("V_COSMO_A3", &FluidProfiles::V_COSMO_A3);

    py::class_<ProfileDatabase >(m, "ProfileDatabase")
        .def("get_profile", &ProfileDatabase::get_profile)
        .def("normalize_identifier", &ProfileDatabase::normalize_identifier)
        ;

    py::class_<VirginiaTechProfileDatabase, ProfileDatabase >(m, "VirginiaTechProfileDatabase")
        .def(py::init<const std::string &, const std::string &>())
        .def("add_profile", &VirginiaTechProfileDatabase::add_profile)
        .def("to_JSON", &VirginiaTechProfileDatabase::to_JSON)
        ;

    py::class_<DelawareProfileDatabase, ProfileDatabase >(m, "DelawareProfileDatabase")
        .def(py::init<const std::string &, const std::string &>())
        .def("add_profile", &DelawareProfileDatabase::add_profile)
        .def("to_JSON", &DelawareProfileDatabase::to_JSON);

    py::class_<AbstractCOSMOModel >(m, "AbstractCOSMOModel")
        .def("get_lngamma_comb", &AbstractCOSMOModel::get_lngamma_comb)
        .def("get_sigma_profiles", &AbstractCOSMOModel::get_sigma_profiles)
        .def("get_mutable_combinatorial_constants", &AbstractCOSMOModel::get_mutable_combinatorial_constants, py::return_value_policy::reference)
        .def("get_psigma_mix", &AbstractCOSMOModel::get_psigma_mix)
        ;

    typedef COSMO1::EigenArray EigenArray1;
    py::class_<COSMO1, AbstractCOSMOModel>(m, "COSMO1")
        .def(py::init<const std::vector<std::string> &, const ProfileDatabase &>())
        .def("get_Gamma", &COSMO1::get_Gamma)
         .def("get_lngamma_resid", py::overload_cast<double, const EigenArray1&>(&COSMO1::get_lngamma_resid, py::const_))
         .def("get_lngamma_resid", py::overload_cast<std::size_t, double, const EigenArray1&>(&COSMO1::get_lngamma_resid, py::const_))
        .def("get_lngamma", &COSMO1::get_lngamma)
        .def("get_mutable_COSMO_constants", &COSMO1::get_mutable_COSMO_constants, py::return_value_policy::reference)
        ;

    typedef COSMO3::EigenArray EigenArray;
    py::class_<COSMO3, AbstractCOSMOModel>(m, "COSMO3")
        .def(py::init<const std::vector<std::string> &, const ProfileDatabase &>())
        .def("get_Gamma", &COSMO3::get_Gamma)
        .def("get_lngamma_resid", py::overload_cast<double, const EigenArray&>(&COSMO3::get_lngamma_resid, py::const_))
        .def("get_lngamma_resid", py::overload_cast<std::size_t, double, const EigenArray&>(&COSMO3::get_lngamma_resid, py::const_))
        .def("get_lngamma_disp", &COSMO3::get_lngamma_disp)
        .def("get_lngamma", &COSMO3::get_lngamma)
        .def("get_AA", &COSMO3::get_AA)
        .def("get_mutable_COSMO_constants", &COSMO3::get_mutable_COSMO_constants, py::return_value_policy::reference)
        ;
}

PYBIND11_MODULE(cCOSMO, m) {
    m.doc() = "COSMO-SAC implementation";
    init_COSMO(m);
}