#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "util/data.h"
#include "util/Wrapper.h"
#include "util/Residue_Representation.h"
#include "util/Entropy_Matrix.h"
#include "util/types.h"
#include "util/My_Error.cpp"

namespace py = pybind11;

class PyWrapper : public Wrapper {
public:
    /* Inherit the constructors */
    using Wrapper::Wrapper;

    /* Trampoline (need one for each virtual function) */

    std::vector <int> getNumbers() override {
        PYBIND11_OVERLOAD_PURE(
            std::vector <int>, /* Return type */
            Wrapper,      /* Parent class */
            getNumbers,          /* Name of function in C++ (must match Python name) */
                 /* Argument(s) */
        );
    }

};

int add(int i, int j, bool lin, float dt0, float dt1, float dt2, char const * file1, char const * file2) {

    double C1 = 1;
    double C2 = 1;
    double C = 1;
    double C0 = 1;
    double Cc = 0;

    if (!lin) {
        unsigned int it0 = 0;
        unsigned int it1 = 0;
        unsigned int it2 = 0;

        if (std::find(data::dts.begin(), data::dts.end(), dt0) != data::dts.end()) {

            it0 = std::find(data::dts.begin(), data::dts.end(), dt0) - data::dts.begin();

        } else {

            double diff = std::abs(data::dts[it0] - dt0);

            for (size_t it = 0; it < data::dts.size(); it++) {
                if (diff > std::abs(data::dts[it] - dt0) ) {

                    it0 = (unsigned int) it;
                    diff = std::abs(data::dts[it] - dt0);

                }
            }

            if (dt0 < 1) { C0 = data::approx(dt0, it0);} // TODO: parabolic approximation
            else { C0 = data::approx(dt0, it0);}

        }

        if (std::find(data::dts.begin(), data::dts.end(), dt1) != data::dts.end()) {

            it1 = std::find(data::dts.begin(), data::dts.end(), dt1) - data::dts.begin();

        } else {

            double diff = std::abs(data::dts[it1] - dt1);

            for (size_t it = 0; it < data::dts.size(); it++) {
                if (diff > std::abs(data::dts[it] - dt1) ) {

                    it1 = (unsigned int) it;
                    diff = std::abs(data::dts[it] - dt1);

                }
            }

            if (dt1 < 1) { C1 = data::approx(dt1, it1);} // TODO: parabolic approximation
            else { C1 = data::approx(dt1, it1);}

        }

        if (std::find(data::dts.begin(), data::dts.end(), dt2) != data::dts.end()) {

            it2 = std::find(data::dts.begin(), data::dts.end(), dt2) - data::dts.begin();

        } else {

            double diff = std::abs(data::dts[it2] - dt2);

            for (size_t it = 0; it < data::dts.size(); it++) {
                if (diff > std::abs(data::dts[it] - dt2) ) {

                    it2 = (unsigned int) it;
                    diff = std::abs(data::dts[it] - dt2);

                }
            }

            if (dt2 < 1) { C2 = data::approx(dt2, it2);} // TODO: parabolic approximation
            else { C2 = data::approx(dt2, it2);}
        }

        C = data::noise[it2]/data::noise[it1];
        Cc = data::noise[it0]/data::noise[it1];
    }

    double Cm1, Cm2;

    if (!lin) { Cm1 = ((Cc*C0 - C2*C)/(C1 - C2*C)); Cm2 = ((C1-C0*Cc)/(C1 - C2*C));}
    else { Cm1 = (dt2-dt0)/(dt2-dt1); Cm2 = -(dt1-dt0)/(dt2-dt1);}

    Interface wrapper(file1, file2, Cm1, Cm2);
    std::vector <int> numbers = wrapper.getNumbers();

    return i + j;
}

PYBIND11_MODULE(_test, m) {
    m.def("add", &add, "A function that adds two numbers");

    py::class_<Wrapper, PyWrapper /* <--- trampoline*/>(m, "Wrapper")
        .def(py::init<char const *, char const *, double, double>())
        .def("getNumbers", &Wrapper::getNumbers);

    py::class_<Interface, Wrapper>(m, "Interface")
        .def(py::init<char const *, char const *, double, double>());
}

