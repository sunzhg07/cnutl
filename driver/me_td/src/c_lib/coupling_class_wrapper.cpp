//test_class_wrapper.cpp
#include <boost/python.hpp>
#include<boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/to_python_converter.hpp>
#include "coupling_coeff.h"

void ceshi(){
    std::cout << "ceshi" << std::endl;
}

BOOST_PYTHON_MODULE(Coupling_Module)
{
    using namespace boost::python;
    // 导出类
    class_<Coupling>("Coupling")                            //如果默认构造函数没有参数，可以省略
        //.def(init<int>())                               //其他构造函数
        .def("add", &Coupling::add)                            //成员函数
        .def("frac", &Coupling::factorial)                            //成员函数
        .def("CG", &Coupling::CGcoeff)                            //成员函数
        .def("SixJ", &Coupling::SixJ)                            //成员函数
        .def("NinJ", &Coupling::NinJ)                            //成员函数
        //.def_readwrite("publicVal", &A::publicVal)      //数据成员，当然是公共的
    ;
    //def("printA", &printA);
    //def("addA", &addA);
    def("ceshi",&ceshi);
}
