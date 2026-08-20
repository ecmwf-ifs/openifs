// (C) Copyright 2003- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <iostream>
#include <cstdlib>
#include <exception>

namespace my_namespace {

class MyException : public std::exception {
public:
    std::string msg_;
    MyException(const std::string& msg) : std::exception(), msg_(msg) {}
    virtual ~MyException() noexcept {}
    virtual const char* what() const noexcept { return msg_.c_str(); }
};

void function_2 () {
    throw MyException("problem in function_2");
}

void function_1 () {
    function_2();
}

} // namespace my_namespace

extern "C" {
    void my_namespace__function_1() { my_namespace::function_1(); }
}
