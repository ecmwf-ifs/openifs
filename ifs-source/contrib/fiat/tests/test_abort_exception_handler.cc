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
