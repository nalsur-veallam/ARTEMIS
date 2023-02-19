#ifndef MY_ERROR_CPP
#define MY_ERROR_CPP

#include <stdexcept>
#include <string>

class My_Error : public std::runtime_error {
public:
  My_Error(const std::string &msg = "") : std::runtime_error(msg) {}
};

#endif
