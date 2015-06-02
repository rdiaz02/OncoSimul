#include <iostream>
//#include <exception>
#include <stdexcept>
#include <sstream>
// using namespace std;

class myEx: public std::runtime_error {
public:
  myEx(const std::string &s) :
    std::runtime_error(s) {}
};


void f1() {
  throw myEx("coco");
  std::cout << "\n hola\n" << std::endl;
}

int main() {
  int i = 0;
  try {
    f1();
    i = i + 1;
  } catch (myEx& e) {
    std::cout << e.what() << std::endl;
  }
  std:: cout << "\n i = " << i << std::endl;
      
}
  





