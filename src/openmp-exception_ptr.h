#ifndef OPENMP_EXCEPTION_PTR_H
#define OPENMP_EXCEPTION_PTR_H
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdexcept>
// #include <exception>

// see https://stackoverflow.com/q/11828539/5861244
class openmp_exception_ptr {
  // std::exception_ptr Ptr = nullptr;
  bool is_set = false;
public:
  void rethrow_if_error(){
#ifdef _OPENMP
    // if(this->Ptr)
    //   std::rethrow_exception(this->Ptr);
    if(is_set)
      throw std::runtime_error("Some exception occured. Further details cannot be provided because of https://stackoverflow.com/q/66362932/5861244.");
#endif
  }

  template <typename Function, typename... Parameters>
  void run(Function f, Parameters... params)
  {
#ifdef _OPENMP
    try
    {
      f(params...);
    }
    catch (...)
    {
#pragma omp atomic write
      is_set = true;
// #pragma omp critical(openmp_exception_ptr)
//       {
//         if(!is_set){
//           this->Ptr = std::current_exception();
//           is_set = true;
//         }
//       }
    }
#else
   f(params...);
#endif
  }
};

#endif // OPENMP_EXCEPTION_PTR_H
