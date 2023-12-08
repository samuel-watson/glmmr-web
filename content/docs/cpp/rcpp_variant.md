---
title: Rcpp design with multiple template classes
weight: 1
---

# Rcpp with multiple template classes
One of the difficulties with building this package was trying to keep the C++ code relatively tidy and readable. The main functionality is contained in the glmmr headers, but an Rcpp interface is required to translate to R. However, there are several template classes for the different types of model. While these classes all have the same named functions, one has to be able to choose between the appropriate classes at run time in R. There are also multiple different return types for different functions. This can lead to very verbose code, requiring large switch or conditional statements duplicating the function for each type of class. As an alternative we can use `std::variant` and `std::visit` to make the code more compact and readable. As this took me a while to work out, I thought it would be useful to explain here in case anyone else finds a similar issue.

A useful feature of Rcpp is `Rcpp::XPtr`, a smart pointer, which means we can instantiate a class object in R and store the pointer to pass to other functions. For example, if we have a class `A` then we could create a new pointer to it with a function:
```
using namespace Rcpp;
SEXP new_A(){
    XPtr<A> ptr(new A());
    return A;
}
```
which we can then call in R. The pointer can then be passed to another function to do something else:
```
SEXP fn_A(SEXP ptr){
    Xptr<A> aptr(ptr);
    double result = aptr->fn();
    return wrap(result);
}
```
which will return the double `result` to R. However, let's say we instead have a templated class `B` with multiple types, but all with a function `fn()`. One option in Rcpp may be to take an argument that identifies which class type we have and then using a switch or conditional statement:
```
SEXP fn_A(SEXP ptr, int type){
    switch(type){
        case 1:
        {
            Xptr<A<type1> > aptr(ptr);
            double result = aptr->fn();
            return wrap(result);
            break;
        }
        case 2:
        {
            Xptr<A<type2> > aptr(ptr);
            double result = aptr->fn();
            return wrap(result);
            break;
        }
        ...
    }
    
}
```
but this can become very long and unnecessarily repetitive, especially if we have a large number of functions. 

## Using std::variant
As an alternative, we can use `std::variant` and `std::visit` to make the code for each function much more compact and prevent too much repetition. My approach has been to use a `struct` that can take and store an `Xptr` to any of the relevant classes. In a header file (such as `myproject.h`) we can include:
```
#pragma once

enum class Type {
    type1 = 0,
    type2 = 1, 
    type3 = 2
}

struct objectType
{
    std::variant<int, XPtr<B<Type1> >, XPtr<B<Type2> >, XPtr<B<Type3> > > ptr;
    objectType(SEXP xptr, Type type) : ptr(0) {
        switch(Type){
            case type1:
                ptr = XPtr<B<Type1> >(xptr);
                break;
            case type2:
                ptr = XPtr<B<Type2> >(xptr);
                break;
            case type3:
                ptr = XPtr<B<Type3> >(xptr);
                break;
        }
    }
}

using returnType = std::variant<int, double, std::string>;

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

```
We've included an int in the variant so that the variant can be initialised properly (an `XPtr` can't be initialised with `nullptr` so we can't default initialise the variant as far as I can tell) and assuming that `Type1`, `Type2`, and `Type3` name a class type. We have also defined a variant for the return type and an overloaded class for the functor.

Now, we can use the following template for Rcpp functions in our .cpp file for functions with and without a return value: 
```
#include "myproject.h"

// [[Rcpp::export]]
SEXP fn(SEXP xp, int type = 0){
  objectType obj(xp,static_cast<Type>(type));
  auto functor = overloaded {
    [](int) {  return returnType(0);}, 
    [](auto ptr){return returnType(ptr->fn());}
  };
  auto S = std::visit(functor,obj.ptr);
  return wrap(std::get<int>(S));
}

// [[Rcpp::export]]
void fn2(SEXP xp, int type = 0){
  objectType obj(xp,static_cast<Type>(type));
  auto functor = overloaded {
    [](int) {}, 
    [](auto ptr){ptr->fn2();}
  };
  auto S = std::visit(functor,obj.ptr);
}
```
We still need to provide an argument specifying the appropriate class, however, this can be wrapped in an appropriate class (like S4) on the R side that selects this properly to prevent errors and crashes. This structure can then be used for all interfaces with the templated classes and with all return types, we only need to change the function and return type. We can pass arguments to the function the standard way for lambda, for example, to provide access to objects in the function's scope by reference we can use `[&](auto ptr){return returnType(ptr->fn(arg));}`. I'm sure there are further refinements one could make to this scheme, but I found it very useful to massively reduce the size of the C++ files without repeating large volumes of code.

This is much more compact, but still repetitive as the only thing changing between function calls is the return type and function name. We could make this a little more compact, for example, by defining in our header:
```
template<typename T>
struct Fn
{
  objectType obj;
  Fn(SEXP xp, int type = 0) : model(xp,static_cast<Type>(type)) {};
  template<class Visitor>
  constexpr T operator()(Visitor vis){
    auto S = std::visit(vis,obj.ptr);
    return std::get<T>(S);
  }
};
```
and then our function would look like this
```
// [[Rcpp::export]]
SEXP Model__P(SEXP xp, int type = 0){
  auto functor = overloaded {
    [](int) {  return returnType(0);},
    [](auto ptr){return returnType(ptr->model.linear_predictor.P());}
  };
  Fn<int> func(xp,type);
  return wrap(func(functor));
}
```
However, this only removes one line of code per function and doesn't deal with the case of passing arguments. Function pointers also don't help much as we would have to get the pointer for each function template. Will update if I find a better solution!
 