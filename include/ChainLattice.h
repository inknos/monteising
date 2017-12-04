#ifndef CHAIN_LATTICE_H
#define CHAIN_LATTICE_H
 
template <typename T>
void bar(T t) {}

void foo2() {}

template <typename Car, typename... Cdr>
void foo2(Car car, Cdr... cdr)
{
  bar(car);
  foo2(cdr...);
}
#endif
