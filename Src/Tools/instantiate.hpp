/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/
#ifndef INSTANTIATE_H_
#define INSTANTIATE_H_

#define INSTANTIATE_DIMS(macro_target, type) \
  macro_target(type, 1u)                     \
  macro_target(type, 2u)                     \
  macro_target(type, 3u)

#ifdef instantiate
  #define instantiateTYPE(type) INSTANTIATE_DIMS(instantiate, type)
  instantiateTYPE(float)
  instantiateTYPE(double)
  instantiateTYPE(long double)
  instantiateTYPE(std::complex<float>)
  instantiateTYPE(std::complex<double>)
  instantiateTYPE(std::complex<long double>)
  #undef instantiateTYPE
  #undef instantiate
#endif

#ifdef instantiate_real
  #define instantiateREAL(type) INSTANTIATE_DIMS(instantiate_real, type)
  instantiateREAL(float)
  instantiateREAL(double)
  instantiateREAL(long double)
  #undef instantiateREAL
  #undef instantiate_real
#endif

#undef INSTANTIATE_DIMS
#endif
