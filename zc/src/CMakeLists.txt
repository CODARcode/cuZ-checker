set (zc_sources
  DynamicByteArray.c       ZC_CompareData.c           ZC_DataProperty_float.c  ZC_ReportGenerator.c     ZC_rw.c
  DynamicDoubleArray.c     ZC_CompareData_double.cpp  ZC_DataSetHandler.c      ZC_conf.c                ZC_util.c
  DynamicFloatArray.c      ZC_CompareData_float.cpp   ZC_gnuplot.c             dictionary.c
  DynamicIntArray.c        ZC_DataProperty.c          ZC_Hashtable.c           ZC_latex.c               iniparser.c
  ZC_ByteToolkit.c         ZC_DataProperty_double.c   ZC_quicksort.c           zc.c
  zc_derivativeMetrics.cpp ZC_ssim.c                  cuZC_entry.cu            cuZC_ssim.cu             timingGPU.cu
  cuZC_typeOne.cu          cuZC_typeTwo.cu            cuZC_typeThree.cu        cuZC_derivatives.cu
)
set_source_files_properties(ZC_CompareData_double.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(ZC_CompareData_float.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(zc_derivativeMetrics.cpp PROPERTIES LANGUAGE CUDA)

# TBA: ZC_R_math.c // R

if (FFTW_FOUND)
  list (APPEND zc_sources ZC_FFTW3_math.c)
  list (APPEND zc_dependencies ${FFTW_LIBRARIES})
endif ()

if (R_FOUND)
  list (APPEND zc_dependencies ${R_LIBRARY_BASE} ${R_LIBRARY_BLAS} ${R_LIBRARY_LAPACK} ${R_LIBRARY_READLINE})
endif ()

if (Boost_FOUND)
  list (APPEND zc_dependencies zserver)
endif ()

add_library (zc SHARED ${zc_sources})
#target_compile_features (zc PUBLIC cxx_std_11)
set_property(TARGET zc PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (zc
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
#target_link_libraries (zc ${zc_dependencies} m)

install (TARGETS zc LIBRARY DESTINATION lib)
