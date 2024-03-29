add_executable (analyzeDataProperty analyzeDataProperty.c)
target_link_libraries (analyzeDataProperty zc)

add_executable (analyzeDataProperty_multivars analyzeDataProperty_multivars.c)
target_link_libraries (analyzeDataProperty_multivars zc)

add_executable (compareDataSets compareDataSets.cu)
set_property(TARGET compareDataSets PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (compareDataSets
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries (compareDataSets zc)

add_executable (test_cuderivative test_cuderivative.cu)
target_link_libraries (test_cuderivative zc)

add_executable (generateGNUPlot generateGNUPlot.c)
target_link_libraries (generateGNUPlot zc)

add_executable (generateReport generateReport.c)
target_link_libraries (generateReport zc)

add_executable (modifyZCConfig modifyZCConfig.c)
target_link_libraries (modifyZCConfig zc)

#add_executable (runOfflineCase runOfflineCase.c)
#target_link_libraries (runOfflineCase zc)

add_executable (testRscript_double testRscript_double.c)
target_link_libraries (testRscript_double zc)

add_executable (testRscript_float testRscript_float.c)
target_link_libraries (testRscript_float zc)

add_executable (testRscript_readBinFloat testRscript_readBinFloat.c)
target_link_libraries (testRscript_readBinFloat zc)

install (TARGETS analyzeDataProperty analyzeDataProperty_multivars compareDataSets
  generateGNUPlot generateReport modifyZCConfig
         RUNTIME DESTINATION bin)

if (MPI_FOUND)
  add_executable (analyzeDataProperty_online analyzeDataProperty_online.c)
  target_link_libraries (analyzeDataProperty_online zc)

  add_executable (compareDataSets_online compareDataSets_online.c)
  target_link_libraries (compareDataSets_online zc)

  if (SZ_FOUND)
    add_executable (heatdis heatdis.c)
    target_link_libraries (heatdis zc ${SZ_LIBRARIES} ${ZLIB_LIBRARIES})
  endif ()
endif ()

#if (SZ_FOUND)
#  add_executable (testdouble_CompDecomp testdouble_CompDecomp.c)
#  target_link_libraries (testdouble_CompDecomp zc)
#endif ()
