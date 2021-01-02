# Install script for directory: /homes/xyu/cuZ-checker/cuZ-checker/zc/include

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/DynamicByteArray.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_ByteToolkit.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_Hashtable.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_latex.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/dictionary.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/DynamicDoubleArray.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_CompareData.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_ReportGenerator.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_quicksort.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/iniparser.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/DynamicFloatArray.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_DataProperty.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_conf.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_rw.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/zc.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/DynamicIntArray.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_DataSetHandler.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_gnuplot.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_util.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/der.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/matrix.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/m_basic.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/m_der.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/m_der_metrics.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/m_norm.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/norm.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/others.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ssim.hpp"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/ZC_ssim.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/cuZC_entry.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/cuZC_ssim.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/timingGPU.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/cuZC_typeOne.h"
    "/homes/xyu/cuZ-checker/cuZ-checker/zc/include/cuZC_typeTwo.h"
    )
endif()

