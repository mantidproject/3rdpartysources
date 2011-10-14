#Description: Native build for google mock
#Author: Owen Arnold @ Tessella/ISIS

cmake_minimum_required(VERSION 2.6)
INCLUDE_DIRECTORIES ( include gtest gtest/include )

SET(PROJECT_SOURCES src/gmock-cardinalities.cc src/gmock-internal-utils.cc src/gmock-matchers.cc src/gmock-printers.cc src/gmock-spec-builders.cc src/gmock.cc gtest/src/gtest-all.cc)

ADD_LIBRARY(gmock STATIC ${PROJECT_SOURCES})