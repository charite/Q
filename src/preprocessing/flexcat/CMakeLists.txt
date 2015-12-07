cmake_minimum_required (VERSION 2.8.2)
project (Flexcat)
message (STATUS "Configuring flexcat")

set (SEQAN_APP_VERSION "1.0.0")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package(CXX11 REQUIRED)

# Disable flexcat if no CXX11 available.
if (NOT CXX11_FOUND)
    message (STATUS "  C++11 compiler required to build flexcat")
    return ()
endif (NOT CXX11_FOUND)


find_package(CXX14 REQUIRED)

# Disable flexcat if no CXX11 available.
if (NOT CXX14_FOUND)
    message (STATUS "  C++14 compiler required to build flexcat")
    return ()
endif (NOT CXX14_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package (SeqAn REQUIRED)

find_package(Boost 1.36.0)
if(NOT Boost_FOUND)
    message (STATUS "  flexcat cannot be built without Boost.")
    return ()
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable(test_demultiplex        test_demultiplex.cpp demultiplex.h)
add_executable(test_trimming           test_trimming.cpp read_trimming.h)
add_executable(test_adapter            test_adapter.cpp adapter_trimming.h)
add_executable(test_general_processing test_general_processing.cpp general_processing.h)

foreach (TEST demultiplex trimming adapter general_processing)
    target_link_libraries (test_${TEST} ${SEQAN_LIBRARIES})
endforeach ()

add_library (flexlib
			 flexlib.cpp
			 flexlib.h
			 flexi_program.h
			 ptc.h
             read.h
			 read_writer.h
			 semaphore.h
             demultiplex.h
			 argument_parser.h
             read_trimming.h
             adapter_trimming.h
             general_processing.h
             helper_functions.h
			 general_stats.h
             )
target_link_libraries (flexlib ${SEQAN_LIBRARIES})

set (KNIME_FLEX_NODES flexcatQC flexcatAR flexcatDMulti flexcatFilter)

foreach (TARGET ${KNIME_FLEX_NODES} flexcat)
    add_executable (${TARGET} flexcat.cpp flexcat.h)
    target_link_libraries (${TARGET} flexlib ${SEQAN_LIBRARIES})
endforeach ()

set_target_properties (flexcatQC PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FlexiProgram::QUALITY_CONTROL")
set_target_properties (flexcatAR PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FlexiProgram::ADAPTER_REMOVAL")
set_target_properties (flexcatDMulti PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FlexiProgram::DEMULTIPLEXING")
set_target_properties (flexcatFilter PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FlexiProgram::FILTERING")

# Set CXX flags 
if (NOT MSVC)
set (CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS}")
endif (NOT MSVC)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${CXX14_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (flexcat)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install flexcat in ${PREFIX}/bin directory
install (TARGETS flexcat
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/flexcat for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (flexcat)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable flexcat in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} ${KNIME_FLEX_NODES} CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:flexcat")
  set (CPACK_PACKAGE_NAME "flexcat")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flexcat - Flexible barcode and adapter removal for NGS")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Benjamin Menkuec <benjamin@menkuec.de>")

  seqan_configure_cpack_app (flexcat "flexcat")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:flexcat")

