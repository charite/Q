# ===========================================================================
#                  Nexus-Tools
# ===========================================================================

# ----------------------------------------------------------------------------
# Load SeqAn build system
# ----------------------------------------------------------------------------

project (NexusTools)

if(LINUX_STATIC)
    cmake_minimum_required (VERSION 3.1)
    message("building with linux static libs")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " -static -static-libgcc -static-libstdc++")
    set(Boost_USE_STATIC_LIBS ON)
    set(SEQAN_LIBRARIES ${SEQAN_LIBRARIES} "-Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
else()
    cmake_minimum_required (VERSION 2.8.7)
endif(LINUX_STATIC)

set (CMAKE_INCLUDE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/include
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/util/cmake
    ${CMAKE_MODULE_PATH})
include (SeqAnContribs)
set (SEQAN_USE_SEQAN_BUILD_SYSTEM TRUE CACHE INTERNAL "Use SeqAn build system." FORCE)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}../seqan" CACHE INTERNAL "Root source directory." FORCE)
option(LINUX_STATIC "build with static libs" OFF)

# include(include/seqan/CMakeLists.txt)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)
include (SeqAnCtdSetup)

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

add_subdirectory(nexcat)
add_subdirectory(flexcat)

if(LINUX_STATIC)
	SET_TARGET_PROPERTIES(nexcat PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexlib PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexcat PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexcatAR PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexcatDMulti PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexcatFilter PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(flexcatQC PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_demultiplex PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_general_processing PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_trimming PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_adapter PROPERTIES LINK_SEARCH_END_STATIC 1)	
endif(LINUX_STATIC)
