cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(fetch_git_repo_sub LANGUAGES NONE)

include(ExternalProject)

ExternalProject_Add(
  googletest
  SOURCE_DIR "/Users/JG/Documents/Software/OpenCAP/OpenCAP/opencap/build/_deps/googletest-src"
  BINARY_DIR "/Users/JG/Documents/Software/OpenCAP/OpenCAP/opencap/build/_deps/googletest-build"
  GIT_REPOSITORY
    https://github.com/google/googletest.git
  GIT_TAG
    release-1.8.1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  TEST_COMMAND ""
  )
