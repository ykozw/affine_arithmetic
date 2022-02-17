cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
include(FetchContent)

# ソリューションエクスプローラーにcmake用のフォルダを追加する
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake")

# 本体プロジェクト
project(AA VERSION 1.0)
set(CMAKE_CXX_STANDARD 23)

# ソースファイル
set(SRC_FILES
    main.cpp
    plot.py
    CMakeLists.txt
    )

# pybind11をダウンロードして配置
FetchContent_Populate(
  pybind11
  URL      https://github.com/pybind/pybind11/archive/refs/tags/v2.9.0.zip
)
include_directories(${pybind11_SOURCE_DIR}/include)
add_subdirectory(${pybind11_SOURCE_DIR}/)

# プロジェクト作成
pybind11_add_module(AA ${SRC_FILES})
set_target_properties(AA PROPERTIES SUFFIX ".pyd")

# インストール
INSTALL(
  TARGETS AA LIBRARY 
  DESTINATION ${CMAKE_SOURCE_DIR}/bin
)

# HACK: ビルドと同時にinstallする
add_custom_command(
  TARGET AA 
    POST_BUILD COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --config $<CONFIG>
)