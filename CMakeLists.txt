cmake_minimum_required (VERSION 3.20)
project (upc-reader)

#sources for LIB_READER
set (SRCS
  src/Plot.cxx
  src/PlotAnalysis.cxx
  src/Fit.cxx
  src/CrossSectionMaker.cxx
  )

#name of LIB_READER
set (LIB_READER upc-reader)

#binary executables
set (BIN PlotsManager)

#UPC picoDsts library
set (UPC_INC include)

# --------------------------------------------------------
# HERE: possibly declare additional paths to other libraries

#set (OTHER_LIB_INC /path/to/other/library/include)

# --------------------------------------------------------

#general compiler options
set (CMAKE_CXX_COMPILER /usr/bin/g++)   #HERE change according to your compiler
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)


#option -Wl,--no-as-needed is necessary only with root5
set (CMAKE_CXX_FLAGS "-Wall -g -Wl,--no-as-needed")

include_directories (include ${UPC_INC})
# --------------------------------------------------------
# HERE: if declared additional libraries, add them here and link to them in the end of this file

#include_directories (${OTHER_LIB_INC})

# --------------------------------------------------------


set (ROOT630 /star/u/mvranovsk/root630/bin/root-config) # HERE change according to your ROOT version installation

#compiler flags
execute_process(COMMAND ${ROOT630} --cflags OUTPUT_VARIABLE ROOT_FLAGS_CMD)
string(REPLACE "\n" "" ROOT_FLAGS "${ROOT_FLAGS_CMD}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_FLAGS}")

#libraries
set(ROOT_DEPS Core EG Hist Physics RIO Tree MathCore RooFit RooFitCore)


#include dir
execute_process(COMMAND ${ROOT630} --libdir OUTPUT_VARIABLE ROOT_LIBDIR_CMD)
string(REPLACE "\n" "" ROOT_LIBDIR "${ROOT_LIBDIR_CMD}")
link_directories (${ROOT_LIBDIR})


#end of ROOT section

#create library LIB_READER
add_library (${LIB_READER} SHARED ${SRCS})
target_link_libraries(${LIB_READER} ${ROOT_DEPS} )

#create executables in BIN
foreach(IBIN ${BIN})
  add_executable(${IBIN} src/${IBIN}.cxx)
  target_link_libraries(${IBIN} ${LIB_READER} ${ROOT_DEPS} )
endforeach(IBIN ${BIN})
