file(GLOB UT_common "common/src/*cu")
file(GLOB UT_decoding "UTDecoding/src/*cu")
file(GLOB UT_consolidate "consolidate/src/*cu")
file(GLOB CompassUT_tracking "compassUT/src/*cu")

include_directories(${CMAKE_SOURCE_DIR}/main/include)
include_directories(${CMAKE_SOURCE_DIR}/stream/gear/include)
include_directories(${CMAKE_SOURCE_DIR}/stream/setup/include)
include_directories(../velo/common/include)
include_directories(../event_model/velo/include)
include_directories(../event_model/UT/include)
include_directories(../event_model/common/include)
include_directories(../utils/sorting/include)
include_directories(../utils/binary_search/include)
include_directories(common/include)
include_directories(UTDecoding/include)
include_directories(consolidate/include)
include_directories(compassUT/include)

if(TARGET_DEVICE STREQUAL "CPU")
  foreach(source_file ${UT_common} ${UT_decoding} ${UT_consolidate} ${CompassUT_tracking})
    set_source_files_properties(${source_file} PROPERTIES LANGUAGE CXX)
  endforeach(source_file)
endif()

allen_add_library(UTCommon STATIC
  ${UT_common}
)

allen_add_library(UT STATIC
  ${UT_decoding}
  ${UT_consolidate}
  ${CompassUT_tracking}
)

target_link_libraries(UT UTCommon)
