include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/main/include)
include_directories(${CMAKE_SOURCE_DIR}/stream/gear/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/UT/common/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/UT/PrVeloUT/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/UT/UTDecoding/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/SciFi/common/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/velo/common/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/event_model/common/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/event_model/velo/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/event_model/UT/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/event_model/SciFi/include)
include_directories(${CMAKE_SOURCE_DIR}/cuda/SciFi/PrForward/include)

file(GLOB x86MomentumForward "src/*cpp" )

cuda_add_library(x86MomentumForward STATIC
  ${x86MomentumForward}
)

if ( ROOT_FOUND )
  target_compile_definitions(x86MomentumForward PUBLIC WITH_ROOT)
  target_link_libraries(x86MomentumForward
    ${ROOT_LIBRARIES}
  )
  target_include_directories(x86MomentumForward PUBLIC ${ROOT_INCLUDE_DIRS})
endif()
