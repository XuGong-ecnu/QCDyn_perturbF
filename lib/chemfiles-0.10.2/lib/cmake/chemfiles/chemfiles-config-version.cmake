set(PACKAGE_VERSION "0.10.2")

if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)
  set(PACKAGE_VERSION_EXACT TRUE)
endif()

# Assume true until shown otherwise
set(PACKAGE_VERSION_COMPATIBLE TRUE)

if(PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
endif()

if(0 EQUAL 0)
    if(NOT PACKAGE_FIND_VERSION_MINOR EQUAL "10")
        set(PACKAGE_VERSION_COMPATIBLE FALSE)
    endif()
else()
    if(NOT PACKAGE_FIND_VERSION_MAJOR EQUAL "0")
        set(PACKAGE_VERSION_COMPATIBLE FALSE)
    endif()
endif()


# if the installed or the using project don't have CMAKE_SIZEOF_VOID_P set, ignore it:
if(NOT "${CMAKE_SIZEOF_VOID_P}" STREQUAL "" AND NOT "8" STREQUAL "")
  # check that the installed version has the same 32/64bit-ness
  # as the one which is currently searching:
  if(NOT CMAKE_SIZEOF_VOID_P STREQUAL "8")
    math(EXPR installedBits "8 * 8")
    set(PACKAGE_VERSION "${PACKAGE_VERSION} (${installedBits}bit)")
    set(PACKAGE_VERSION_UNSUITABLE TRUE)
  endif()
endif()
