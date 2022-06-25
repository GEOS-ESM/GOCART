# UFS/GOCART interface
#
# UFS porting of ecBuild package.
# See original version at: https://github.com/JCSDA-internal/ecbuild (tag: 3.3.2.jcsda3)

macro (ecbuild_install_project)

  install (EXPORT ${PROJECT_NAME}-targets
           DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

endmacro ()
