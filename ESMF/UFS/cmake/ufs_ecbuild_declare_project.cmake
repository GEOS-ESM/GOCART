# UFS/GOCART interface
#
# UFS porting of ecBuild package.
# See original version at: https://github.com/JCSDA-internal/ecbuild (tag: 3.3.2.jcsda3)

macro (ecbuild_declare_project)

  include (GNUInstallDirs)
  set  (PROJECT_TARGETS_FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake")
  file (REMOVE ${PROJECT_TARGETS_FILE})

endmacro ()
