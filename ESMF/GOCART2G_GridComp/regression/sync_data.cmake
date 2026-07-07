# sync_data.cmake
# Downloads GOCART2G regression data from S3.
# Required inputs (via -D): LOCAL_DIR, WORK_DIR, ESMA_SYNC_DATA_SCRIPT

if(DEFINED ENV{LOCAL_REGRESSION_DATA_DIR} AND NOT "$ENV{LOCAL_REGRESSION_DATA_DIR}" STREQUAL "")
  message(STATUS "LOCAL_REGRESSION_DATA_DIR is set -- skipping S3 sync")
  return()
endif()

if(EXISTS "${LOCAL_DIR}")
  message(STATUS "Regression data already present: ${LOCAL_DIR} -- skipping sync")
  return()
endif()

set(AWS_CLI "aws")
set(S3_URI "s3://gmao-usw2-si-team/regression-data/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp")
set(API_KEY_ENV "AWS_API_KEY_GMAO_SITEAM_S3")
set(CREDENTIALS_URL "https://llvsm4u7ij.execute-api.us-east-1.amazonaws.com/credentials")
set(REQUIRED FALSE)

include("${ESMA_SYNC_DATA_SCRIPT}")
