#------------------------------------------------------------------------------
# DGtal Configuration file for the install target
#------------------------------------------------------------------------------

#--- Headers installation
install(DIRECTORY "${PROJECT_SOURCE_DIR}/src/" 
  DESTINATION "${INSTALL_INCLUDE_DIR}/"
  FILES_MATCHING PATTERN "*.*h")

#-- specific install for Config.h.in
install(DIRECTORY "${PROJECT_BINARY_DIR}/src/DGtal/" 
  DESTINATION "${INSTALL_INCLUDE_DIR}/DGtal/" 
  FILES_MATCHING PATTERN "*.*h")
INCLUDE(${CMAKE_MODULE_PATH}/DGtalConfigGenerator.cmake)



if(WITH_VISU3D_OGRE)
  #--- Ogre viewer ressources installation
  install(DIRECTORY "${PROJECT_SOURCE_DIR}/src/DGtal/io/viewers/OGRE/Ressources" 
    DESTINATION " ${DGTAL_CMAKE_DIR}/OGRE/"
    FILES_MATCHING PATTERN "*.*")
  
  install(FILES
    "${PROJECT_SOURCE_DIR}/cmake/FindOgreDGtal.cmake"
    DESTINATION "${DGTAL_CMAKE_DIR}" COMPONENT dev)
  
  install (CODE "set(INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR})
  configure_file(  ${PROJECT_SOURCE_DIR}/src/DGtal/io/viewers/OGRE/RequirementsInstall.h.in
    ${INSTALL_INCLUDE_DIR}/DGtal/io/viewers/OGRE/Requirements.h )")
  install(CODE  "configure_file(  ${PROJECT_SOURCE_DIR}/src/DGtal/io/viewers/OGRE/Ressources/plugins.cfg.in
    ${DGTAL_CMAKE_DIR}/OGRE/Ressources/plugins.cfg )")
  IF (APPLE)
    install(CODE  "configure_file(  ${PROJECT_SOURCE_DIR}/src/DGtal/io/viewers/OGRE/Ressources/info.plist.in
    ${DGTAL_CMAKE_DIR}/OGRE/Ressources/info.plist.in )")
  ENDIF(APPLE)
endif(WITH_VISU3D_OGRE)


