

if (  WITH_VISU3D_QGLVIEWER )
  SET(DGTAL_EXAMPLES_QGL_SRC
    digitalSetFromVol
    digitalSetFromPointList
    )
  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC}) 
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} DGtalIO DGtal
      ${DGtalLibDependencies})   
  ENDFOREACH(FILE)
  
add_subdirectory(viewers)
endif ( WITH_VISU3D_QGLVIEWER   )

add_subdirectory(boards)

