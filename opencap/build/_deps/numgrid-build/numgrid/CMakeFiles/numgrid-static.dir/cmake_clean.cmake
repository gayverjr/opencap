file(REMOVE_RECURSE
  "../../../lib/libnumgrid.a"
  "../../../lib/libnumgrid.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/numgrid-static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
