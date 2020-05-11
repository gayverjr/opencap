file(REMOVE_RECURSE
  "../../../lib/libnumgrid.dylib"
  "../../../lib/libnumgrid.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/numgrid-shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
