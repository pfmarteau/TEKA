file(REMOVE_RECURSE
  "libkats.pdb"
  "libkats.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang)
  include(CMakeFiles/kats.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
