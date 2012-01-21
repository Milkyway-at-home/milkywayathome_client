

macro(add_opa_test name)
  add_executable(${name} ${name}.c opa_test.h)
  target_link_libraries(${name} opa)
  if(OPA_HAVE_LIBPTHREAD)
    target_link_libraries(${name} pthread)
  endif()

  add_test(NAME ${name} COMMAND ${name})
endmacro()


add_opa_test(sanity)
add_opa_test(test_primitives)
add_opa_test(test_barriers)
add_opa_test(test_queue)



