context("C++")
test_that("Catch unit tests pass", {
  skip_on_os("solaris")
  expect_cpp_tests_pass("mdgc")
})
