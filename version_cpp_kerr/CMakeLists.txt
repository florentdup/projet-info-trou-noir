
PYTHON_ADD_MODULE(kerr main.cpp)
FILE(COPY kerr.py DESTINATION .)
ADD_TEST(NAME version_cpp_kerr COMMAND ${PYTHON_EXECUTABLE} kerr.py)
