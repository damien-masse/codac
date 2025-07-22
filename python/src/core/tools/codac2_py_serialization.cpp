/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <streambuf>
#include <istream>
#include <codac2_serialization.h>
#include "codac2_py_serialization_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

// py_istreambuf: A custom stream buffer that adapts a Python file-like object (with .read())
// into a C++ std::streambuf. This allows C++ code to use std::istream to read data directly
// from Python-managed files, enabling reuse of existing codac2::deserialize(std::istream&, T&) logic.
class py_istreambuf : public std::streambuf
{
  py::object py_file;
  inline static const size_t buffer_size = 4096;
  char buffer[buffer_size];

  public:

    py_istreambuf(py::object py_file_) : py_file(py_file_)
    {
      setg(buffer, buffer, buffer);
    }

  protected:

    int_type underflow() override
    {
      if(gptr() < egptr())
        return traits_type::to_int_type(*gptr());

      py::bytes data_bytes = py_file.attr("read")(buffer_size);
      std::string data_str = data_bytes;

      if(data_str.empty())
        return traits_type::eof();

      std::memcpy(buffer, data_str.data(), data_str.size());
      setg(buffer, buffer, buffer + data_str.size());
      return traits_type::to_int_type(*gptr());
    }
};

void export_serialization(py::module& m)
{
  m.def("serialize", [](py::object py_file, const Vector& x)
    {
      std::ostringstream oss(std::ios::binary);
      serialize(oss, x);
      std::string buffer = oss.str();
      py_file.attr("write")(py::bytes(buffer));
    },
    VOID_SERIALIZE_OSTREAM_REF_CONST_EIGEN_MATRIX_TRC_REF,
    "f"_a, "x"_a);

  m.def("deserialize", [](py::object py_file, Vector& x)
    {
      py_istreambuf buf(py_file);
      std::istream stream(&buf);
      deserialize(stream, x);
    },
    VOID_DESERIALIZE_ISTREAM_REF_EIGEN_MATRIX_TRC_REF,
    "f"_a, "x"_a);
}