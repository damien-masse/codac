# -*- coding: utf-8 -*-
""" doxygen2docstring.py

  This script converts XML file generated by Doxygen (C++)
  into docstring format (Python). The output is a C++ file
  containing char* strings of docstring values, available as
  constants. They can be used with pybind11 for sharing C++
  documentation to Python binding.
  
  Project: codac
  Author: Simon Rohou
  Date: 2020
"""

import os
import sys
import re
import xml.etree.ElementTree as ET

def normalize_template_label(str_template_label):

  return str_template_label \
    .replace("< ", "_") \
    .replace("<", "_") \
    .replace(", ", "") \
    .replace(",", "") \
    .replace(" >", "") \
    .replace(">", "")

def normalize_label(str_label):

  str_label = str_label \
    .replace("or<<", "OROUT") \
    .replace("or>>", "ORIN")

  return normalize_template_label(str_label) \
    .replace("constexpr ", "") \
    .replace("codac::", "") \
    .replace("codac2::", "") \
    .replace("std::", "") \
    .replace("const ", "const_") \
    .replace(" ", "_") \
    .replace("::", "_") \
    .replace("\"\"_", "LITT") \
    .replace("~", "TILD_") \
    .replace("operator+", "OPERATORPLUS") \
    .replace("operator-", "OPERATORMINUS") \
    .replace("operator*", "OPERATORMUL") \
    .replace("operator/", "OPERATORDIV") \
    .replace("operator[]", "OPERATORCOMPO") \
    .replace("operator()", "OPERATORCALL") \
    .replace("operator&=", "OPERATORANDEQ") \
    .replace("operator|=", "OPERATOROREQ") \
    .replace("operator&", "OPERATORAND") \
    .replace("operator|", "OPERATOROR") \
    .replace("operator==", "OPERATOREQ") \
    .replace("operator!=", "OPERATORNEQ") \
    .replace("operator=", "OPERATORAFF") \
    .replace("&", "_REF") \
    .replace("*", "_PTR") \
    .replace("=", "EQ") \
    .replace("size_t", "SIZET") \
    .replace("...", "_VARIADIC") \
    .replace("(", "_") \
    .replace(")", "_") \
    .replace("__REF", "_REF") \
    .replace("__PTR", "_PTR") \
    .replace("-1", "MINUSONE")

def docstring_varname(memberdef):

  txt = normalize_label(memberdef.find(".//definition").text)

  for param in memberdef.findall("param"):
    if param.find("type") != None:
      txt += "_" + normalize_label("".join(param.find("type").itertext()))

  if memberdef.get('const') != None and memberdef.attrib['const'] == "yes":
    txt += "_const"

  return txt.upper()


def get_originate_file(m):

    location_file = m.find(".//location").get("file").split("/")[-1]
    return location_file.replace(".h", "_docs.h").replace("codac2_", "codac2_py_")


if not os.path.exists(sys.argv[1]):
  print("-- /!\\ Enable to build doc files (generate XML Doxygen files first)")
  quit()

files = os.listdir(sys.argv[1])
for xml_doc in files:

  # Parsing XML file from Doxygen output (XML mode enabled)
  tree = ET.parse(sys.argv[1] + "/" + xml_doc)
  root = tree.getroot()


  # Main documentation (ex: comment of the class)

  class_compound = root.find("./compounddef")

  if class_compound == None or (class_compound.get('kind') != None and \
    (class_compound.attrib['kind'] == "file" or class_compound.attrib['kind'] == "dir")):
    continue

  class_name = class_compound.find("./compoundname").text

  with open(sys.argv[2] + "/" + get_originate_file(class_compound), 'a', encoding='utf-8') as f:

    print("/// Class " + class_name, file=f)
    #print("const char* "
    #   + normalize_label(class_name + "_MAIN").upper()
    #   + " = R\"Docstring documentation will be available in next release.\";", file=f)
    print("#define "
       + normalize_label(class_name + "_MAIN").upper()
       + " \"Docstring documentation will be available in next release.\"", file=f)


  # Members documentation (ex: methods, functions, etc.)

  for memberdef in root.findall(".//memberdef"):

    with open(sys.argv[2] + "/" + get_originate_file(memberdef), 'a', encoding='utf-8') as f:

      if memberdef.find("./definition") == None:
        continue

      # Doc string
      print("// " + memberdef.find("./definition").text, end="", file=f)
      argstring = memberdef.find("./argsstring").text
      if argstring:
        print(argstring, end="", file=f)
      print("", file=f)

      #if memberdef.get('id') != None and (memberdef.attrib['id'] == "codac2___interval__operations_8h_1a64ecd476860d3067b13595c4b3022edd" or memberdef.attrib['id'] == "classcodasc2_1_1_interval_matrix_template___1a8b52d938a0484425639aa86603b03009"):
      #  print("start")
      #  print("const char* " + docstring_varname(memberdef) + " = \"Docstring documentation will be available in next release.\"")
      #  exit(1)

      if memberdef.get('kind') != None and memberdef.attrib['kind'] == "friend" \
          and memberdef.get('inline') != None and memberdef.attrib['inline'] != "yes":
        print("// Friend member is not documented here.\n", file=f)
        continue

      #print("const char* "
      #   + docstring_varname(memberdef)
      #   + " = R\"Docstring documentation will be available in next release.\";", file=f)

      print("#define "
         + docstring_varname(memberdef)
         + " \"Docstring documentation will be available in next release.\"", file=f)

      print("\n", file=f)