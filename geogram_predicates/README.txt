This is Predicates Pluggable Software Module, extracted
from GEOGRAM source tree. It contains a standalone .cpp/.h
pair that can be used in any program and that does not have
any dependency. 

It may also contain an example program that can be compiled by using:
  g++ --std=c++11 Predicates_psm.cpp Predicates_example.cpp -o Predicates_example
(or gcc if it is plain C, as in OpenNL)

Some examples may require additional compilation flags (see comments at the beginning
of the example source, e.g. Delaunay_example.cpp).
