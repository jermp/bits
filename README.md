Bits: Succinct utilities to handle bits
---------------------------------------

A C++ collection of basic data structures, including:

- bitvectors with rank/select queries;
- compact vectors (vectors of fixed-size integers);
- Elias-Fano sequences;
- integer codes (unary, binary, gamma, delta, Rice).


### Inclusion in other projects

### Tests

Tests are written using [doctest](https://github.com/doctest/doctest). To run the tests, first compile the tests:

	mkdir build
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=Debug
	make
	make test
	
To run a specific test, use the tool `-tc` of doctest and specify the test case name. For example:

	./test_elias_fano -tc=access