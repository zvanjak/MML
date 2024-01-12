# Scalar and vector field operations

- grad
- div
- curl

~~~c++
TEST_CASE("Matrix_default_ctor_init_to_zero", "[simple]") {
    MML::Matrix a(2,2);

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}
~~~

