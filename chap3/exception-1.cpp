#include <exception>
#include <iostream>

sd::exception_ptr eptr;

void foo()
{
	try {
		throw ...;
	}
	catch(...) {
		eptr = std::current_exception();
	}
}

void bar()
{
	if (eptr != nullptr) {
		std::rethrow_exception(eptr);
	}
}
