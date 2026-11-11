#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/ADJacobians.h>
#endif

#include <cstdint>
#include <sstream>
#include <utility>
#include <vector>

using namespace MML;
using namespace MML::AD;

namespace {
	void CheckStale(const ADVar& value) {
		CHECK_THROWS_AS(value.value(), ArgumentError);
		CHECK_THROWS_AS(value.adjoint(), ArgumentError);
		CHECK_THROWS_AS(value.index(), ArgumentError);
		CHECK_THROWS_AS(value.tape(), ArgumentError);
		CHECK_THROWS_AS(value.backward(), ArgumentError);
		CHECK_THROWS_AS(+value, ArgumentError);
		CHECK_THROWS_AS(-value, ArgumentError);
	}
}

TEST_CASE("Tape clear rejects old variables constants and outputs after index reuse", "[ad-tape-contract]") {
	Tape tape;
	ADVar x(&tape, 2);
	const auto constant = ADVar::constant(&tape, 7);
	const auto y = x * x;
	tape.clear();
	CHECK(tape.empty());
	CHECK_FALSE(tape.isBackwardComplete());
	ADVar replacement(&tape, 5);
	const auto newConstant = ADVar::constant(&tape, 11);
	const auto result = replacement * replacement;
	CheckStale(x);
	CheckStale(constant);
	CheckStale(y);
	CHECK(replacement.value() == 5);
	CHECK(newConstant.value() == 11);
	result.backward();
	CHECK(replacement.adjoint() == 10);
}

TEST_CASE("Stale recording operands fail before any graph mutation", "[ad-tape-contract]") {
	Tape tape;
	ADVar stale(&tape, 2);
	tape.clear();
	ADVar live(&tape, 5);
	const auto size = tape.size();
	CHECK_THROWS_AS(stale + live, ArgumentError);
	CHECK_THROWS_AS(live - stale, ArgumentError);
	CHECK_THROWS_AS(stale * 3, ArgumentError);
	CHECK_THROWS_AS(3 / stale, ArgumentError);
	CHECK_THROWS_AS(stale / live, ArgumentError);
	CHECK_THROWS_AS(stale / Real(0), ArgumentError);
	CHECK_THROWS_AS(stale /= Real(0), ArgumentError);
	CHECK_THROWS_AS(stale += live, ArgumentError);
	CHECK_THROWS_AS(live += stale, ArgumentError);
	CHECK_THROWS_AS(live *= stale, ArgumentError);
	CHECK_THROWS_AS(live /= stale, ArgumentError);
	CHECK_THROWS_AS(sin(stale), ArgumentError);
	CHECK_THROWS_AS(log(stale), ArgumentError);
	CHECK_THROWS_AS(sqrt(stale), ArgumentError);
	CHECK_THROWS_AS(pow(stale, 0), ArgumentError);
	CHECK_THROWS_AS(pow(2, stale), ArgumentError);
	CHECK_THROWS_AS(ADVar::makeUnary(stale, ADOp::NEG, -2), ArgumentError);
	CHECK_THROWS_AS(ADVar::makeBinary(ADVar(3), stale, ADOp::ADD, 5), ArgumentError);
	CHECK_THROWS_AS(ADVar::makeBinary(stale, ADVar(3), ADOp::ADD, 5), ArgumentError);
	std::ostringstream text;
	CHECK_THROWS_AS(text << stale, ArgumentError);
	CHECK(tape.size() == size);
	CHECK(live.value() == 5);
}

TEST_CASE("Cleared empty tapes reject handles before any raw entry access", "[ad-tape-contract]") {
	Tape tape;
	ADVar old(&tape, 2);
	tape.clear();
	// The legacy index accessor is safe even when its entry no longer exists.
	REQUIRE_THROWS_AS(old.index(), ArgumentError);
	CheckStale(old);
	CHECK_THROWS_AS(ADVar(nullptr, 2), ArgumentError);
}

TEST_CASE("Tape move construction invalidates source handles even after source refill", "[ad-tape-contract]") {
	Tape source;
	ADVar x(&source, 3);
	const auto y = x * x;
	const int outputIndex = y.index();
	y.backward();
	Tape moved(std::move(source));
	CHECK(source.empty());
	CHECK_FALSE(source.isBackwardComplete());
	CHECK(moved.isBackwardComplete());
	ADVar fresh(&source, 7);
	const auto freshOutput = fresh * fresh;
	CheckStale(x);
	CheckStale(y);
	CHECK(moved[0].value == 3);
	CHECK(moved[outputIndex].value == 9);
	moved.backward(outputIndex);
	CHECK(moved[0].adjoint == 6);
	freshOutput.backward();
	CHECK(fresh.adjoint() == 14);
}

TEST_CASE("Tape copy assignment invalidates destination handles and preserves source handles", "[ad-tape-contract]") {
	Tape source;
	ADVar x(&source, 3);
	const auto y = x * x;
	y.backward();
	Tape destination;
	ADVar oldInput(&destination, 7);
	const auto oldOutput = oldInput * oldInput;
	destination = source;
	CHECK(destination.isBackwardComplete());
	CheckStale(oldInput);
	CheckStale(oldOutput);
	CHECK(x.value() == 3);
	CHECK(y.value() == 9);
	CHECK(x.adjoint() == 6);
	destination.clear();
	CHECK(y.value() == 9);
	y.backward();
	CHECK(x.adjoint() == 6);
}

TEST_CASE("Tape move assignment invalidates both source and replaced destination handles", "[ad-tape-contract]") {
	Tape source;
	ADVar x(&source, 3);
	const auto y = x * x;
	const int outputIndex = y.index();
	y.backward();
	Tape destination;
	ADVar oldInput(&destination, 7);
	const auto oldOutput = oldInput * oldInput;
	destination = std::move(source);
	CHECK(destination.isBackwardComplete());
	CHECK(source.empty());
	CHECK_FALSE(source.isBackwardComplete());
	ADVar fresh(&source, 8);
	const auto freshOutput = fresh * fresh;
	CheckStale(x);
	CheckStale(y);
	CheckStale(oldInput);
	CheckStale(oldOutput);
	destination.backward(outputIndex);
	CHECK(destination[0].adjoint == 6);
	freshOutput.backward();
	CHECK(fresh.adjoint() == 16);
}

TEST_CASE("Tape copies own independent graph storage without invalidating the original", "[ad-tape-contract]") {
	Tape source;
	ADVar x(&source, 3);
	const auto y = x * x;
	Tape copy(source);
	CHECK(copy.size() == source.size());
	CHECK_FALSE(copy.isBackwardComplete());
	copy.backward(y.index());
	CHECK(copy[0].adjoint == 6);
	CHECK(x.adjoint() == 0);
	CHECK_FALSE(source.isBackwardComplete());
	ADVar extra(&copy, 7);
	copy.clear();
	CHECK(x.value() == 3);
	y.backward();
	CHECK(x.adjoint() == 6);
}

TEST_CASE("Tape self copy and self move preserve handles graph and completed sweep", "[ad-tape-contract]") {
	Tape tape;
	ADVar x(&tape, 3);
	const auto y = x * x;
	y.backward();
	const auto size = tape.size();
	Tape* same = &tape;
	tape = *same;
	CHECK(tape.size() == size);
	CHECK(tape.isBackwardComplete());
	CHECK(x.adjoint() == 6);
	tape = std::move(*same);
	REQUIRE(tape.size() == size);
	CHECK(tape.isBackwardComplete());
	CHECK(x.value() == 3);
	CHECK(y.value() == 9);
	y.backward();
	CHECK(x.adjoint() == 6);
}

TEST_CASE("Ordinary tape reallocation preserves variables outputs and copied handles", "[ad-tape-contract]") {
	Tape tape;
	ADVar x(&tape, 2);
	const auto copied = x;
	const auto moved = std::move(x);
	const auto square = x * x;
	auto output = square;
	bool reallocated = false;
	for (int i = 0; i < 512; ++i) {
		const auto address = reinterpret_cast<std::uintptr_t>(&tape[0]);
		output += x;
		reallocated = reallocated || address != reinterpret_cast<std::uintptr_t>(&tape[0]);
	}
	CHECK(reallocated);
	CHECK(x.value() == 2);
	CHECK(copied.value() == 2);
	CHECK(moved.value() == 2);
	CHECK(square.value() == 4);
	CHECK(output.value() == 4 + 2 * 512);
	CHECK(x.tape() == &tape);
	output.backward();
	CHECK(copied.adjoint() == 4 + 512);
	CHECK(moved.adjoint() == x.adjoint());
	output.backward();
	CHECK(x.adjoint() == 4 + 512);
}

TEST_CASE("Repeated clear and replacement never resurrect an old matching index", "[ad-tape-contract]") {
	Tape tape;
	ADVar old(&tape, 2);
	for (int i = 0; i < 8; ++i) {
		tape.clear();
		ADVar fresh(&tape, Real(i + 3));
		CHECK_THROWS_AS(old.value(), ArgumentError);
		Tape incoming;
		ADVar input(&incoming, Real(i + 7));
		tape = std::move(incoming);
		CHECK_THROWS_AS(old.value(), ArgumentError);
		CHECK_THROWS_AS(fresh.value(), ArgumentError);
		CHECK(tape[0].value == i + 7);
	}
}

TEST_CASE("Reverse helpers reject callbacks that recycle recorded input handles", "[ad-tape-contract]") {
	const auto stale = [](auto& input) {
		Tape* tape = input[0].tape();
		const auto old = input[0];
		tape->clear();
		ADVar replacement(tape, 5);
		return old;
	};
	const auto rebound = [](auto& input) {
		Tape* tape = input[0].tape();
		tape->clear();
		ADVar replacement(tape, 5);
		return replacement * replacement;
	};
	CHECK_THROWS_AS(CalcGradientReverseAD(stale, std::vector<Real>{2}), ArgumentError);
	CHECK_THROWS_AS(CalcGradientReverseAD(rebound, std::vector<Real>{2}), ArgumentError);
	CHECK_THROWS_AS(CalcJacobianReverseAD([&](auto& input) { return std::vector<ADVar>{stale(input)}; }, std::vector<Real>{2}), ArgumentError);
	CHECK_THROWS_AS(CalcJacobianReverseAD([&](auto& input) { return std::vector<ADVar>{rebound(input)}; }, Vector<Real>{2}), ArgumentError);
}

TEST_CASE("Untracked constants remain independent of every tape lifetime transition", "[ad-tape-contract]") {
	const ADVar constant(7);
	const auto alsoConstant = ADVar::constant(nullptr, 11);
	Tape tape;
	ADVar input(&tape, 2);
	tape.clear();
	Tape other(std::move(tape));
	tape = other;
	CHECK(constant.value() == 7);
	CHECK(constant.adjoint() == 0);
	CHECK(constant.index() == -1);
	CHECK(constant.tape() == nullptr);
	CHECK(alsoConstant.value() == 11);
	CHECK_THROWS_AS(constant.backward(), ArgumentError);
	const auto unusedInputs = [](auto& input) {
		input[0].tape()->clear();
		return ADVar(7);
	};
	const auto gradient = CalcGradientReverseAD(unusedInputs, std::vector<Real>{2});
	REQUIRE(gradient.size() == 1);
	CHECK(gradient[0] == 0);
	const auto jacobian = CalcJacobianReverseAD([&](auto& values) { return std::vector<ADVar>{unusedInputs(values)}; }, std::vector<Real>{2});
	CHECK(jacobian(0, 0) == 0);
}
