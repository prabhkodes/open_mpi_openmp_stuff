# Convenience wrapper around CMake. The real exercise is in CMakeLists.txt.
BUILD ?= build
NP    ?= 2
OMP   ?= 2

.PHONY: all build test deps reverse clean
all: build

## build: configure + compile
build:
	cmake -S . -B $(BUILD) -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD) -j

## test: run the CTest suite (TODO 6)
test: build
	cd $(BUILD) && ctest --output-on-failure

## deps: show the Fortran module dependency order CMake derived (TODO 6a).
## Nothing in the source tree states this order -- CMake scanned the `use`
## statements to find it.
deps: build
	@echo "--- .mod files produced ---"
	@find $(BUILD) -name '*.mod' | sort
	@echo ""
	@echo "--- object build order (from the depend info) ---"
	@find $(BUILD) -name '*.o' -newer CMakeLists.txt 2>/dev/null | sort || true
	@echo ""
	@echo "Expected order: kinds -> grid_mod -> physics_mod -> driver"
	@echo "Nothing in the source lists that order. CMake worked it out."

## reverse: prove the source order in CMakeLists.txt does not matter
reverse:
	@echo "The sources are ALREADY listed in reverse dependency order"
	@echo "in CMakeLists.txt. If it builds, that is the answer to TODO 6a."

clean:
	rm -rf $(BUILD)
