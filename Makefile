PLATFORM ?= Linux

ifeq ($(PLATFORM), Windows)
	COMPILER = x86_64-w64-mingw32-g++
	BIN_SUBFOLDER_NAME = mingw
	MAIN_OUT_NAME = main.exe
	CFLAGS = -I$$VCPKG_ROOT/packages/openblas_x64-mingw-dynamic/include/openblas \
			 -L$$VCPKG_ROOT/packages/openblas_x64-mingw-dynamic/lib -lopenblas
	EXTERNAL_DEPENDENCIES = bin/$(BIN_SUBFOLDER_NAME)/libopenblas.dll

else ifeq ($(PLATFORM), Linux)
	COMPILER = g++
	BIN_SUBFOLDER_NAME = linux
	MAIN_OUT_NAME = main
	CFLAGS = -Wl,-rpath='$$ORIGIN/' \
			 -I$$VCPKG_ROOT/packages/openblas_x64-linux-dynamic/include/openblas \
			 -L$$VCPKG_ROOT/packages/openblas_x64-linux-dynamic/lib -lopenblas
	
	VERSIONS := $(wildcard libopenblas.so.[0-9])
	EXTERNAL_DEPENDENCIES = bin/$(BIN_SUBFOLDER_NAME)/libopenblas.so.0 \
						    bin/$(BIN_SUBFOLDER_NAME)/$(VERSIONS)
else
	$(error non-supported platform specified)
endif

ifeq ($(FLAG), debug)
	ifeq ($(PLATFORM), Windows)
		$(error debug not supported for mingw)
	endif
	CFLAGS += -g3
	MAIN_OUT_NAME := $(MAIN_OUT_NAME)_debug
endif

main: src/main.cpp $(EXTERNAL_DEPENDENCIES)
	$(COMPILER) src/main.cpp -static-libstdc++ -o bin/$(BIN_SUBFOLDER_NAME)/$(MAIN_OUT_NAME) $(CFLAGS)

# External lib files
bin/linux/libopenblas.so.0:
	cp $$VCPKG_ROOT/packages/openblas_x64-linux-dynamic/lib/libopenblas.so.0 bin/$(BIN_SUBFOLDER_NAME)/
bin/linux/$(VERSIONS):
	cp $$VCPKG_ROOT/packages/openblas_x64-linux-dynamic/lib/$(VERSIONS) bin/$(BIN_SUBFOLDER_NAME)/
bin/mingw/libopenblas.dll:
	cp $$VCPKG_ROOT/packages/openblas_x64-mingw-dynamic/bin/libopenblas.dll bin/$(BIN_SUBFOLDER_NAME)/