# Standard platform if not specified as input.
PLATFORM ?= Linux

# Flags used for all compilers.
CFLAGS = -Llib -Iinclude -std=c++20 -static-libstdc++

ifeq ($(PLATFORM), Windows)
	COMPILER = x86_64-w64-mingw32-g++
	BIN_SUBFOLDER_NAME = mingw
	MAIN_OUT_NAME = main.exe
	DYNAMIC_EXTENSION = .dll

else ifeq ($(PLATFORM), Linux)
	COMPILER = g++
	BIN_SUBFOLDER_NAME = linux
	MAIN_OUT_NAME = main
	CFLAGS += -Wl,-rpath='$$ORIGIN'
	DYNAMIC_EXTENSION = .so
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

# Main file for both linux and windows.
bin/$(BIN_SUBFOLDER_NAME)/$(MAIN_OUT_NAME): src/main.cpp lib/libSystemModels$(DYNAMIC_EXTENSION)
	$(COMPILER) $(CFLAGS) src/main.cpp -lSystemModels -o bin/$(BIN_SUBFOLDER_NAME)/$(MAIN_OUT_NAME) \
	-Lexternal/openblas_x64-$(BIN_SUBFOLDER_NAME)-dynamic/lib -lopenblas

# SystemModels windows.
lib/libSystemModels.dll: src/SystemModels.cpp bin/mingw/libopenblas.dll
	$(COMPILER) $(CFLAGS) -fPIC -shared src/SystemModels.cpp -o lib/libSystemModels.dll \
	-Iexternal/openblas_x64-mingw-dynamic/include/openblas -Lexternal/openblas_x64-mingw-dynamic/lib -lopenblas
	cp -a lib/libSystemModels.dll bin/mingw

# SystemModels linux.
lib/libSystemModels.so: src/SystemModels.cpp bin/linux/libopenblas.so.0
	$(COMPILER) $(CFLAGS) -Wl,-soname,libSystemModels.so.0 -fPIC -shared src/SystemModels.cpp -o bin/linux/libSystemModels.so.0.0 \
	-Iexternal/openblas_x64-linux-dynamic/include/openblas -Lexternal/openblas_x64-linux-dynamic/lib -lopenblas
	ln -sfr bin/linux/libSystemModels.so.0.0 bin/linux/libSystemModels.so.0
	ln -sfr bin/linux/libSystemModels.so.0 lib/libSystemModels.so

# External lib files
bin/linux/libopenblas.so.0:
	cp -a external/openblas_x64-linux-dynamic/lib/libopenblas.so.0 external/openblas_x64-linux-dynamic/lib/libopenblas.so.0.3 bin/$(BIN_SUBFOLDER_NAME)/
bin/mingw/libopenblas.dll:
	cp -a /external/openblas_x64-mingw-dynamic/bin/libopenblas.dll bin/$(BIN_SUBFOLDER_NAME)/