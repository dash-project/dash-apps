include ../make.defs

KMER_LENGTH = 19
KMER_PACKED_LENGTH = $(shell echo $$((($(KMER_LENGTH)+3)/4)))

SRC_EXT = cc
BUILD_PATH = build
BIN_PATH = $(BUILD_PATH)/bin
BIN_NAME = dash

INCLUDES = $(INC) -Iinclude

CXXFLAGS += -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH)

SOURCES = $(shell find src -name '*.$(SRC_EXT)' | sort -k 1nr | cut -f2-)
OBJECTS = $(SOURCES:src/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)

release: CXXFLAGS +=-O2 -Wall -Wextra
release: dirs all

debug: CXXFLAGS += -g -O0 -Wall -Wextra
debug: dirs all

.PHONY: all
all: $(BIN_PATH)/$(BIN_NAME)
	@echo "Making symlink: $(BIN_NAME) -> $<"
	@$(RM) $(BIN_NAME)
	@ln -s $(BIN_PATH)/$(BIN_NAME) $(BIN_NAME)

# Creation of the executable
$(BIN_PATH)/$(BIN_NAME): $(OBJECTS)
	@echo "Linking: $@"
	$(MPICXX) $(OBJECTS) -o $@ $(LIB)

$(BUILD_PATH)/%.o: src/%.$(SRC_EXT)
	@echo "Compiler flags: $(CXXFLAGS)"
	@echo "Compiling: $< -> $@"
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@

.PHONY: clean
clean:
	@echo "Deleting $(BIN_NAME) symlink"
	@$(RM) $(BIN_NAME)
	@echo "Deleting objects and binaries"
	@$(RM) -r $(OBJECTS)
	@$(RM) -r $(BIN_PATH)/$(BIN_NAME)

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJECTS))
	@mkdir -p $(BIN_PATH)

.PHONY: printenv
printenv :
	@echo "MPICXX        = $(MPICXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"
