CXX := g++
CXXFLAGS := -Wno-parentheses -Wall -Wextra -std=c++20 -Iinclude -Iexternal

SRCDIR := src src/user

SRCS = $(wildcard $(patsubst %,%/*.cpp,$(SRCDIR)))
# Collect targets
TARGETS = $(SRCS)
TARGETS := $(subst src/,,$(TARGETS))           # Remove src/ 
TARGETS := $(subst user/,user_,$(TARGETS))     # user/ -> user_
TARGETS := $(patsubst %.cpp, bin/%, $(TARGETS)) 

$(info SOURCES = $(SRCS))
$(info TARGETS = $(TARGETS))

all: $(TARGETS) 

# Build
bin/%: src/%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/user_%: src/user/%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

bin:
	mkdir -p bin

clean:
	rm -rf bin build

# Select target
run: bin/$(TARGET)
	./bin/$(TARGET) $(ARGS)

.PHONY: all clean run

