CXX := g++
CXXFLAGS := -Wno-parentheses -Wall -Wextra -std=c++20 -Iinclude -Iexternal
LDFLAGS = -lboost_unit_test_framework


SRCDIR := src src/user
SRCS = $(wildcard $(patsubst %,%/*.cpp,$(SRCDIR)))
# Collect targets
TARGETS = $(SRCS)
TARGETS := $(subst src/,,$(TARGETS))           # Remove src/ 
TARGETS := $(subst user/,user_,$(TARGETS))     # user/ -> user_
TARGETS := $(patsubst %.cpp, bin/%, $(TARGETS)) 

TEST_SRCS = $(wildcard tests/*.cpp)
TEST_TARGETS = $(patsubst tests/%.cpp, bin/tests/%, $(TEST_SRCS))

$(info SOURCES = $(SRCS))
$(info TARGETS = $(TARGETS))

all: $(TARGETS) $(TEST_TARGETS)

# Build default
bin/%: src/%.cpp | bin
	$(CXX) $(CXXFLAGS) -o $@ $<

# Build untracked user projects
bin/user_%: src/user/%.cpp | bin
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/tests/%: tests/%.cpp | bin/tests
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< 
bin:
	mkdir -p bin/tests bin/user

clean:
	rm -rf bin build

# Select target
run: bin/$(TARGET)
	./bin/$(TARGET) $(ARGS)

test: $(TEST_TARGETS)
	@for test in $(TEST_TARGETS); do \
		echo "Running $$test"; \
		$$test; \
	done

.PHONY: all clean run test

