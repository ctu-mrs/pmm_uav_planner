uniq = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))

CXXFLAGS:=$(call uniq,$(CXXFLAGS))

clean_make: clean $(TARGET)

bin: $(TARGET)

TARGET_OBJS=$(addprefix $(OBJ_DIR)/,$(OBJS))
TEST_TARGET_OBJS=$(addprefix $(OBJ_DIR)/,$(OBJS_TESTS))

$(TEST_TARGET): create_directories $(OBJS) $(OBJS_TESTS)
	echo "making $(TEST_TARGET)"
	$(CXX) -c tests/$(TEST_TARGET).cpp $(CXXFLAGS) $(CPPFLAGS) -o $(OBJ_DIR)/$(TEST_TARGET).o
	$(CXX) $(CXXFLAGS) -o $@ $(TARGET_OBJS) $(TEST_TARGET_OBJS) $(OBJ_DIR)/$(TEST_TARGET).o $(LDFLAGS)

$(OBJS): %.o: src/%.cpp
	echo "making $<"
	$(CXX) -c $< $(CXXFLAGS) $(CPPFLAGS) -o $(OBJ_DIR)/$@

$(OBJS_TESTS): %.o: tests/%.cpp
	echo "making $<"
	$(CXX) -c $< $(CXXFLAGS) $(CPPFLAGS) -o $(OBJ_DIR)/$@

$(TARGET): create_directories $(OBJS) $(TEST_TARGET)
	echo "making $(TARGET)"
	$(CXX) -c src/$(TARGET).cpp $(CXXFLAGS) $(CPPFLAGS) -o $(OBJ_DIR)/$(TARGET).o
	$(CXX) $(CXXFLAGS) -o $@ $(TARGET_OBJS) $(OBJ_DIR)/$(TARGET).o $(LDFLAGS)

create_directories:
	echo "create dircetory $(OBJ_DIR)"
	mkdir -p $(OBJ_DIR)

clean:
	$(RM) $(TARGET)
	$(RM) $(OBJ_DIR)/*
	@for i in $(SUBDIRS) ;\
        do \
        echo "cleaning" all "in $(CURRENT_DIR)/$$i..."; \
        $(MAKE) -C $$i clean; \
        done

clean_all: clean clean_dependencies

