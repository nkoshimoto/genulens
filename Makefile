BUILD_DIR ?= build
CMAKE     ?= cmake
CTEST     ?= ctest

.PHONY: all genulens pre_gapmoe python test clean configure

all: genulens

configure:
	$(CMAKE) -S . -B $(BUILD_DIR)

genulens: configure
	$(CMAKE) --build $(BUILD_DIR) --target genulens
	cp $(BUILD_DIR)/genulens ./genulens

pre_gapmoe: configure
	$(CMAKE) --build $(BUILD_DIR) --target calc_rho_profile calc_mass_dist calc_murel_dist
	mkdir -p pre_gapmoe
	cp $(BUILD_DIR)/pre_gapmoe/calc_rho_profile pre_gapmoe/calc_rho_profile
	cp $(BUILD_DIR)/pre_gapmoe/calc_mass_dist pre_gapmoe/calc_mass_dist
	cp $(BUILD_DIR)/pre_gapmoe/calc_murel_dist pre_gapmoe/calc_murel_dist

python: configure
	$(CMAKE) --build $(BUILD_DIR) --target genulens_python

test: genulens pre_gapmoe
	$(CMAKE) --build $(BUILD_DIR) --target test_core
	$(CMAKE) --build $(BUILD_DIR) --target genulens_python || true
	$(CTEST) --test-dir $(BUILD_DIR) --output-on-failure

clean:
	rm -rf $(BUILD_DIR)
	rm -f genulens *.o *~
	rm -f pre_gapmoe/*.o pre_gapmoe/calc_rho_profile pre_gapmoe/calc_mass_dist pre_gapmoe/calc_murel_dist
