# Development Tasks: Distributional Regression with Time-Varying Parameters

## Prerequisites
- Read CLAUDE.md thoroughly for package standards and workflows
- Use `devtools::load_all()` before testing any changes
- **MANDATORY**: Use code-reviewer agent for ALL proposed R code changes BEFORE making edits
- Run `Rscript -e "devtools::check()"` after significant changes
- Follow one sub-task at a time workflow: complete task → mark as done → ask "Ready for next sub-task?" → wait for approval

## Overview
After analyzing mgcv, gratia, and ffc implementations, this revised plan leverages existing multi-parameter infrastructure rather than rebuilding it. The focus is on extending ffc's `fts()` syntax and forecasting to work seamlessly with mgcv's distributional families.

## Key Insights from Package Analysis
- **mgcv already handles list formulae** via `family$nlp`, storing multiple linear predictors in `pred.formula` and `lpi`
- **gratia provides multi-parameter sampling** through `fitted_samples()` and parameter-specific predictions using `lpmat`
- **ffc's modular design** separates formula interpretation (`interpret_ffc()`), basis evaluation (`basis_fn()`), and model fitting (`ffc_gam()`)

## Core Challenge
Enable `fts()` syntax within list formulae for distributional regression:
```r
# Target syntax:
ffc_gam(list(y ~ fts(time, k=10), ~ fts(time, k=5)), 
        data=data, family=mgcv::gaulss())
```

---

## Task 1: Formula Interpretation Extension (60 min total)

### 1.1: List Formula Detection and Validation (15 min) ✅ **COMPLETED**
**Goal**: Add validation-only logic to `ffc_gam()` for list formulae compatibility

**Tasks**:
1. ✅ Add validation after line 94 in `ffc_gam()` to check list formula/family compatibility
2. ✅ Validate list length matches `family$nlp` if multi-parameter family  
3. ✅ Add informative error messages for mismatched list/family combinations
4. ✅ Add warning for single formula with multi-parameter family
5. ✅ **Note**: Actual list formula processing will be handled in Task 1.2 via `interpret_ffc()` extension

**Testing**:
- ✅ Test validation catches family without `nlp` property
- ✅ Test validation catches length mismatches  
- ✅ Test warning for single formula with multi-parameter family
- ✅ Verify single formulae still work unchanged

**Files**:
- ✅ Modified: `R/ffc_gam.R` (validation logic only)
- ✅ Updated: `tests/testthat/test-list-formulae.R` (fast unit tests)

### 1.2: Multi-Parameter interpret_ffc() (15 min) ✅ **COMPLETED**
**Goal**: Extend `interpret_ffc()` to natively handle list formulae containing `fts()` terms

**Tasks**:
1. ✅ Modify `interpret_ffc()` in `R/interpret_ffc.R` to detect and process list formulae
2. ✅ Process each formula element separately for `fts()` terms within single function call
3. ✅ Return same structure as single formula but with list of interpreted formulae
4. ✅ Preserve parameter-specific smoothing specifications and coordinate data consistently
5. ✅ Ensure `ffc_gam()` can use returned structure directly without modification

**Testing**:
- ✅ Test with list formulae containing `fts()` terms
- ✅ Verify each parameter gets correct smooth specifications  
- ✅ Check non-`fts()` terms pass through unchanged
- ✅ Test integration with `ffc_gam()` end-to-end

**Files**:
- ✅ Modified: `R/interpret_ffc.R` (in-place DRY implementation)
- ✅ Updated: `tests/testthat/test-list-formulae.R` (existing tests sufficient)

### 1.3: Basis Function Coordination (15 min) ✅ **COMPLETED**
**Goal**: Ensure basis functions are correctly set up for each parameter

**Tasks**:
1. ✅ Verify basis setup functions handle list of formulae (already implemented)
2. ✅ Ensure each parameter gets appropriate basis functions evaluated (already implemented)  
3. ✅ Handle different `k` values per parameter (already implemented via independent processing)
4. ✅ Add tests to verify different `k` per parameter works correctly

**Testing**:
- ✅ Test basis evaluation with different `k` per parameter
- ✅ Verify each parameter processes independently
- ✅ Check basis functions have correct dimensions

**Files**:
- ✅ Verified: `R/interpret_ffc.R` (existing list formula handling sufficient)
- ✅ Updated: `tests/testthat/test-list-formulae.R` (tests added)

### 1.4: Model Object Structure Preservation (15 min) ✅ **COMPLETED**
**Goal**: Ensure ffc_gam objects preserve mgcv's multi-parameter structure

**Tasks**:
1. ✅ Verify `ffc_gam()` preserves `pred.formula` and `lpi` from mgcv (already preserved automatically)
2. ✅ Parameter-specific metadata accessible via mgcv's existing properties (`family$nlp`, `family$names`)
3. ✅ Class `ffc_gam_multi` already added for distributional models (lines 209-214)
4. ✅ Backward compatibility maintained for single-parameter models

**Testing**:
- ✅ Multi-parameter models retain mgcv structure (mgcv handles automatically)
- ✅ Parameter metadata accessible via existing mgcv properties
- ✅ Single-parameter models unchanged

**Files**:
- ✅ Verified: `R/ffc_gam.R` (existing implementation sufficient - mgcv preserves all required properties)

---

## Task 2: Time-Varying Coefficient Extraction (75 min total)

### 2.1: Multi-Parameter fts() Method (15 min) ✅ **COMPLETED**
**Goal**: Add parameter-aware coefficient naming for distributional models

**Tasks**:
1. ✅ Add `parameter_id` argument to `dyn_to_spline()` function
2. ✅ Update coefficient naming to include parameter prefix for distributional models
3. ✅ Pass parameter index from `interpret_ffc()` list formula processing
4. ✅ Maintain backward compatibility for single-parameter models

**Testing**:
- ✅ Coefficient names now include parameter information (e.g., "param1_fts_bs_...")
- ✅ Single-parameter models unchanged (parameter_id = NULL)
- ✅ Users can distinguish coefficients by parameter for identical predictors

**Files**:
- ✅ Modified: `R/interpret_ffc.R` (added parameter_id support and naming logic)

### 2.2: Parameter-Specific Coefficient Utilities (15 min) ✅ **COMPLETED**
**Goal**: Provide access to parameter-specific coefficients

**Tasks**:
1. ✅ Parameter-prefixed coefficient names from Task 2.1 provide user-friendly access
2. ✅ mgcv's `pred.formula` provides parameter names
3. ✅ mgcv's `lpi` provides coefficient indices mapping
4. ✅ No additional utility functions needed - existing infrastructure sufficient

**Testing**:
- ✅ Parameter names accessible via `names(object$pred.formula)`
- ✅ Coefficient indices accessible via `object$lpi`
- ✅ Parameter-specific coefficients identifiable by name prefix

**Files**:
- ✅ No new files needed (existing mgcv structure + Task 2.1 naming sufficient)

### 2.3: Time Index Coordination (15 min) ✅ **COMPLETED**
**Goal**: Determine if time coordination is needed across parameters

**Tasks**:
1. ✅ Analyzed time coordination requirements - not needed for distributional models
2. ✅ Different parameters (mu, sigma) should naturally evolve at different rates
3. ✅ Independent parameter processing is correct and beneficial design
4. ✅ No time coordination constraints needed - parameters should be flexible

**Testing**:
- ✅ Each parameter processes independently as designed
- ✅ No validation needed for mismatched time specifications (this is desirable)
- ✅ Flexibility maintained for different temporal resolutions per parameter

**Files**:
- ✅ No changes needed (current independent processing is optimal)

### 2.4: fts Object Enhancement (15 min) ✅ **COMPLETED**
**Goal**: Enhance fts objects to carry parameter information

**Tasks**:
1. ✅ Added `.parameter` column to `fts_ts` objects with standardized parameter names
2. ✅ Created `extract_parameter_from_basis()` helper function with clean naming system
3. ✅ Standardized parameter names: "location", "scale", "shape" for all distributional families
4. ✅ Full backward compatibility - single-parameter models show "location"

**Testing**:
- ✅ Test parameter extraction for distributional families (location, scale)
- ✅ Test single-parameter model fallback to "location"
- ✅ Test invalid parameter number handling  
- ✅ All 40 tests passing with standardized naming

**Files**:
- ✅ Modified: `R/fts_coefs.R` (added parameter column with location/scale/shape naming)
- ✅ Modified: `tests/testthat/test-list-formulae.R` (comprehensive parameter extraction tests)

### 2.5: Integration Testing (15 min) ✅ **COMPLETED**
**Goal**: Test full coefficient extraction pipeline

**Tasks**:
1. ✅ Test end-to-end coefficient extraction with gaulss()
2. ✅ Verify twlss() and betar() work correctly
3. ✅ Test with different `fts()` specifications per parameter
4. ✅ Check integration with existing fts methods

**Testing**:
- ✅ Run integration tests with multiple families
- ✅ Test different fts specifications
- ✅ Verify no regressions in single-parameter models

**Files**:
- ✅ Updated: `tests/testthat/test-list-formulae.R` (comprehensive integration tests added)

### 2.6: Fix Core Infrastructure Issues (30 min) ✅ **COMPLETED**
**Goal**: Resolve dimension mismatch errors in distributional model fitting

**Tasks**:
1. ✅ Debug and identify root cause of `data.frame(...): arguments imply differing number of rows` errors
2. ✅ Fix `update_mod_data()` function to handle list formulae correctly
3. ✅ Replace `terms.formula()` term labels with `all.vars()` for variable extraction  
4. ✅ Add comprehensive helper functions for formula handling
5. ✅ Ensure backward compatibility with single-parameter models

**Testing**:
- ✅ Verify distributional models can be fitted successfully
- ✅ Test dimension compatibility in data frame operations
- ✅ Confirm no regressions in existing functionality

**Files**:
- ✅ Modified: `R/ffc_gam.R` (fixed `update_mod_data()` function)
- ✅ Modified: `R/utils.R` (added `validate_formula_input()` and `get_primary_formula()` helpers)
- ✅ Modified: `R/validations.R` (updated `convert_re_to_factors()` and `validate_response_in_data()`)

### 2.7: Full Test Suite Validation (20 min) ✅ **COMPLETED**
**Goal**: Run complete test suite and resolve any remaining failures or warnings

**Tasks**:
1. ✅ Use r-test-runner agent to execute full test suite (`devtools::test()`)
2. ✅ Analyze all test failures, warnings, and errors comprehensively
3. ✅ Fix any remaining issues with distributional regression integration
4. ✅ Address coefficient extraction issues for models without `fts()` terms
5. ✅ Ensure all tests pass with no warnings

**Testing**:
- ✅ Run complete test suite across all test files
- ✅ Verify all existing functionality remains intact
- ✅ Check performance and memory usage acceptable
- ✅ Confirm package loads and builds cleanly

**Files**:
- ✅ Modified: `R/fts_coefs.R` (fixed pattern matching from `:fts_` to `fts_`)
- ✅ Modified: `R/interpret_ffc.R` (added semantic naming with `get_parameter_prefix()` helper)
- ✅ Updated: `tests/testthat/test-list-formulae.R` (updated all tests for semantic naming)

---

## Task 3: Forecasting Pipeline Extension (90 min total)

### 3.1: Multi-Parameter Forecast Detection (15 min)
**Goal**: Extend `forecast.ffc_gam()` to detect and route multi-parameter models

**Tasks**:
1. Modify `forecast.ffc_gam()` in `R/forecast.R` to detect distributional families
2. Add parameter detection using utilities from Task 2.2
3. Route to appropriate forecasting method based on family type
4. Maintain existing interface and arguments

**Testing**:
- Test detection of multi-parameter families
- Verify correct routing to forecasting methods
- Check backward compatibility

**Files**:
- Modify: `R/forecast.R`
- Modify: `tests/testthat/test-forecast.R`

### 3.2: Parameter-Specific Time Series Forecasting (15 min)
**Goal**: Forecast time-varying coefficients for each parameter separately

**Tasks**:
1. Add `forecast_parameter_coefficients()` function
2. Extract fts object for specific parameter
3. Apply time series forecasting to parameter coefficients
4. Return forecasted coefficients with uncertainty

**Testing**:
- Test forecasting for individual parameters
- Verify uncertainty propagation
- Check forecasted coefficient structure

**Files**:
- Modify: `R/forecast.R`
- Create: `tests/testthat/test-parameter-forecasting.R`

---

## Task 4: Enhanced Prediction Interface (60 min total)

### 4.1: Multi-Parameter Prediction Method (15 min)
**Goal**: Extend `predict.ffc_gam()` to support parameter selection

**Tasks**:
1. Add `parameter` argument to `predict.ffc_gam()`
2. Use mgcv's `lpmat` for parameter-specific predictions
3. Return list of predictions when `parameter=NULL`
4. Maintain backward compatibility

**Testing**:
- Test parameter-specific predictions
- Test prediction list output
- Verify mgcv lpmat integration

**Files**:
- Modify: `R/predict.R`
- Modify: `tests/testthat/test-predict.R`


### 4.2: Output Format Standardization (15 min)
**Goal**: Create consistent output formats across functions

**Tasks**:
1. Add output formatting utilities
2. Support both wide and long format outputs
3. Add proper class attributes and metadata
4. Ensure consistency across predict/forecast functions

**Testing**:
- Test output format consistency
- Verify metadata preservation
- Check format conversion utilities

**Files**:
- Modify: `R/distributional_utils.R`
- Modify: `tests/testthat/test-distributional-utils.R`

---

## Task 5: Documentation and Testing (75 min total)

### 5.1: Function Documentation (15 min)
**Goal**: Add comprehensive roxygen2 documentation

**Tasks**:
1. Document all modified functions with roxygen2
2. Add distributional regression examples
3. Update existing documentation for new parameters
4. Include cross-references between related functions

**Testing**:
- Build documentation and check examples
- Verify cross-references work
- Test examples run successfully

**Files**:
- All modified R files

### 5.2: Integration Testing Suite (15 min)
**Goal**: Create comprehensive integration tests

**Tasks**:
1. Test complete workflow with gaulss() family
2. Test with twlss() and betar() families
3. Test forecasting pipeline end-to-end
4. Verify no regressions in existing functionality

**Testing**:
- Run complete integration tests
- Check all existing tests pass
- Verify performance acceptable

**Files**:
- Create: `tests/testthat/test-distributional-integration-complete.R`

### 5.3: Error Handling and Edge Cases (15 min)
**Goal**: Add robust error handling throughout

**Tasks**:
1. Add validation for list formulae and family compatibility
2. Handle edge cases in coefficient extraction
3. Add informative error messages for common mistakes
4. Test error conditions comprehensively

**Testing**:
- Test error conditions and messages
- Verify graceful handling of edge cases
- Check validation catches problems early

**Files**:
- All modified R files
- Enhanced test files

### 5.4: Performance Testing (15 min)
**Goal**: Ensure acceptable performance with distributional models

**Tasks**:
1. Benchmark fitting time vs single-parameter models
2. Test memory usage with large datasets
3. Profile forecasting performance
4. Optimize bottlenecks if needed

**Testing**:
- Run performance benchmarks
- Check memory usage reasonable
- Verify scalability

**Files**:
- Create: `tests/testthat/test-performance.R`

### 5.5: Vignette Creation (15 min)
**Goal**: Create tutorial for distributional regression

**Tasks**:
1. Create `vignettes/distributional-regression.Rmd`
2. Include financial volatility example with gaulss()
3. Show environmental data example with twlss()
4. Demonstrate parameter-specific forecasting

**Testing**:
- Build vignette successfully
- Verify examples run and produce sensible results
- Check educational value

**Files**:
- Create: `vignettes/distributional-regression.Rmd`

---

## Implementation Priority Order

1. **Task 1 (Formula Interpretation)** - Foundation for everything else
2. **Task 2 (Coefficient Extraction)** - Essential for time-varying functionality  
3. **Task 3 (Forecasting)** - Core new functionality
4. **Task 4 (Prediction)** - User interface completion
5. **Task 5 (Documentation/Testing)** - Quality assurance

## Final Validation Checklist

- [ ] All tests pass: `devtools::test()`
- [ ] Package check clean: `devtools::check()`
- [ ] Examples and vignettes build successfully
- [ ] No breaking changes to existing functionality
- [ ] Performance acceptable for typical use cases
- [ ] Code reviewed and approved

---

## Key Success Metrics

1. **Functional**: Can fit `ffc_gam(list(y ~ fts(time), ~ fts(time)), family=gaulss())`
2. **Extraction**: Can extract parameter-specific coefficients with `fts(model, parameter="mu")`
3. **Forecasting**: Can forecast multi-parameter models with `forecast(model, h=10)`
4. **Integration**: Existing single-parameter functionality unchanged
5. **Documentation**: Clear examples and tutorial available

## Relevant Files Modified/Created

### Core R Files
- `R/ffc_gam.R` - List formula handling, multi-parameter detection
- `R/interpret_ffc.R` - Multi-parameter formula interpretation  
- `R/fts.R` - Parameter-specific coefficient extraction
- `R/forecast.R` - Multi-parameter forecasting pipeline
- `R/predict.R` - Parameter-specific prediction methods
- `R/distributional_utils.R` - New utilities for multi-parameter support
- `R/posterior_samples.R` - New posterior sampling methods

### Test Files  
- `tests/testthat/test-distributional-*.R` - Comprehensive test suite
- Modified existing test files for backward compatibility

### Documentation
- `vignettes/distributional-regression.Rmd` - Tutorial vignette
- Updated man pages for all modified functions
