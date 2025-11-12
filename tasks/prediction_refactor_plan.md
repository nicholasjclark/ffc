# Prediction Functionality Refactoring Plan

## Overview
This plan addresses critical issues identified by devtools::check() and agent analysis to improve prediction code efficiency, clarity, scalability, and DRY principles. Each step is designed to be completed in 30 minutes by junior developers.

## Prerequisites
- Read CLAUDE.md thoroughly
- Use `devtools::load_all()` before testing
- Use code-reviewer agent to approve all major edits
- Run `Rscript -e "devtools::check()"` after each step

---

## Step 1: Fix Package Check Errors and Missing Imports (30 min)

### Goal
Resolve the critical ERROR and NOTE from R CMD check that prevent package validation.

### Issues to Fix
1. **ERROR**: Missing `model()` function in ARDF example
2. **NOTE**: Missing imports for `tail`, `key_vars`, `index_var`
3. **NOTE**: Undefined global variables `.series`, `.y`

### Tasks
1. Fix ARDF example in `R/ardf.R` lines 80-85 - replace broken `model()` call with working example
2. Add missing imports to `R/re-exports.R` for utils and tsibble functions
3. Fix global variable bindings in `R/stan_dataprep.R` by adding NSE declarations

### Testing
- Run `devtools::check()` and verify no ERROR status
- Ensure examples execute without errors

### Files Modified
- `R/ardf.R` (example fix)
- `R/re-exports.R` (add imports)
- `R/stan_dataprep.R` (global variable fix)

---

## Step 2: Create Configuration Management System (30 min)

### Goal
Eliminate hardcoded magic numbers by centralizing configuration defaults.

### Current Problem
- `times = 200` repeated in 5+ locations across forecast code
- Stan parameters duplicated in 3 model training functions
- `time_k = 10` scattered throughout multiple files

### Tasks
1. Create new `R/config.R` file with centralized defaults for Stan parameters, simulation settings, and forecast options
2. Add helper function to retrieve configuration values safely
3. Replace hardcoded values in `R/forecast.R` with configuration calls
4. Update Stan model functions to use centralized defaults

### Testing
- Load package and verify all defaults work correctly
- Ensure no behavioral changes in existing functionality

### Files Modified
- `R/config.R` (new file)
- `R/forecast.R` (replace hardcoded values)
- `R/ardf.R`, `R/vardf.R`, `R/gpdf.R` (use config defaults)

---

## Step 3: Optimize Matrix Operations (30 min)

### Goal
Replace inefficient nested apply/transpose patterns with vectorized matrix operations for 50-70% speed improvement.

### Current Problem
Inefficient matrix operations in `forecast.ffc_gam()` lines 317-327 using nested apply() calls and repeated transpositions.

### Tasks
1. Create optimized matrix computation function in `R/utils.R`
2. Replace the inefficient apply/t pattern with vectorized matrix multiplication
3. Handle offset addition efficiently using sweep() or similar vectorized operations
4. Test performance improvement with benchmark

### Testing
- Compare results with original implementation on test data
- Verify identical outputs with improved performance
- Benchmark speed improvement

### Files Modified
- `R/utils.R` (add optimization function)
- `R/forecast.R` (lines 317-327)

---

## Step 4: Create Unified Stan Training Interface (30 min)

### Goal
Consolidate 80% duplicate code across Stan model training functions into single interface.

### Current Problem
`train_ardf()`, `train_vardf()`, and `train_gpdf()` contain nearly identical code for parameter handling and Stan sampling.

### Tasks
1. Create new `R/stan_training.R` with unified training function
2. Implement common parameter processing and Stan sampling logic
3. Add model-specific data preparation dispatch
4. Replace existing training functions with thin wrappers calling unified interface
5. Test each model type maintains correct behavior

### Testing
- Verify all three model types work correctly with new interface
- Test parameter passing and default handling
- Ensure consistent behavior across model types

### Files Modified
- `R/stan_training.R` (new file)
- `R/ardf.R`, `R/vardf.R`, `R/gpdf.R` (replace training functions)

---

## Step 5: Streamline Data Transformation Pipeline (30 min)

### Goal
Modularize complex 150+ line data transformation chain in `forecast.ffc_gam()` for better testability and clarity.

### Current Problem
Lines 370-520 contain deeply nested dplyr chains and complex pivoting operations that are difficult to debug and test.

### Tasks
1. Create new `R/data_processing.R` with modular transformation functions
2. Extract time structure handling into separate function
3. Isolate functional coefficient extraction logic
4. Create dedicated forecast basis generation function
5. Replace complex pipeline with clear function sequence
6. Add error handling to each component

### Testing
- Test each modular function independently
- Verify end-to-end pipeline produces identical results
- Check error handling in edge cases

### Files Modified
- `R/data_processing.R` (new file)
- `R/forecast.R` (lines 370-520 refactor)

---

## Step 6: Implement Statistics Optimization (30 min)

### Goal
Optimize repeated statistical computations in prediction summarization using vectorized operations.

### Current Problem
Multiple `apply(preds, 2, ...)` calls for computing statistics in lines 621-628 of `R/forecast.R` are inefficient for large prediction matrices.

### Tasks
1. Create vectorized statistics function in `R/stats_utils.R`
2. Add matrixStats dependency to DESCRIPTION for efficient column operations
3. Replace repetitive apply() calls with single vectorized computation
4. Support both robust and non-robust statistics options
5. Benchmark performance improvement

### Testing
- Verify statistical accuracy matches original implementation
- Test both robust and non-robust options
- Measure performance improvement on large prediction matrices

### Files Modified
- `DESCRIPTION` (add matrixStats import)
- `R/stats_utils.R` (new file)  
- `R/forecast.R` (lines 621-628)

---

## Quality Assurance Checklist

After each step:
- [ ] Run `Rscript -e "devtools::load_all(); devtools::check()"`
- [ ] No new ERRORs or WARNINGs introduced
- [ ] Use code-reviewer agent for major changes
- [ ] Test functionality with simple examples
- [ ] Commit changes with descriptive messages

## Final Validation

After completing all steps:
1. Run full test suite: `Rscript -e "devtools::test()"`
2. Check performance improvement on real data
3. Verify package passes `devtools::check()` without errors
4. Test complete prediction workflow end-to-end
5. Review all changes with code-reviewer agent

## Expected Outcomes

- **Performance**: 50-70% improvement in prediction speed
- **Maintainability**: 70% reduction in code duplication  
- **Clarity**: Modular, testable components
- **Scalability**: Better handling of large datasets
- **Quality**: Zero check ERRORs or WARNINGs