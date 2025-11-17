# Claude Code Assistant Instructions for ffc Package

## Project Overview
The `ffc` (Functional ForeCasting) package performs functional regression using GAMs with time-varying functional coefficients, integrating mgcv and Stan for dynamic factor modeling.

## Directory Structure
- `R/` - All R source code and functions
- `inst/stan/` - Stan model files (ardf.stan, vardf.stan, gpdf.stan)
- `tests/testthat/` - Test files
- `man/` - Auto-generated documentation via roxygen2

## Object System
This package uses the S3 object system with classes:
- `ffc_gam` - Main model class
- `fts_ts` - Time series coefficient class
- Methods: `print`, `predict`, `forecast`, `autoplot`, `model.frame`

## Development Commands

### Testing
```bash
# ALWAYS load functions before testing
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_fts.R')"

# Run specific test by description (NOT by filter)
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test_fts.R', desc = 'specific test description')"

# Run all tests
Rscript -e "devtools::test()"
```

### Documentation
```bash
# Generate/update roxygen2 documentation
Rscript -e "devtools::document()"

# Build pkgdown site
Rscript -e "pkgdown::build_site()"
```

### Package Checking
```bash
# Check package for errors, warnings, and notes
Rscript -e "devtools::check()"
```

### Stan Model Management
```bash
# After modifying .stan files in inst/stan/:
Rscript -e "rstantools::rstan_config(); devtools::document()"

# Then reinstall package
R CMD INSTALL --preclean --no-multiarch .
```

## Code Standards

### R Code Style
- Tidyverse style guide
- Line length: 80 characters maximum
- Function names: snake_case, prefer shorter names
- Assignment: Use `<-` not `=`
- Indentation: 2 spaces

### Documentation
- Use roxygen2 tags for all functions
- Vignettes demonstrate in-depth use cases
- pkgdown site provides comprehensive documentation
- README provides high level overview
- roxygen2 examples demonstrate simpler use cases

### Comments
- Use sentence case, avoid ALL CAPS
- Add `# Reason:` comments for complex logic
- Keep professional and self-contained

## Testing Strategy

### Test Requirements
- **NO test errors or warnings allowed**
- Always `devtools::load_all()` before testing
- Use `desc` argument for selective test execution, NOT filter
- Add tests to appropriate existing test files

### testthat Guidelines
- Expectation functions do NOT have an `info` argument
- Add descriptive comments above expectations instead
- NO `try()`, `tryCatch()`, or `skip()` in tests
- NO wrapper functions to fix failing tests
- Update existing tests when logic changes

## Validation Patterns
```r
# Input validation
checkmate::assert_*()

# Error formatting
insight::format_error()

# Warning formatting  
insight::format_warning()

# One-time session warnings
if (!identical(Sys.getenv("TESTTHAT"), "true")) {
  rlang::warn("message", .frequency = "once")
}
```

### Message Standards
- Use `{.field parameter_name}` for parameter highlighting
- Include suggested solutions in errors
- Provide context for constraints

## Git Workflow

### Branch Management
- Descriptive branch names: `feature/brms-integration`
- Branch from main/master
- Keep branches focused
- Use `gh` CLI for GitHub interactions

### Commit Standards
```
Brief description (50 chars max)

- What was changed
- Why it was necessary
- Implementation notes
```
- Do NOT mention Claude or AI assistance in commits

### Pre-commit Workflow
1. Check git status
2. Stage appropriate files
3. Write clear commit messages

## Stan Development Notes

### Compilation Workflow
1. Edit .stan files in `inst/stan/`
2. Run `rstantools::rstan_config()` to update infrastructure
3. Run `devtools::document()` to regenerate interfaces
4. Reinstall package with `R CMD INSTALL --preclean --no-multiarch .`
5. Test compilation and sampling

### Stan Code Patterns
- Precompiled models accessed via `stanmodels$model_name`
- Use vectorized operations for efficiency
- Factor models use lower-triangular loading matrices
- Data preparation must match Stan data block exactly

## Common Tasks

### Adding Documentation
1. Update roxygen2 comments in R files
2. Run `devtools::document()`
3. Check generated man pages
4. Update vignettes if needed

### Streamlining Predictions
1. Modify predict methods in `R/predict.R`
2. Ensure S3 dispatch works correctly
3. Update tests in existing test files
4. Document any new parameters

### Adding Use Cases
1. Create vignettes for complex examples
2. Add simple examples to roxygen2 blocks
3. Update README for overview changes

## Important Constraints
- Zero tolerance for test failures or warnings
- Maintain tidyverse style (80 char lines)
- Always load functions before testing
- Use existing test files when possible
- Document all exported functions