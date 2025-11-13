# Contributing to ffc

This outlines how to propose a change to the `ffc` 📦. For a detailed
discussion on contributing to this and other tidyverse packages, please
see the [development contributing guide](https://rstd.io/tidy-contrib)
and our [code review principles](https://code-review.tidyverse.org/).

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the
documentation directly using the GitHub web interface, as long as the
changes are made in the *source* file. This generally means you’ll need
to edit [roxygen2
comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`,
not a `.Rd` file. You can find the `.R` file that generates the `.Rd` by
reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it’s a good idea to first file an
issue and make sure someone from the team agrees that it’s needed. If
you’ve found a bug, please file an issue that illustrates the bug with a
minimal [reprex](https://www.tidyverse.org/help/#reprex) (this will also
help you write a unit test, if needed). See our guide on [how to create
a great issue](https://code-review.tidyverse.org/issues/) for more
advice.

### Pull request process

- Fork the package and clone onto your computer. If you haven’t done
  this before, we recommend using
  `usethis::create_from_github("nicholasjclark/ffc", fork = TRUE)`.

- Install all development dependencies with
  [`devtools::install_dev_deps()`](https://devtools.r-lib.org/reference/install_deps.html),
  and then make sure the package passes R CMD check by running
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html).
  If R CMD check doesn’t pass cleanly, it’s a good idea to ask for help
  before continuing.

- Create a Git branch for your pull request (PR). We recommend using
  `usethis::pr_init("brief-description-of-change")`.

- Make your changes, commit to git, and then create a PR by running
  [`usethis::pr_push()`](https://usethis.r-lib.org/reference/pull-requests.html),
  and following the prompts in your browser. The title of your PR should
  briefly describe the change. The body of your PR should contain
  `Fixes #issue-number`.

- For user-facing changes, add a bullet to the top of `NEWS.md`
  (i.e. just below the first header). Follow the style described in
  <https://style.tidyverse.org/news.html>.

### Code style

- New code should follow the tidyverse [style
  guide](https://style.tidyverse.org). You can use the
  [styler](https://CRAN.R-project.org/package=styler) package to apply
  these styles, but please don’t restyle code that has nothing to do
  with your PR.

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
  [Markdown
  syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html),
  for documentation.

- We use [testthat](https://cran.r-project.org/package=testthat) for
  unit tests. Contributions with test cases included are easier to
  accept.

# Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in
our community a harassment-free experience for everyone, regardless of
age, body size, visible or invisible disability, ethnicity, sex
characteristics, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance,
race, caste, color, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open,
welcoming, diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

- Demonstrating empathy and kindness toward other people
- Being respectful of differing opinions, viewpoints, and experiences
- Giving and gracefully accepting constructive feedback
- Accepting responsibility and apologizing to those affected by our
  mistakes, and learning from the experience
- Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

- The use of sexualized language or imagery, and sexual attention or
  advances of any kind
- Trolling, insulting or derogatory comments, and personal or political
  attacks
- Public or private harassment
- Publishing others’ private information, such as a physical or email
  address, without their explicit permission
- Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our
standards of acceptable behavior and will take appropriate and fair
corrective action in response to any behavior that they deem
inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, and will
communicate reasons for moderation decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also
applies when an individual is officially representing the community in
public spaces. Examples of representing our community include using an
official e-mail address, posting via an official social media account,
or acting as an appointed representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may
be reported to the community leaders responsible for enforcement at
<codeofconduct@posit.co>. All complaints will be reviewed and
investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security
of the reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in
determining the consequences for any action they deem in violation of
this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior
deemed unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders,
providing clarity around the nature of the violation and an explanation
of why the behavior was inappropriate. A public apology may be
requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, for a specified period of
time. This includes avoiding interactions in community spaces as well as
external channels like social media. Violating these terms may lead to a
temporary or permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards,
including sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No
public or private interaction with the people involved, including
unsolicited interaction with those enforcing the Code of Conduct, is
allowed during this period. Violating these terms may lead to a
permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of
individuals.

**Consequence**: A permanent ban from any sort of public interaction
within the community.

## Attribution

This Code of Conduct is adapted from the [Contributor
Covenant](https://www.contributor-covenant.org), version 2.1, available
at
<https://www.contributor-covenant.org/version/2/1/code_of_conduct.html>.

Community Impact Guidelines were inspired by \[Mozilla’s code of conduct
enforcement ladder\]\[<https://github.com/mozilla/inclusion>\].

For answers to common questions about this code of conduct, see the FAQ
at <https://www.contributor-covenant.org/faq>. Translations are
available at <https://www.contributor-covenant.org/translations>.
