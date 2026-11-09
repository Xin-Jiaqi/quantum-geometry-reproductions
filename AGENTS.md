# Repository collaboration rules

- User preference (2026-11-09): do not enable GitHub Actions email notifications.
  Report relevant check outcomes in the conversation or project worklog instead.
- Research archives, figures and worklogs should not trigger unrelated model
  regressions. Preserve checks for model code, tests and numerical test fixtures;
  do not disable validation or security checks merely to suppress email.
- A skipped workflow is not a passed workflow. Keep local and remote validation
  claims separate, and batch routine result uploads when practical.
