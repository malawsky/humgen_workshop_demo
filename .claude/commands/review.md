Run an adversarial review of the work since the last commit.

Use the **reviewer** subagent. Show it `git diff` and `git diff --staged`, have
it list concrete problems (each with a file and line) following its brief, and
end with the single most important issue to fix first. Do not change any code
in this command — only review.
