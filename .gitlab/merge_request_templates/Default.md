> When using the template, not all fields may be relevant to every merge request. Checklist items that are not relevant should be removed when creating the merge request. Example: if AI is not used, the checklist item to verify it in the DISCLAIMER.md is not needed.

# Merge Request Summary

- Description: `<One sentence description>`
- Related issues: `<Comma delimited list of related issues or N/A, for closing use Closes #<number>>`
- Reviewer: `<Set reviewer with @username>`
- Assignee: `<Set reviewer with @username>`
- Merger: `<Tag who should merge the request uppon approval (generally the reviwer)>`
- Tag: `<New tag version or N/A>`
  - Creator: `<Tag who should create the tag after merge or N/A>`
  - Release [y/n]: `<y/n>`
    - Creator: `<Tag who should create the release or N/A>`

# Detailed List of Changes

`<List of changes>`

## Considerations

`<Describe how this could impact future usage, development, etc.>`

## References

- AI Usage [y/n]: `<y/n>`
  - AI Model Used: `<Model name and version or N/A>`
  - Model usage description: `<Description of how model was used or N/A>`
- Other relevant citations, references, links, etc.:
  - `<List of references or N/A>`

<br><br>

---

> Documentation below is for the reviewer to check off as they perform the review.

---

# Reviewer Checklists

## Pre-Merge Checklist

> NOTE: Not all checks are applicable to every M.R. Reviews may use the --strikethrough-- anotation to mark irrelevant checklist items

- [ ] Fetched the merge request and perform relevant checks:
  > fetch m.r. with `git fetch upstream merge-requests/<merge request number>/head:<merge request source branch>`
  - [ ] Tests were added for code changes
  - [ ] Tests can run locally
- [ ] Pipeline checks pass
- [ ] All commits and changes reviewed for no presence of:
  - Personally identifiable information (PII)
  - Absolute file system paths
  - Internal server host names or IP addresses
  - Usernames/passwords
- [ ] Important changes added to CHANGELOG.md
- [ ] Relevant AI usage added to DISCLAIMER.md
- [ ] Important scientific references added to README.md
- [ ] Relevant [semantic version](https://semver.org/) change documented in relevant files (e.g. pyproject.toml, package.json, code.json)
- [ ] M.R. Approved

### Post-Merge Checklist

- [ ] Create relevant tags/releases if assigned
- [ ] Perform relevant deploy(s) if assigned
- [ ] Verify that pipelines run upon merge pass
  - Create issue for pipelines that do not pass upon merge
