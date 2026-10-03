# Agent skills for kmerseek

Instructions for any coding agent working in this repo (Claude Code, Codex, Cursor,
Copilot, or a person). Each folder holds one skill: a `SKILL.md` with a `name` and a
`description` saying when to use it, plus any reference files or scripts it needs. The
format follows the open Agent Skills layout, so tools that read `.agents/skills/` pick
them up directly; for any other tool, point it at the `SKILL.md` by path.

| Skill | Use it when |
|---|---|
| [`small-prs`](small-prs/SKILL.md) | Planning, splitting, describing, reviewing or answering review on a PR |
| [`rustacean-review`](rustacean-review/SKILL.md) | Reviewing Rust code, a diff or a PR, or before merging Rust |
| [`are-you-sure`](are-you-sure/SKILL.md) | Before calling any code, command, commit or PR done |
| [`writing-style`](writing-style/SKILL.md) | Writing any prose: docs, PR bodies, commit messages, doc comments |
| [`clear-figures`](clear-figures/SKILL.md) | Drawing any figure or diagram, including the figure a sequence test needs |
| [`publication-figures`](publication-figures/SKILL.md) | Writing matplotlib code (fonts, sizes, colours, export) |

These are copies of Olga's personal skills. When a rule changes, change it here in a PR,
so every agent and every contributor reads the same version.
